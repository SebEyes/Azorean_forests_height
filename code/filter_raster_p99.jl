#!/usr/bin/env julia
#
# filter_raster_p99.jl
#
# Clamp the merged, per-island vegetation-height rasters to their 99th
# percentile: for every island mosaic, compute P99 over all valid (non-NoData)
# pixels and write a copy in which every value strictly greater than P99 is
# replaced by NoData ("NA").
#
# Both passes read the rasters in horizontal row windows, so even the 1.4 GB
# São Miguel mosaic is never held in memory at once. The originals are kept
# untouched; filtered mosaics are written as <CODE>_vegetation_height_p99.tif.
#
# Example:
#   julia --threads=4 code/filter_raster_p99.jl \
#       --raster-dir data/vegetation_height \
#       --outdir data/vegetation_height
#   julia code/filter_raster_p99.jl --island COR --outdir data/vegetation_height/qc

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
Pkg.instantiate()

using ArchGDAL
using Printf

const DEFAULT_BINWIDTH    = 0.01   # m; histogram resolution used for P99
const DEFAULT_MAXHEIGHT   = 120.0  # m; larger values go to the overflow bin
const DEFAULT_BLOCK_ROWS  = 256
const DEFAULT_SUFFIX      = "_p99"
const NODATA_FALLBACK     = Float32(-9999)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

function parse_cli()
    opts = Dict{String,String}()
    flags = Set{String}()
    i = 1
    while i <= length(ARGS)
        a = ARGS[i]
        if startswith(a, "--")
            if a == "--help"
                push!(flags, "help")
            elseif occursin('=', a)
                k, v = split(a[3:end], '=', limit=2)
                opts[k] = v
            else
                k = a[3:end]
                i += 1
                i <= length(ARGS) || error("Missing value for --$k")
                opts[k] = ARGS[i]
            end
        end
        i += 1
    end
    return opts, flags
end

function help()
    println("""
    Usage: julia code/filter_raster_p99.jl [options]

    Options:
      --raster-dir DIR    Directory with merged island mosaics
                          (default: data/vegetation_height)
      --outdir DIR        Directory for the filtered mosaics and the summary CSV
                          (default: same as --raster-dir)
      --suffix STR        Suffix appended to the output filename (default: $(DEFAULT_SUFFIX))
      --island CODE       Process only one island (SMG SMA TER GRA SJO PIC FAI FLO COR)
      --bin-width M       Histogram resolution for the P99 estimate (default: $DEFAULT_BINWIDTH)
      --max-height M      Top of the histogram; larger values go to the overflow bin
                          (default: $DEFAULT_MAXHEIGHT)
      --block-rows N      Rows per read/write window (default: $DEFAULT_BLOCK_ROWS)
      --help              Show this help

    For each island the summary CSV (p99_filter_summary.csv) reports the P99
    threshold, the number of valid pixels and the number/percentage of pixels
    set to NoData.
    """)
end

function island_code(filename::AbstractString)
    m = match(r"^([A-Z]{3})_vegetation_height\.tif$", basename(filename))
    return m === nothing ? nothing : m.captures[1]
end

# ---------------------------------------------------------------------------
# Raster metadata and windowed scans
# ---------------------------------------------------------------------------

function band_info(path::AbstractString)
    ds = ArchGDAL.read(path)
    try
        band = ArchGDAL.getband(ds, 1)
        nd = ArchGDAL.getnodatavalue(band)
        return (
            width  = ArchGDAL.width(ds),
            height = ArchGDAL.height(ds),
            gt     = collect(ArchGDAL.getgeotransform(ds)),
            wkt    = ArchGDAL.getproj(ds),
            ndval  = nd === nothing ? NODATA_FALLBACK : Float32(nd),
        )
    finally
        ArchGDAL.destroy(ds)
    end
end

"""
    scan_histogram(path, info, binwidth, maxheight, blockrows)

First pass: accumulate a fine histogram of all valid pixel values. The last
bin is an overflow bin holding every value >= maxheight. Returns the histogram
and the number of NoData pixels.
"""
function scan_histogram(path::AbstractString, info, binwidth::Float64,
                        maxheight::Float64, blockrows::Int)
    nbins = ceil(Int, maxheight / binwidth)
    hist = zeros(Int, nbins + 1)          # +1: overflow bin
    nodata = 0

    ds = ArchGDAL.read(path)
    try
        for rowstart in 1:blockrows:info.height
            nrows = min(blockrows, info.height - rowstart + 1)
            a = ArchGDAL.read(ds, 1, 0, rowstart - 1, info.width, nrows)
            @inbounds for v in a
                if v == info.ndval || !isfinite(v)
                    nodata += 1
                elseif v >= maxheight
                    hist[end] += 1
                else
                    bin = clamp(floor(Int, Float64(v) / binwidth) + 1, 1, nbins)
                    hist[bin] += 1
                end
            end
        end
    finally
        ArchGDAL.destroy(ds)
    end
    return hist, nodata
end

"Linear-interpolated quantile of a histogram whose first bin starts at zero."
function hist_quantile(h::Vector{Int}, p::Float64, binwidth::Float64)
    n = sum(h)
    n == 0 && return NaN
    target = p * n
    running = 0
    for (i, c) in enumerate(h)
        running += c
        if running >= target
            prev = running - c
            frac = c == 0 ? 0.0 : (target - prev) / c
            return (i - 1 + frac) * binwidth
        end
    end
    return length(h) * binwidth
end

"""
    write_filtered(path, outpath, threshold, info, blockrows)

Second pass: copy the raster, replacing every valid value strictly above
`threshold` with the NoData value. Returns the number of masked pixels.
"""
function write_filtered(path::AbstractString, outpath::AbstractString,
                        threshold::Float64, info, blockrows::Int)
    masked = 0
    ds = ArchGDAL.read(path)
    mkpath(dirname(abspath(outpath)))
    options = ["COMPRESS=DEFLATE", "TILED=YES", "BLOCKXSIZE=256",
               "BLOCKYSIZE=256", "PREDICTOR=3"]

    try
        ArchGDAL.create(outpath; driver=ArchGDAL.getdriver("GTiff"),
                        width=info.width, height=info.height, nbands=1,
                        dtype=Float32, options=options) do ods
            ArchGDAL.setgeotransform!(ods, info.gt)
            ArchGDAL.setproj!(ods, info.wkt)
            ArchGDAL.setnodatavalue!(ArchGDAL.getband(ods, 1), info.ndval)

            for rowstart in 1:blockrows:info.height
                nrows = min(blockrows, info.height - rowstart + 1)
                a = ArchGDAL.read(ds, 1, 0, rowstart - 1, info.width, nrows)
                out = copy(a)
                @inbounds for i in eachindex(a)
                    v = a[i]
                    if v != info.ndval && isfinite(v) && v > threshold
                        out[i] = info.ndval
                        masked += 1
                    end
                end
                ArchGDAL.write!(ods, out, 1, 0, rowstart - 1, info.width, nrows)
            end
        end
    finally
        ArchGDAL.destroy(ds)
    end
    return masked
end

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

fmt(x, digits=4) = isnan(x) ? "" : @sprintf("%.*f", digits, x)

function write_csv(path::AbstractString, header::Vector{String},
                   rows::Vector{Vector{String}})
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join(row, ","))
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    opts, flags = parse_cli()
    "help" in flags && (help(); return)

    rasterdir = get(opts, "raster-dir", joinpath("data", "vegetation_height"))
    outdir    = get(opts, "outdir", rasterdir)
    suffix    = get(opts, "suffix", DEFAULT_SUFFIX)
    binwidth  = parse(Float64, get(opts, "bin-width", string(DEFAULT_BINWIDTH)))
    maxheight = parse(Float64, get(opts, "max-height", string(DEFAULT_MAXHEIGHT)))
    blockrows = parse(Int, get(opts, "block-rows", string(DEFAULT_BLOCK_ROWS)))
    only      = haskey(opts, "island") ? uppercase(strip(opts["island"])) : nothing

    binwidth > 0 && maxheight > 0 && blockrows > 0 || error("bin-width, max-height and block-rows must be positive")
    isdir(rasterdir) || error("Raster directory not found: $rasterdir")
    mkpath(outdir)

    paths = [p for p in sort(readdir(rasterdir; join=true))
             if (c = island_code(p)) !== nothing && (only === nothing || c == only)]
    isempty(paths) && error("No merged island rasters found in $rasterdir")

    println("P99 filtering of merged vegetation-height rasters")
    println("  input     : $(abspath(rasterdir))")
    println("  output    : $(abspath(outdir))  (suffix $suffix)")
    println("  P99 basis : all valid pixels, $(binwidth) m histogram bins")
    println("  rasters   : $(length(paths))\n")
    println("Island   P99 (m)   Valid px     Masked px   Masked %")

    header = ["island", "input_raster", "output_raster", "p99_m",
              "n_valid_pixels", "n_nodata_pixels", "n_masked_pixels", "pct_masked"]
    rows = Vector{Vector{String}}()

    for p in paths
        code = island_code(p)
        info = band_info(p)
        hist, nodata = scan_histogram(p, info, binwidth, maxheight, blockrows)
        nvalid = sum(hist)
        p99 = hist_quantile(hist, 0.99, binwidth)

        outname = code * "_vegetation_height" * suffix * ".tif"
        outpath = joinpath(outdir, outname)
        masked = write_filtered(p, outpath, p99, info, blockrows)
        pct = nvalid == 0 ? 0.0 : 100.0 * masked / nvalid

        @printf("%-6s  %8.3f  %11d  %11d  %8.4f%%\n",
                code, p99, nvalid, masked, pct)
        push!(rows, [code, basename(p), outname, fmt(p99), string(nvalid),
                     string(nodata), string(masked), fmt(pct, 6)])
    end

    summary = joinpath(outdir, "p99_filter_summary.csv")
    write_csv(summary, header, rows)
    println("\nSummary written to $summary")
end

(abspath(PROGRAM_FILE) == @__FILE__) && main()
