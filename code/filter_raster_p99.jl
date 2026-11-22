#!/usr/bin/env julia
#
# filter_raster_p99.jl
#
# Clamp the merged, per-island vegetation-height rasters to a chosen upper
# percentile (default P99): for every island mosaic, compute the percentile
# over all valid (non-NoData) pixels and write a copy in which every value
# strictly greater than the threshold is replaced by NoData ("NA").
#
# Both passes read the rasters in horizontal row windows, so even the 1.4 GB
# São Miguel mosaic is never held in memory at once. The originals are kept
# untouched; filtered mosaics are written as
# <CODE>_vegetation_height_<percentile>.tif (e.g. _p99, _p999). With
# --interpolate, the removed pixels are filled by bilinear interpolation from
# the surrounding valid pixels; original NoData pixels are never filled.
#
# Example:
#   julia code/filter_raster_p99.jl --percentile 0.999 --interpolate \
#       --raster-dir data/vegetation_height --outdir data/vegetation_height
#   julia code/filter_raster_p99.jl --island COR --outdir data/vegetation_height/qc

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
Pkg.instantiate()

using ArchGDAL
using Printf

const DEFAULT_PERCENTILE  = 0.99
const DEFAULT_BINWIDTH    = 0.01   # m; histogram resolution used for the percentile
const DEFAULT_MAXHEIGHT   = 120.0  # m; larger values go to the overflow bin
const DEFAULT_BLOCK_ROWS  = 256
const NODATA_FALLBACK     = Float32(-9999)

"Short label for a percentile, e.g. 0.99 -> p99, 0.999 -> p999."
percent_label(p::Float64) = "p" * replace(@sprintf("%g", 100p), "." => "")

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
            if a in ("--help", "--interpolate")
                push!(flags, a[3:end])
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
      --percentile P      Upper percentile, e.g. 0.99 or 0.999 (default: $DEFAULT_PERCENTILE)
      --interpolate       Fill pixels removed as outliers by bilinear interpolation
                          (original NoData pixels are left as NoData)
      --suffix STR        Suffix for output filenames (default: _<percentile>, e.g. _p999)
      --summary FILE      Summary CSV filename (default: <percentile>_filter_summary.csv)
      --island CODE       Process only one island (SMG SMA TER GRA SJO PIC FAI FLO COR)
      --bin-width M       Histogram resolution for the percentile (default: $DEFAULT_BINWIDTH)
      --max-height M      Top of the histogram; larger values go to the overflow bin
                          (default: $DEFAULT_MAXHEIGHT)
      --block-rows N      Rows per read/write window (default: $DEFAULT_BLOCK_ROWS)
      --help              Show this help

    For each island the summary CSV reports the threshold, the number of valid
    pixels and the number/percentage of pixels set to NoData.
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

# ---------------------------------------------------------------------------
# Bilinear filling of the gaps created by the outlier filter
# ---------------------------------------------------------------------------

"""
    horizontal_brackets!(leftidx, rightidx, a, valid)

For every cell of the current window, record the nearest valid column to its
left and right (0 when none). `a` is laid out as [col, row].
"""
function horizontal_brackets!(leftidx, rightidx, a, valid)
    ncols, nrows = size(a)
    @inbounds for r in 1:nrows
        last = 0
        for c in 1:ncols
            leftidx[c, r] = last
            valid[c, r] && (last = c)
        end
        last = ncols + 1
        for c in ncols:-1:1
            rightidx[c, r] = last == ncols + 1 ? 0 : last
            valid[c, r] && (last = c)
        end
    end
    return nothing
end

"Linearly interpolate valid values around `c` within row `r`; NaN if impossible."
function row_value(a, valid, leftidx, rightidx, c, r)
    @inbounds begin
        valid[c, r] && return Float64(a[c, r])
        li, ri = leftidx[c, r], rightidx[c, r]
        (li == 0 || ri == 0) && return NaN
        dl, dr = c - li, ri - c
        return (dr * Float64(a[li, r]) + dl * Float64(a[ri, r])) / (dl + dr)
    end
end

"""
    bilinear_value(a, valid, leftidx, rightidx, c, r)

Bilinear fill: interpolate horizontally in the nearest rows above and below
that have valid values around column `c`, then interpolate vertically between
those two estimates. Falls back to the nearest available row; NaN if none.
"""
function bilinear_value(a, valid, leftidx, rightidx, c, r)
    up = NaN
    dup = 0
    @inbounds for rr in (r - 1):-1:1
        if valid[c, rr] || (leftidx[c, rr] != 0 && rightidx[c, rr] != 0)
            up = row_value(a, valid, leftidx, rightidx, c, rr)
            if !isnan(up)
                dup = r - rr
                break
            end
        end
    end

    down = NaN
    ddown = 0
    @inbounds for rr in (r + 1):size(a, 2)
        if valid[c, rr] || (leftidx[c, rr] != 0 && rightidx[c, rr] != 0)
            down = row_value(a, valid, leftidx, rightidx, c, rr)
            if !isnan(down)
                ddown = rr - r
                break
            end
        end
    end

    if !isnan(up) && !isnan(down)
        return (ddown * up + dup * down) / (dup + ddown)
    elseif !isnan(up)
        return up
    elseif !isnan(down)
        return down
    else
        return NaN
    end
end

"""
    write_filtered(path, outpath, threshold, info, blockrows;
                   interpolate=false, halo=64)

Second pass: copy the raster. Values above `threshold` are either set to
NoData or, when `interpolate=true`, replaced by a bilinear estimate from
surrounding valid pixels. Original NoData pixels are never filled.
Returns (n_masked, n_filled, n_unfilled).
"""
function write_filtered(path::AbstractString, outpath::AbstractString,
                        threshold::Float64, info, blockrows::Int;
                        interpolate::Bool=false, halo::Int=64)
    masked = filled = unfilled = 0
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
                r0 = max(1, rowstart - halo)
                r1 = min(info.height, rowstart + nrows - 1 + halo)
                a = ArchGDAL.read(ds, 1, 0, r0 - 1, info.width, r1 - r0 + 1)
                firstrow = rowstart - r0 + 1
                out = similar(a, Float32, info.width, nrows)

                if interpolate
                    valid = similar(a, Bool)
                    @inbounds for i in eachindex(a)
                        v = a[i]
                        valid[i] = v != info.ndval && isfinite(v) && v <= threshold
                    end
                    leftidx = zeros(Int32, size(a))
                    rightidx = similar(leftidx)
                    horizontal_brackets!(leftidx, rightidx, a, valid)
                end

                @inbounds for j in 1:nrows
                    wr = firstrow + j - 1
                    for c in 1:info.width
                        v = a[c, wr]
                        if v == info.ndval || !isfinite(v)
                            out[c, j] = info.ndval
                        elseif v <= threshold
                            out[c, j] = v
                        else
                            masked += 1
                            if interpolate
                                est = bilinear_value(a, valid, leftidx, rightidx, c, wr)
                                if isnan(est)
                                    out[c, j] = info.ndval
                                    unfilled += 1
                                else
                                    out[c, j] = Float32(est)
                                    filled += 1
                                end
                            else
                                out[c, j] = info.ndval
                            end
                        end
                    end
                end
                ArchGDAL.write!(ods, out, 1, 0, rowstart - 1, info.width, nrows)
            end
        end
    finally
        ArchGDAL.destroy(ds)
    end
    return masked, filled, unfilled
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

    rasterdir  = get(opts, "raster-dir", joinpath("data", "vegetation_height"))
    outdir     = get(opts, "outdir", rasterdir)
    percentile = parse(Float64, get(opts, "percentile", string(DEFAULT_PERCENTILE)))
    label      = percent_label(percentile)
    suffix     = get(opts, "suffix", "_" * label)
    summaryname = get(opts, "summary", label * "_filter_summary.csv")
    binwidth   = parse(Float64, get(opts, "bin-width", string(DEFAULT_BINWIDTH)))
    maxheight  = parse(Float64, get(opts, "max-height", string(DEFAULT_MAXHEIGHT)))
    blockrows  = parse(Int, get(opts, "block-rows", string(DEFAULT_BLOCK_ROWS)))
    only       = haskey(opts, "island") ? uppercase(strip(opts["island"])) : nothing
    interpolate = "interpolate" in flags

    0 < percentile < 1 || error("percentile must be strictly between 0 and 1")
    binwidth > 0 && maxheight > 0 && blockrows > 0 || error("bin-width, max-height and block-rows must be positive")
    isdir(rasterdir) || error("Raster directory not found: $rasterdir")
    mkpath(outdir)

    paths = [p for p in sort(readdir(rasterdir; join=true))
             if (c = island_code(p)) !== nothing && (only === nothing || c == only)]
    isempty(paths) && error("No merged island rasters found in $rasterdir")

    println("Upper-percentile filtering of merged vegetation-height rasters")
    println("  input        : $(abspath(rasterdir))")
    println("  output       : $(abspath(outdir))  (suffix $suffix)")
    println("  threshold    : $(label) of all valid pixels, $(binwidth) m histogram bins")
    println("  outlier NAs  : $(interpolate ? "bilinearly interpolated" : "left as NoData")")
    println("  rasters      : $(length(paths))\n")
    println("Island   $(rpad(label * " (m)", 9))Valid px     Masked px   Masked %")

    header = ["island", "input_raster", "output_raster", label * "_m",
              "n_valid_pixels", "n_nodata_pixels", "n_masked_pixels", "pct_masked",
              "n_interpolated", "n_unfilled", "interpolated"]
    rows = Vector{Vector{String}}()

    for p in paths
        code = island_code(p)
        info = band_info(p)
        hist, nodata = scan_histogram(p, info, binwidth, maxheight, blockrows)
        nvalid = sum(hist)
        threshold = hist_quantile(hist, percentile, binwidth)

        outname = code * "_vegetation_height" * suffix * ".tif"
        outpath = joinpath(outdir, outname)
        masked, filled, unfilled = write_filtered(p, outpath, threshold, info,
                                                  blockrows; interpolate=interpolate)
        pct = nvalid == 0 ? 0.0 : 100.0 * masked / nvalid

        @printf("%-6s  %8.3f  %11d  %11d  %8.4f%%\n",
                code, threshold, nvalid, masked, pct)
        push!(rows, [code, basename(p), outname, fmt(threshold), string(nvalid),
                     string(nodata), string(masked), fmt(pct, 6),
                     string(filled), string(unfilled), string(interpolate)])
    end

    summary = joinpath(outdir, summaryname)
    write_csv(summary, header, rows)
    println("\nSummary written to $summary")
end

(abspath(PROGRAM_FILE) == @__FILE__) && main()
