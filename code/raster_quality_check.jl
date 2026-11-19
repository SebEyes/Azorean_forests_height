#!/usr/bin/env julia
#
# raster_quality_check.jl
#
# Check the distribution of the merged, per-island vegetation-height rasters
# and locate outliers, separately for every island.
#
# The mosaics produced by compute_archipelago.jl (one GeoTIFF per island) are
# scanned in row blocks, so even the 1.4 GB São Miguel raster is never loaded
# into memory at once.  For each island the script accumulates a fine histogram
# of all valid pixels, robust summary statistics (quantiles, Tukey fences),
# counts above physically meaningful height thresholds, the location of the
# most extreme pixels, and a per-row-block outlier summary for mapping.
#
# Outputs (default: data/vegetation_height/qc):
#   height_distribution_by_island.csv  summary statistics, one row per island
#   height_histogram.csv               histogram counts, one row per island/bin
#   height_outlier_blocks.csv          per-row-block outlier counts and maxima
#   height_extreme_pixels.csv          pixel locations above --extreme-threshold
#   height_distribution_report.html    self-contained report with SVG histograms
#
# Example:
#   julia --threads=6 code/raster_quality_check.jl \
#       --raster-dir data/vegetation_height \
#       --outdir data/vegetation_height/qc
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
Pkg.instantiate()

using ArchGDAL
using Printf

const ISLAND_NAMES = Dict(
    "SMG" => "São Miguel",
    "SMA" => "Santa Maria",
    "TER" => "Terceira",
    "GRA" => "Graciosa",
    "SJO" => "São Jorge",
    "PIC" => "Pico",
    "FAI" => "Faial",
    "FLO" => "Flores",
    "COR" => "Corvo",
)

const DEFAULT_BINWIDTH       = 0.25   # m
const DEFAULT_MAXHEIGHT      = 120.0  # m; larger values go to the overflow bin
const DEFAULT_EXTREME        = 40.0   # m; pixel locations written to CSV
const DEFAULT_MAX_EXTREME    = 1_000_000
const DEFAULT_BLOCK_ROWS     = 256
const THRESHOLDS             = [20.0, 25.0, 30.0, 35.0, 40.0, 50.0]
const FOREST_MIN_HEIGHT      = 2.0    # m, class-5 ("high vegetation") threshold
const NODATA_FALLBACK        = Float32(-9999)

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
    Usage: julia --threads=N code/raster_quality_check.jl [options]

    Options:
      --raster-dir DIR         Directory containing the merged island mosaics
                               (default: data/vegetation_height)
      --outdir DIR             Directory for CSV tables and HTML report
                               (default: data/vegetation_height/qc)
      --island CODE            Restrict the analysis to one island
                               (SMG SMA TER GRA SJO PIC FAI FLO COR, or a name)
      --bin-width M            Histogram bin width (default $DEFAULT_BINWIDTH)
      --max-height M           Top of the histogram, in m (default $DEFAULT_MAXHEIGHT)
      --extreme-threshold M    Write pixel locations above this height to CSV
                               (default $DEFAULT_EXTREME)
      --max-extreme N          Maximum number of extreme pixels written (default $DEFAULT_MAX_EXTREME)
      --block-rows N           Rows per read window (default $DEFAULT_BLOCK_ROWS)
      --help                   Show this help

    The number of windows processed in parallel is set by JULIA_NUM_THREADS or
    `julia --threads=N`.
    """)
end

function resolve_island(s::AbstractString)
    u = uppercase(strip(s))
    haskey(ISLAND_NAMES, u) && return u
    for (code, name) in ISLAND_NAMES
        lowercase(s) == lowercase(name) && return code
    end
    error("Unknown island '$s'. Use one of: $(join(sort(collect(keys(ISLAND_NAMES))), ", "))")
end

island_code(filename::AbstractString) =
    (m = match(r"^([A-Z]{3})_vegetation_height\.tif$", basename(filename))) === nothing ?
    "UNK" : m.captures[1]

# ---------------------------------------------------------------------------
# Windowed scan of the merged rasters
# ---------------------------------------------------------------------------

struct ExtremePixel
    height::Float32
    col::Int
    row::Int
    x::Float64
    y::Float64
end

struct BlockResult
    island::String
    raster::String
    rowstart::Int              # 1-based row of the first line in this window
    nrows::Int
    ytop::Float64
    ybottom::Float64
    nodata::Int
    negative::Int
    zero::Int
    hist::Vector{Int}          # bin i covers [(i-1)*w, i*w)
    overflow::Int              # values >= maxheight
    maxheight::Float32
    above::Vector{Int}         # counts above each fixed threshold
    extreme::Vector{ExtremePixel}
end

"""
    process_block(raster, island, rowstart, nrows, binwidth, nbins,
                  maxheight, extremethresh)

Open one merged raster and scan one horizontal window of `nrows` lines.
The window is read as a whole but never larger than a few tens of MB.
"""
function process_block(raster::AbstractString, island::String, rowstart::Int,
                       nrows::Int, binwidth::Float64, nbins::Int,
                       maxheight::Float64, extremethresh::Float64)
    ds = ArchGDAL.read(raster)
    hist = zeros(Int, nbins)
    extreme = ExtremePixel[]
    above = zeros(Int, length(THRESHOLDS))
    nodata = negative = zero = overflow = 0
    maxv = Float32(-Inf)

    try
        band = ArchGDAL.getband(ds, 1)
        nd = ArchGDAL.getnodatavalue(band)
        ndval = nd === nothing ? NODATA_FALLBACK : Float32(nd)
        gt = ArchGDAL.getgeotransform(ds)
        a = ArchGDAL.read(ds, 1, 0, rowstart - 1, ArchGDAL.width(ds), nrows)
        nx, ny = size(a)
        rowend = rowstart + ny - 1
        ytop = gt[4] + (rowstart - 1) * gt[6]
        ybottom = gt[4] + rowend * gt[6]

        @inbounds for row in 1:ny, col in 1:nx
            v = a[col, row]
            if v == ndval
                nodata += 1
                continue
            end
            v > maxv && (maxv = v)
            for k in eachindex(THRESHOLDS)
                v > THRESHOLDS[k] && (above[k] += 1)
            end
            if v < 0
                negative += 1
            elseif v == 0
                zero += 1
                hist[1] += 1
            elseif v >= maxheight
                overflow += 1
            else
                bin = floor(Int, Float64(v) / binwidth) + 1
                hist[bin] += 1
            end
            if v > extremethresh
                grow = rowstart - 1 + row
                x = gt[1] + (col - 0.5) * gt[2] + (grow - 0.5) * gt[3]
                y = gt[4] + (col - 0.5) * gt[5] + (grow - 0.5) * gt[6]
                push!(extreme, ExtremePixel(v, col, grow, x, y))
            end
        end
        return BlockResult(island, basename(raster), rowstart, ny, ytop, ybottom,
                           nodata, negative, zero, hist, overflow,
                           isfinite(maxv) ? maxv : Float32(NaN), above, extreme)
    finally
        ArchGDAL.destroy(ds)
    end
end

# ---------------------------------------------------------------------------
# Statistics derived from histograms
# ---------------------------------------------------------------------------

"Add the overflow bin to a histogram so that it covers the whole value range."
full_hist(h::Vector{Int}, overflow::Int) = [h; overflow]

"Quantile from a histogram whose first bin starts at `lower`."
function hist_quantile(h::Vector{Int}, p::Float64, binwidth::Float64, lower::Float64)
    n = sum(h)
    n == 0 && return NaN
    target = p * n
    running = 0
    for (i, c) in enumerate(h)
        running += c
        if running >= target
            prev = running - c
            frac = c == 0 ? 0.0 : (target - prev) / c
            return lower + (i - 1 + frac) * binwidth
        end
    end
    return lower + length(h) * binwidth
end

function hist_mean(h::Vector{Int}, binwidth::Float64, lower::Float64)
    n = sum(h)
    n == 0 && return NaN
    s = 0.0
    @inbounds for (i, c) in enumerate(h)
        s += (lower + (i - 0.5) * binwidth) * c
    end
    return s / n
end

function hist_std(h::Vector{Int}, binwidth::Float64, lower::Float64)
    n = sum(h)
    n == 0 && return NaN
    μ = hist_mean(h, binwidth, lower)
    s = 0.0
    @inbounds for (i, c) in enumerate(h)
        d = lower + (i - 0.5) * binwidth - μ
        s += d * d * c
    end
    return sqrt(s / n)
end

"Count of pixels strictly above `t`, using bins whose lower edge is >= t."
function hist_count_above(h::Vector{Int}, t::Float64, binwidth::Float64, lower::Float64=0.0)
    length(h) == 0 && return 0
    k = floor(Int, (t - lower) / binwidth) + 1
    k >= length(h) && return h[end]
    k <= 0 && return sum(h)
    return sum(@view h[k:end])
end

function summarize_island(code::String, blocks::Vector{BlockResult}, binwidth::Float64,
                          maxheight::Float64)
    hist = zeros(Int, length(first(blocks).hist))
    overflow = nodata = negative = zero = 0
    above = zeros(Int, length(THRESHOLDS))
    maxv = Float32(-Inf)
    for b in blocks
        hist .+= b.hist
        overflow += b.overflow
        nodata += b.nodata
        negative += b.negative
        zero += b.zero
        above .+= b.above
        isfinite(b.maxheight) && b.maxheight > maxv && (maxv = b.maxheight)
    end
    h = full_hist(hist, overflow)
    n = sum(h)

    q = [hist_quantile(h, p, binwidth, 0.0) for p in
         (0.01, 0.05, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.999, 0.9999)]
    iqr = q[7] - q[3]
    fence15 = q[7] + 1.5 * iqr
    fence30 = q[7] + 3.0 * iqr

    # "Forest" = high-vegetation pixels (> 2 m), i.e. class 5 in the source data.
    firstforest = floor(Int, FOREST_MIN_HEIGHT / binwidth) + 1
    fh = firstforest > length(h) ? Int[] : h[firstforest:end]
    forestn = sum(fh)
    forestq = forestn == 0 ? fill(NaN, 4) :
        [hist_quantile(fh, p, binwidth, (firstforest - 1) * binwidth) for p in
         (0.5, 0.90, 0.99, 0.999)]
    forestfence = forestn == 0 ? NaN : forestq[2] + 1.5 * (forestq[2] - forestq[1])

    return (
        code = code, name = get(ISLAND_NAMES, code, code), nblocks = length(blocks),
        n_valid = n, n_nodata = nodata, n_negative = negative, n_zero = zero,
        pct_zero = n == 0 ? NaN : 100 * zero / n,
        mean = hist_mean(h, binwidth, 0.0), std = hist_std(h, binwidth, 0.0),
        p01 = q[1], p05 = q[2], p25 = q[3], p50 = q[4], p75 = q[5],
        p90 = q[6], p95 = q[7], p99 = q[8], p999 = q[9], p9999 = q[10],
        max = isfinite(maxv) ? Float64(maxv) : NaN,
        fence15 = fence15, fence30 = fence30,
        n_gt_fence15 = hist_count_above(h, fence15, binwidth),
        n_gt_fence30 = hist_count_above(h, fence30, binwidth),
        above = above,
        forest_n = forestn,
        forest_mean = forestn == 0 ? NaN : hist_mean(fh, binwidth, (firstforest - 1) * binwidth),
        forest_p50 = forestq[1], forest_p90 = forestq[2],
        forest_p99 = forestq[3], forest_p999 = forestq[4],
        forest_fence15 = forestfence,
        forest_n_gt_fence15 = forestn == 0 ? 0 :
            hist_count_above(fh, forestfence, binwidth, (firstforest - 1) * binwidth),
        hist = h,
    )
end

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

fmt(x, digits=4) = isnan(x) ? "" : @sprintf("%.*f", digits, x)

"Write a CSV file from a header and rows of already-formatted strings."
function write_csv(path::AbstractString, header::Vector{String}, rows::Vector{Vector{String}})
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join(row, ","))
        end
    end
    return nothing
end

function summary_rows(stats)
    header = ["island", "island_name", "n_blocks", "n_valid_pixels", "n_nodata_pixels",
              "n_negative_pixels", "n_zero_pixels", "pct_zero", "mean_m", "std_m",
              "p01_m", "p05_m", "p25_m", "p50_m", "p75_m", "p90_m", "p95_m",
              "p99_m", "p999_m", "p9999_m", "max_m",
              "tukey_fence_1.5iqr_m", "tukey_fence_3iqr_m", "n_gt_tukey_1.5iqr",
              "n_gt_tukey_3iqr", "n_gt_20m", "n_gt_25m", "n_gt_30m", "n_gt_35m",
              "n_gt_40m", "n_gt_50m", "forest_n_gt_2m", "forest_mean_m",
              "forest_p50_m", "forest_p90_m", "forest_p99_m", "forest_p999_m",
              "forest_tukey_fence_m", "forest_n_gt_tukey"]
    rows = Vector{Vector{String}}()
    for s in stats
        push!(rows, [
            s.code, s.name, string(s.nblocks), string(s.n_valid), string(s.n_nodata),
            string(s.n_negative), string(s.n_zero), fmt(s.pct_zero, 3),
            fmt(s.mean), fmt(s.std), fmt(s.p01), fmt(s.p05), fmt(s.p25),
            fmt(s.p50), fmt(s.p75), fmt(s.p90), fmt(s.p95), fmt(s.p99),
            fmt(s.p999), fmt(s.p9999), fmt(s.max), fmt(s.fence15), fmt(s.fence30),
            string(s.n_gt_fence15), string(s.n_gt_fence30),
            [string(s.above[k]) for k in eachindex(THRESHOLDS)]...,
            string(s.forest_n), fmt(s.forest_mean), fmt(s.forest_p50),
            fmt(s.forest_p90), fmt(s.forest_p99), fmt(s.forest_p999),
            fmt(s.forest_fence15), string(s.forest_n_gt_fence15),
        ])
    end
    return header, rows
end

function histogram_rows(stats, binwidth)
    header = ["island", "bin_lower_m", "bin_upper_m", "count", "cumulative", "is_overflow"]
    rows = Vector{Vector{String}}()
    for s in stats
        cum = 0
        for (i, c) in enumerate(s.hist)
            lower = (i - 1) * binwidth
            cum += c
            push!(rows, [s.code, fmt(lower, 4), fmt(lower + binwidth, 4),
                         string(c), string(cum), "0"])
        end
        lower = (length(s.hist) - 1) * binwidth
        cum += s.hist[end]
        push!(rows, [s.code, fmt(lower, 4), "Inf", string(s.hist[end]), string(cum), "1"])
    end
    return header, rows
end

function block_rows(blocks)
    header = ["island", "raster", "row_start", "row_end", "y_top", "y_bottom",
              "n_nodata", "n_valid", "max_m",
              "n_gt_20m", "n_gt_25m", "n_gt_30m", "n_gt_35m", "n_gt_40m", "n_gt_50m"]
    rows = Vector{Vector{String}}()
    for b in sort(blocks; by=b -> (b.island, b.rowstart))
        n = sum(b.hist) + b.overflow
        push!(rows, [b.island, b.raster, string(b.rowstart), string(b.rowstart + b.nrows - 1),
                     fmt(b.ytop, 2), fmt(b.ybottom, 2), string(b.nodata), string(n),
                     fmt(Float64(b.maxheight), 3), [string(x) for x in b.above]...])
    end
    return header, rows
end

function extreme_rows(blocks, extremethresh, maxextreme)
    header = ["island", "raster", "col", "row", "x_m", "y_m", "height_m"]
    all = Tuple{String,String,ExtremePixel,Float32}[]
    for b in blocks, p in b.extreme
        push!(all, (b.island, b.raster, p, p.height))
    end
    sort!(all; by=x -> -x[4])
    truncated = length(all) > maxextreme
    truncated && resize!(all, maxextreme)
    rows = Vector{Vector{String}}()
    for (island, raster, p, _) in all
        push!(rows, [island, raster, string(p.col), string(p.row),
                     fmt(p.x, 2), fmt(p.y, 2), fmt(Float64(p.height), 3)])
    end
    return header, rows, length(all), truncated
end

# ---------------------------------------------------------------------------
# HTML / SVG report
# ---------------------------------------------------------------------------

"SVG histogram on a log10(count+1) scale, with reference lines."
function histogram_svg(s, binwidth, plotmax)
    width, height = 900.0, 320.0
    left, right, top, bottom = 60.0, 20.0, 25.0, 45.0
    pw, ph = width - left - right, height - top - bottom
    nbins = min(length(s.hist), ceil(Int, plotmax / binwidth))
    counts = s.hist[1:nbins]
    maxlog = log10(maximum(counts) + 1)
    parts = String[]

    for k in 0:min(6, floor(Int, maxlog))
        y = top + ph - (k / max(1, maxlog)) * ph
        push!(parts, "<line x1='$left' y1='$(fmt(y,1))' x2='$(left+pw)' y2='$(fmt(y,1))' stroke='#e5e5e5'/>")
        push!(parts, "<text x='$(left-8)' y='$(fmt(y+4,1))' text-anchor='end' font-size='11' fill='#666'>10^$k</text>")
    end

    bw = pw / nbins
    for (i, c) in enumerate(counts)
        h = (log10(c + 1) / max(1, maxlog)) * ph
        x = left + (i - 1) * bw
        y = top + ph - h
        push!(parts, "<rect x='$(fmt(x,2))' y='$(fmt(y,2))' width='$(fmt(bw+0.3,2))' height='$(fmt(h,2))' fill='#4c78a8' stroke='none'><title>$(fmt((i-1)*binwidth,2))–$(fmt(i*binwidth,2)) m: $c px</title></rect>")
    end

    push!(parts, "<line x1='$left' y1='$(top+ph)' x2='$(left+pw)' y2='$(top+ph)' stroke='#333'/>")
    tickstep = plotmax <= 30 ? 5.0 : 10.0
    t = 0.0
    while t <= plotmax
        x = left + (t / plotmax) * pw
        push!(parts, "<line x1='$(fmt(x,1))' y1='$(top+ph)' x2='$(fmt(x,1))' y2='$(top+ph+5)' stroke='#333'/>")
        push!(parts, "<text x='$(fmt(x,1))' y='$(top+ph+20)' text-anchor='middle' font-size='11' fill='#333'>$t</text>")
        t += tickstep
    end
    push!(parts, "<text x='$(left+pw/2)' y='$(height-5)' text-anchor='middle' font-size='12' fill='#333'>height (m)</text>")

    refs = [("Tukey 1.5×IQR = $(fmt(s.fence15,1)) m", s.fence15, "#e45756"),
            ("30 m", 30.0, "#f58518"), ("40 m", 40.0, "#b279a2")]
    for (label, xval, color) in refs
        (0 < xval <= plotmax) || continue
        x = left + (xval / plotmax) * pw
        push!(parts, "<line x1='$(fmt(x,1))' y1='$top' x2='$(fmt(x,1))' y2='$(top+ph)' stroke='$color' stroke-width='1.5' stroke-dasharray='4,3'/>")
        push!(parts, "<text x='$(fmt(x+3,1))' y='$(top+12)' font-size='11' fill='$color'>$label</text>")
    end

    return "<svg viewBox='0 0 $width $height' xmlns='http://www.w3.org/2000/svg' style='width:100%;height:auto;background:#fff'>" * join(parts) * "</svg>"
end

html_escape(s::AbstractString) = replace(String(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")

function write_report(path::AbstractString, stats, blocks, binwidth, extremethresh)
    io = open(path, "w")
    println(io, """<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>
    <title>Vegetation height distribution by island</title>
    <style>
      body { font-family: -apple-system, Segoe UI, Roboto, Arial, sans-serif; color: #222;
             max-width: 1080px; margin: 2em auto; padding: 0 1em; line-height: 1.5; }
      h1, h2 { color: #23395d; }
      table { border-collapse: collapse; width: 100%; font-size: 0.9em; margin: 1em 0; }
      th, td { border: 1px solid #ddd; padding: 5px 8px; text-align: right; }
      th { background: #f2f5fa; }
      td:first-child, th:first-child { text-align: left; }
      .note { color: #555; font-size: 0.92em; }
      .panel { border: 1px solid #e0e0e0; border-radius: 8px; padding: 0.5em 1em 1em;
               margin: 1.5em 0; box-shadow: 0 1px 3px rgba(0,0,0,.05); }
    </style></head><body>""")
    println(io, "<h1>Vegetation-height distribution and outliers, by island</h1>")
    println(io, "<p class='note'>Source: merged per-island mosaics in <code>data/vegetation_height</code>. " *
                "Pixels = all valid (non-NoData) cells of the 1 m vegetation-height rasters. " *
                "“Forest” below means cells above $(FOREST_MIN_HEIGHT) m (high-vegetation class). " *
                "Tukey fences are Q3 ± 1.5×IQR / 3×IQR of all valid cells. " *
                "Histogram bars are on a log10(count+1) scale; overflow cells above the histogram range are not plotted.</p>")

    println(io, "<h2>Summary</h2>")
    println(io, "<table><thead><tr><th>Island</th><th>Valid pixels</th><th>% zero</th>" *
                "<th>Median (m)</th><th>P90 (m)</th><th>P99 (m)</th><th>Max (m)</th>" *
                "<th>Tukey 1.5×IQR (m)</th><th>&gt; fence</th><th>&gt; 30 m</th><th>&gt; 40 m</th></tr></thead><tbody>")
    for s in stats
        println(io, "<tr><td>$(html_escape(s.name)) ($(s.code))</td><td>$(s.n_valid)</td>" *
                    "<td>$(fmt(s.pct_zero,1))</td><td>$(fmt(s.p50,2))</td>" *
                    "<td>$(fmt(s.p90,2))</td><td>$(fmt(s.p99,2))</td><td>$(fmt(s.max,2))</td>" *
                    "<td>$(fmt(s.fence15,2))</td><td>$(s.n_gt_fence15)</td>" *
                    "<td>$(s.above[findfirst(==(30.0), THRESHOLDS)])</td>" *
                    "<td>$(s.above[findfirst(==(40.0), THRESHOLDS)])</td></tr>")
    end
    println(io, "</tbody></table>")

    println(io, "<h2>Distributions</h2>")
    plotmax = 60.0
    for s in stats
        println(io, "<div class='panel'><h3>$(html_escape(s.name)) ($(s.code))</h3>")
        println(io, histogram_svg(s, binwidth, plotmax))
        println(io, "<p class='note'>Median $(fmt(s.p50,2)) m · P99 $(fmt(s.p99,2)) m · " *
                    "P99.9 $(fmt(s.p999,2)) m · max $(fmt(s.max,2)) m · " *
                    "Tukey fences $(fmt(s.fence15,2)) / $(fmt(s.fence30,2)) m · " *
                    "pixels &gt; 30 m: $(s.above[findfirst(==(30.0), THRESHOLDS)]) · " *
                    "&gt; 40 m: $(s.above[findfirst(==(40.0), THRESHOLDS)]) · " *
                    "&gt; 50 m: $(s.above[findfirst(==(50.0), THRESHOLDS)]) · " *
                    "negative: $(s.n_negative)</p>")

        islblocks = filter(b -> b.island == s.code, blocks)
        k30 = findfirst(==(30.0), THRESHOLDS)
        k40 = findfirst(==(40.0), THRESHOLDS)
        ranked = sort(islblocks; by=b -> -b.above[k30])
        top = ranked[1:min(5, length(ranked))]
        println(io, "<table><thead><tr><th>Row blocks with most pixels &gt; 30 m</th>" *
                    "<th>&gt; 30 m</th><th>&gt; 40 m</th><th>max (m)</th></tr></thead><tbody>")
        for b in top
            println(io, "<tr><td>rows $(b.rowstart)–$(b.rowstart+b.nrows-1) " *
                        "(y $(fmt(b.ytop,0)) … $(fmt(b.ybottom,0)))</td>" *
                        "<td>$(b.above[k30])</td><td>$(b.above[k40])</td>" *
                        "<td>$(fmt(Float64(b.maxheight),2))</td></tr>")
        end
        println(io, "</tbody></table></div>")
    end

    println(io, "<p class='note'>Pixel-level locations of all cells above $(extremethresh) m " *
                "(UTM 26N / PTRA08, EPSG:5015) are in <code>height_extreme_pixels.csv</code>. " *
                "Per-window counts are in <code>height_outlier_blocks.csv</code>.</p>")
    println(io, "</body></html>")
    close(io)
    return nothing
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    opts, flags = parse_cli()
    "help" in flags && (help(); return)

    rasterdir = get(opts, "raster-dir", joinpath("data", "vegetation_height"))
    outdir = get(opts, "outdir", joinpath("data", "vegetation_height", "qc"))
    binwidth = parse(Float64, get(opts, "bin-width", string(DEFAULT_BINWIDTH)))
    maxheight = parse(Float64, get(opts, "max-height", string(DEFAULT_MAXHEIGHT)))
    extremethresh = parse(Float64, get(opts, "extreme-threshold", string(DEFAULT_EXTREME)))
    maxextreme = parse(Int, get(opts, "max-extreme", string(DEFAULT_MAX_EXTREME)))
    blockrows = parse(Int, get(opts, "block-rows", string(DEFAULT_BLOCK_ROWS)))
    onlyisland = haskey(opts, "island") ? resolve_island(opts["island"]) : nothing

    binwidth > 0 && maxheight > binwidth && blockrows > 0 ||
        error("--bin-width, --max-height and --block-rows must be positive")
    isdir(rasterdir) || error("Not a directory: $rasterdir")

    rasters = sort([joinpath(rasterdir, f) for f in readdir(rasterdir)
                    if endswith(lowercase(f), "_vegetation_height.tif")])
    onlyisland === nothing || filter!(p -> island_code(p) == onlyisland, rasters)
    isempty(rasters) && error("No merged island rasters found in $(abspath(rasterdir))")

    nbins = ceil(Int, maxheight / binwidth)
    println("Scanning $(length(rasters)) merged rasters with $(Threads.nthreads()) thread(s): " *
            join([island_code(p) for p in rasters], ", "))
    flush(stdout)

    # One task per read window; each task opens the raster independently, which
    # keeps GDAL state thread-local and bounds memory to one window per thread.
    tasks = Tuple{String,String,Int,Int}[]
    for p in rasters
        ds = ArchGDAL.read(p)
        h = ArchGDAL.height(ds)
        code = island_code(p)
        ArchGDAL.destroy(ds)
        for rowstart in 1:blockrows:h
            push!(tasks, (code, p, rowstart, min(blockrows, h - rowstart + 1)))
        end
    end
    blocks = Vector{BlockResult}(undef, length(tasks))
    Threads.@threads for i in eachindex(tasks)
        code, p, rowstart, nrows = tasks[i]
        blocks[i] = process_block(p, code, rowstart, nrows, binwidth, nbins,
                                  maxheight, extremethresh)
    end

    codes = sort(unique([b.island for b in blocks]))
    stats = [summarize_island(c, filter(b -> b.island == c, blocks), binwidth, maxheight)
             for c in codes]

    mkpath(outdir)
    h, r = summary_rows(stats)
    write_csv(joinpath(outdir, "height_distribution_by_island.csv"), h, r)
    h, r = histogram_rows(stats, binwidth)
    write_csv(joinpath(outdir, "height_histogram.csv"), h, r)
    h, r = block_rows(blocks)
    write_csv(joinpath(outdir, "height_outlier_blocks.csv"), h, r)
    h, r, nextreme, truncated = extreme_rows(blocks, extremethresh, maxextreme)
    write_csv(joinpath(outdir, "height_extreme_pixels.csv"), h, r)
    write_report(joinpath(outdir, "height_distribution_report.html"), stats, blocks,
                 binwidth, extremethresh)

    println("\nPer-island summary (merged rasters):")
    println("Island  Valid px      Median    P99    Max   >30 m   >40 m  Tukey(1.5IQR)")
    for s in stats
        @printf("%-6s %10d %8.2f %6.2f %6.2f %8d %8d %10.2f\n",
                s.code, s.n_valid, s.p50, s.p99, s.max,
                s.above[findfirst(==(30.0), THRESHOLDS)],
                s.above[findfirst(==(40.0), THRESHOLDS)], s.fence15)
    end
    println("\nWrote 5 files to $(abspath(outdir))")
    println("  extreme pixels: $nextreme" * (truncated ? " (truncated to --max-extreme)" : ""))
    return nothing
end

(abspath(PROGRAM_FILE) == @__FILE__) && main()
