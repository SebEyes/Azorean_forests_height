#!/usr/bin/env julia
#
# laz_to_vegetation_height.jl
#
# Compute a 1 m vegetation-height raster (a canopy height model, CHM) from a
# single classified LAS/LAZ tile of the DGT ACORES-LAZ dataset.
#
# Classes used (see doc/Data_source.pdf and the .qmd in code/Fabio/):
#     2  Terreno          -> ground, used to build the DTM
#     3  Vegetacao Baixa   -> low vegetation     (0 - 0.5 m)
#     4  Vegetacao Media   -> medium vegetation  (0.5 - 2 m)
#     5  Vegetacao Alta    -> high vegetation    (> 2 m)
#
# All other DGT classes (1 Outro, 6 Construcoes, 7 Ruido, 9 Agua, 26 Pontes)
# are ignored. Ground/vegetation class numbers can be overridden on the CLI;
# the per-class histogram printed for each tile makes the selection checkable.
#
# For every 1 m cell:
#   * DTM  = minimum Z of class-2 returns, gaps filled by nearest-neighbour
#            interpolation from surrounding ground cells (bounded radius).
#   * DSMv = maximum Z of classes 3/4/5 returns.
#   * height = max(0, DSMv - DTM).  Cells with no vegetation get 0; cells with
#     no ground data within the fill radius get NoData.
#
# Points are streamed in chunks via PointClouds.jl, so peak memory stays low
# even for the largest tiles (65-221 million points).
#
# Example:
#   julia code/laz_to_vegetation_height.jl \
#       data/azores_laz/2024_COR_...laz \
#       -o data/vegetation_height/tiles/2024_COR_..._vegheight.tif
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
Pkg.instantiate()

using PointClouds
using ArchGDAL

const DEFAULT_NODATA   = Float32(-9999)
const DEFAULT_RESOLUTION = 1.0
const DEFAULT_CHUNK    = 2_000_000
const DEFAULT_MAXDIST  = 10   # cells (metres) for DTM gap filling

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

function parse_cli()
    opts = Dict{String,String}()
    i = 1
    positional = String[]
    while i <= length(ARGS)
        a = ARGS[i]
        if a == "-h" || a == "--help"
            opts["help"] = "1"
        elseif startswith(a, "--")
            k = a[3:end]
            i += 1
            i <= length(ARGS) || error("Missing value for --$k")
            opts[k] = ARGS[i]
        elseif a == "-o"
            i += 1
            i <= length(ARGS) || error("Missing value for -o")
            opts["output"] = ARGS[i]
        else
            push!(positional, a)
        end
        i += 1
    end
    return positional, opts
end

function help()
    println("""
    Usage: julia code/laz_to_vegetation_height.jl INPUT.laz -o OUTPUT.tif [options]

    Options:
      --resolution M     Pixel size in metres (default $DEFAULT_RESOLUTION)
      --chunk-size N     Points per chunk (default $DEFAULT_CHUNK)
      --ground-class N   ASPRS ground class (default 2)
      --veg-classes A B …  Classes treated as vegetation (default 3 4 5)
                         NOTE: quote the value, e.g. --veg-classes "3 4 5"
      --fill-distance M  Max DTM interpolation radius in metres (default $DEFAULT_MAXDIST)
      --nodata VALUE     NoData value (default $DEFAULT_NODATA)
      --epsg CODE        Override the CRS from the LAZ header
      --threads N        Decompression/raster worker threads (default: Julia threads)
      --help             Show this help
    """)
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

"""Grid geometry, aligned to whole metres so adjacent tiles line up exactly."""
function tile_geometry(las, resolution::Float64)
    minx, miny, _ = las.coord_min
    maxx, maxy, _ = las.coord_max
    x0    = floor(minx / resolution) * resolution
    ytop  = ceil(maxy / resolution) * resolution
    ncols = ceil(Int, (maxx - x0) / resolution)
    nrows = ceil(Int, (ytop - miny) / resolution)
    (ncols < 1 || nrows < 1) && error("Invalid tile bounds")
    return x0, ytop, ncols, nrows
end

"""
    _rasterize_block!(ground, veg, pts, ranges, chunks, classcounts_row,
                      sx, sy, sz, ox, oy, oz, x0, ytop, inv_res,
                      ncols, nrows, ground_class, vegset,
                      progress, printlock, npoints)

Decode and scatter one worker's worth of chunks.

Every parameter is passed explicitly rather than captured from a closure:
Julia boxes variables captured by `Threads.@threads`, which made the hot loop
~15x slower.  A separate, fully type-specialised function avoids that.
"""
function _rasterize_block!(ground, veg, pts, ranges, chunks, classcounts_row,
                           sx, sy, sz, ox, oy, oz, x0, ytop, inv_res,
                           ncols, nrows, ground_class, vegset,
                           progress, printlock, npoints)
    nread = nground = nveg = 0
    IO_ = PointClouds.IO
    for ci in chunks
        start, stop = ranges[ci]
        # Decode this chunk once, extracting both coordinates and class.
        raw = collect(IO_.getattrs((Val(:coords), Val(:metadata)), pts, start:stop))
        @inbounds for k in eachindex(raw)
            cr, mr = raw[k]
            c = Int(IO_.meta_class(mr))
            0 <= c <= 255 && (classcounts_row[c + 1] += 1)
            isg = c == ground_class
            isv = !isg && c in vegset
            nread += 1
            (isg || isv) || continue

            x = cr[1] * sx + ox
            y = cr[2] * sy + oy
            zz = Float32(cr[3] * sz + oz)
            col = floor(Int, (x - x0) * inv_res) + 1
            row = floor(Int, (ytop - y) * inv_res) + 1
            (1 <= col <= ncols && 1 <= row <= nrows) || continue

            if isg
                zz < ground[col, row] && (ground[col, row] = zz)
                nground += 1
            else
                zz > veg[col, row] && (veg[col, row] = zz)
                nveg += 1
            end
        end
        Threads.atomic_add!(progress, stop - start + 1)
        lock(printlock) do
            print("\r  read $(progress[]) / $(npoints) points ...")
            flush(stdout)
        end
    end
    return nread, nground, nveg
end

"""
    rasterize_points(laz, las, x0, ytop, ncols, nrows, resolution,
                     chunk, ground_class, veg_classes; nthreads)

Stream points into per-cell minimum ground Z and maximum vegetation Z.

Performance notes
-----------------
* LASzip decompression is the dominant cost and is inherently sequential
  *within* a single stream. We open one independent reader per worker thread
  and let each worker decompress a contiguous block of chunks in parallel
  (`laszip_seek_point` jumps straight to chunk boundaries in LAZ 1.4, so
  seeking is cheap).  This is the CPU analogue of GPU parallelisation.
* `coordinates(las, inds)` and `classification(las, inds)` each seek and
  decompress the whole chunk; a single `getattrs` call with both attributes
  decodes every point exactly once, halving the decode work.
* Each worker writes to its own ground/vegetation buffers, reduced with
  min/max afterwards. This avoids locks in the hot point loop.
"""
function rasterize_points(laz::AbstractString, las, x0::Float64, ytop::Float64,
                          ncols::Int, nrows::Int, resolution::Float64,
                          chunk::Int, ground_class::Int, veg_classes::Vector{Int};
                          nthreads::Int=Threads.nthreads())
    npoints = length(las)
    sx, sy, sz = las.coord_scale
    ox, oy, oz = las.coord_offset
    vegset = Set{Int}(veg_classes)

    nthreads = max(1, min(nthreads, npoints > 0 ? npoints : 1))
    grounds = [fill!(Matrix{Float32}(undef, ncols, nrows), Inf) for _ in 1:nthreads]
    vegs    = [fill!(Matrix{Float32}(undef, ncols, nrows), -Inf) for _ in 1:nthreads]

    # Contiguous blocks of chunks per worker: each reader then streams forward
    # with very few seeks.
    ranges = [(s, min(s + chunk - 1, npoints)) for s in 1:chunk:npoints]
    nchunks = length(ranges)
    counts = zeros(Int, nthreads, 3)   # per worker: read, ground, vegetation
    classcounts = zeros(Int, nthreads, 256)  # per worker: ASPRS class histogram
    progress = Threads.Atomic{Int}(0)
    printlock = ReentrantLock()
    inv_res = 1.0 / resolution

    Threads.@threads for tid in 1:nthreads
        # Independent LASzip reader per worker (thread-safe: separate library
        # instances with separate point buffers).
        wlas = read(laz, LAS)
        c1 = min(nchunks, div((tid - 1) * nchunks, nthreads) + 1)
        c2 = min(nchunks, div(tid * nchunks, nthreads))
        nr, ng, nv = _rasterize_block!(
            grounds[tid], vegs[tid], wlas.points, ranges, c1:c2,
            @view(classcounts[tid, :]),
            sx, sy, sz, ox, oy, oz, x0, ytop, inv_res,
            ncols, nrows, ground_class, vegset,
            progress, printlock, npoints)
        counts[tid, 1] = nr
        counts[tid, 2] = ng
        counts[tid, 3] = nv
    end
    println()
    nread = sum(@view counts[:, 1])
    nground = sum(@view counts[:, 2])
    nveg = sum(@view counts[:, 3])

    # Reduce per-worker buffers.
    ground = grounds[1]
    veg = vegs[1]
    @inbounds for tid in 2:nthreads, i in eachindex(ground)
        ground[i] = min(ground[i], grounds[tid][i])
        veg[i] = max(veg[i], vegs[tid][i])
    end
    return ground, veg, nread, nground, nveg, vec(sum(classcounts; dims=1))
end

"""
Fill NoData ground cells by nearest-neighbour propagation from valid cells
(8-connected BFS), up to `maxdist` cells away. Pure Julia, no external tool.
"""
function fill_dtm(ground::Matrix{Float32}, nodata::Float32, maxdist::Int)
    ncols, nrows = size(ground)
    filled = similar(ground)
    dist = Matrix{Int}(undef, ncols, nrows)
    queue = CartesianIndex{2}[]
    sizehint!(queue, min(ncols * nrows, 1_000_000))

    @inbounds for j in 1:nrows, i in 1:ncols
        v = ground[i, j]
        if isfinite(v)
            filled[i, j] = v
            dist[i, j] = 0
            push!(queue, CartesianIndex(i, j))
        else
            filled[i, j] = nodata
            dist[i, j] = -1
        end
    end

    head = 1
    @inbounds while head <= length(queue)
        idx = queue[head]
        head += 1
        d = dist[idx]
        d >= maxdist && continue
        i, j = idx.I
        v = filled[idx]
        for dj in -1:1, di in -1:1
            (di == 0 && dj == 0) && continue
            ii, jj = i + di, j + dj
            (1 <= ii <= ncols && 1 <= jj <= nrows) || continue
            nid = CartesianIndex(ii, jj)
            dist[nid] < 0 || continue
            filled[nid] = v
            dist[nid] = d + 1
            push!(queue, nid)
        end
    end
    return filled
end

function build_height(veg::Matrix{Float32}, dtm::Matrix{Float32}, nodata::Float32)
    height = similar(veg)
    @inbounds for i in eachindex(height)
        d = dtm[i]
        if d == nodata
            height[i] = nodata
        elseif isfinite(veg[i])
            height[i] = max(Float32(0), veg[i] - d)
        else
            height[i] = Float32(0)
        end
    end
    return height
end

function output_wkt(laz::AbstractString, epsg)
    epsg === nothing || return ProjFormats.epsg_to_wkt(epsg)
    las = read(laz, LAS)
    wkt = getcrs(String, las)
    (wkt === nothing || isempty(wkt)) &&
        error("No CRS in the LAZ header; pass --epsg CODE.")
    return wkt
end

# Avoid an extra Proj dependency: use ArchGDAL's EPSG->WKT via import.
module ProjFormats
using ArchGDAL
function epsg_to_wkt(code::Integer)
    srs = ArchGDAL.importEPSG(Int(code))
    return ArchGDAL.toWKT(srs)
end
end

function write_raster(path::AbstractString, height::Matrix{Float32},
                      x0::Float64, ytop::Float64, resolution::Float64,
                      wkt::String, nodata::Float32)
    ncols, nrows = size(height)
    options = ["COMPRESS=DEFLATE", "TILED=YES", "BLOCKXSIZE=256",
               "BLOCKYSIZE=256", "PREDICTOR=3"]
    ArchGDAL.create(path; driver=ArchGDAL.getdriver("GTiff"),
                    width=ncols, height=nrows, nbands=1,
                    dtype=Float32, options=options) do ds
        ArchGDAL.setgeotransform!(ds, [x0, resolution, 0.0, ytop, 0.0, -resolution])
        ArchGDAL.setproj!(ds, wkt)
        band = ArchGDAL.getband(ds, 1)
        ArchGDAL.setnodatavalue!(band, nodata)
        # ArchGDAL expects the buffer indexed [x+1, y+1] (i.e. ncols × nrows),
        # which is exactly how `height` is laid out. Do NOT permutedims here.
        ArchGDAL.write!(ds, height, 1)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    positional, opts = parse_cli()
    haskey(opts, "help") && (help(); return)
    isempty(positional) && error("Missing input LAZ file (see --help)")
    laz = positional[1]
    isfile(laz) || error("Input file not found: $laz")
    haskey(opts, "output") || error("Missing -o OUTPUT.tif")
    output = opts["output"]

    resolution = haskey(opts, "resolution") ? parse(Float64, opts["resolution"]) : DEFAULT_RESOLUTION
    chunk = haskey(opts, "chunk-size") ? parse(Int, opts["chunk-size"]) : DEFAULT_CHUNK
    ground_class = haskey(opts, "ground-class") ? parse(Int, opts["ground-class"]) : 2
    veg_classes = haskey(opts, "veg-classes") ?
        [parse(Int, s) for s in split(opts["veg-classes"]) if !isempty(s)] : [3, 4, 5]
    filldist = haskey(opts, "fill-distance") ? parse(Int, opts["fill-distance"]) : DEFAULT_MAXDIST
    nodata = haskey(opts, "nodata") ? parse(Float32, opts["nodata"]) : DEFAULT_NODATA
    epsg = haskey(opts, "epsg") ? parse(Int, opts["epsg"]) : nothing
    nthreads = haskey(opts, "threads") ? parse(Int, opts["threads"]) : Threads.nthreads()

    t0 = time()
    println("Computing vegetation height from ", basename(laz),
            " ($nthreads thread(s)) ...")
    las = read(laz, LAS)
    x0, ytop, ncols, nrows = tile_geometry(las, resolution)

    ground, veg, total, nground, nveg, classcounts = rasterize_points(
        laz, las, x0, ytop, ncols, nrows, resolution, chunk,
        ground_class, veg_classes; nthreads=nthreads)
    println("  $(total) points: $nground ground, $nveg vegetation")
    println("  classes used: ground=$ground_class, vegetation=$(join(veg_classes, ","))")
    seen = findall(>(0), classcounts)
    isempty(seen) || println("  class histogram: " *
        join(["$(c - 1)=$(classcounts[c])" for c in seen], " "))
    nground == 0 &&
        @warn "No class-$ground_class ground points in $(basename(laz)); heights will be NoData."
    nveg == 0 &&
        @warn "No points of vegetation classes $(join(veg_classes, ",")) in $(basename(laz))."

    dtm = fill_dtm(ground, nodata, filldist)
    height = build_height(veg, dtm, nodata)
    wkt = output_wkt(laz, epsg)
    mkpath(dirname(abspath(output)))
    write_raster(output, height, x0, ytop, resolution, wkt, nodata)

    ncells = count(isfinite, veg)
    println("  wrote $output ($(ncols)x$(nrows) px @ $(resolution) m, " *
            "$ncells vegetated cells) in $(round(time()-t0, digits=1)) s")
end

# Only run the CLI when executed directly, not when included for testing.
(abspath(PROGRAM_FILE) == @__FILE__) && main()
