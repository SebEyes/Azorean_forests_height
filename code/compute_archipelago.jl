#!/usr/bin/env julia
#
# compute_archipelago.jl
#
# Tile-by-tile pipeline that turns the DGT ACORES-LAZ point clouds into a
# 1 m vegetation-height (canopy height model) raster for each island.
#
# For every tile, in order:
#   1. download the single .laz file into a temporary work directory,
#   2. compute the 1 m vegetation-height raster with laz_to_vegetation_height.jl,
#   3. save that raster as a GeoTIFF in --tiles-dir,
#   4. delete the .laz file immediately (disk use stays bounded),
# and once a whole island is finished:
#   5. merge all of that island's raster tiles into one GeoTIFF (in Julia).
#
# The STAC search / Keycloak login / download code is reused from
# download_azores_laz.jl (which is includable without running its CLI).
#
# Example:
#   julia code/compute_archipelago.jl --island corvo --limit 2 \
#       --tiles-dir data/vegetation_height/tiles \
#       --outdir data/vegetation_height
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
Pkg.instantiate()

using ArchGDAL

include(joinpath(@__DIR__, "download_azores_laz.jl"))

const PROCESS_SCRIPT = joinpath(@__DIR__, "laz_to_vegetation_height.jl")
const DEFAULT_NODATA = Float32(-9999)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

function parse_pipeline_cli()
    opts = Dict{String,String}()
    flags = Set{String}()
    i = 1
    while i <= length(ARGS)
        a = ARGS[i]
        if startswith(a, "--")
            if occursin('=', a)
                k, v = split(a[3:end], '=', limit=2)
                opts[k] = v
            else
                k = a[3:end]
                if k in ("help", "resume", "delete-tiles", "keep-laz")
                    push!(flags, k)
                else
                    i += 1
                    i <= length(ARGS) || error("Missing value for --$k")
                    opts[k] = ARGS[i]
                end
            end
        end
        i += 1
    end
    return opts, flags
end

function pipeline_help()
    println("""
    Usage: julia code/compute_archipelago.jl [options]

    Pipeline options:
      --island CODE       Island to process (SMG SMA TER GRA SJO PIC FAI FLO COR,
                          or a name like "corvo"); default: whole archipelago
      --limit N           Process at most N tiles (useful for a test run)
      --tiles-dir DIR     Per-tile vegetation-height GeoTIFFs
                          (default: data/vegetation_height/tiles)
      --outdir DIR        Per-island merged GeoTIFFs
                          (default: data/vegetation_height)
      --workdir DIR       Temporary directory for the downloaded .laz files
                          (default: a system temp directory, removed on exit)
      --chunk-size N      Points per chunk passed to the processor (default 2000000)
      --threads N         Worker threads per tile (default: all CPU cores)
      --resolution M      Output pixel size in metres (default 1)
      --epsg CODE         Override CRS if a LAZ header has none
      --resume            Skip tiles whose raster already exists
      --delete-tiles      Remove per-tile rasters after the island mosaic is written
      --keep-laz          Do not delete downloaded .laz files (debugging)

    Authentication (same as download_azores_laz.jl):
      --username USER / DGT_USERNAME, --password PASS / DGT_PASSWORD,
      --otp CODE / DGT_OTP, or --cookie VALUE
      --help              Show this help
    """)
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function authenticate!(opts)
    username = get(opts, "username", get(ENV, "DGT_USERNAME", ""))
    password = get(opts, "password", get(ENV, "DGT_PASSWORD", ""))
    if haskey(opts, "cookie")
        COOKIES["connect.sid"] = opts["cookie"]
        verify_login() || error("The supplied connect.sid is not authenticated (or has expired).")
        println("  ✓ Using browser session cookie")
    elseif !isempty(username) && !isempty(password)
        login!(username, password; otp=get(opts, "otp", get(ENV, "DGT_OTP", nothing)))
    else
        error("Provide --username/--password (or DGT_USERNAME/DGT_PASSWORD) or --cookie.")
    end
end

function process_tile(laz::AbstractString, raster::AbstractString, opts)
    project = normpath(joinpath(@__DIR__, ".."))
    # Each tile gets all CPU cores by default; LASzip decompression and the
    # point scatter are both parallelised across worker threads.
    nthreads = get(opts, "threads", "auto")
    cmd = `$(Base.julia_cmd()) --threads=$nthreads --project=$project $PROCESS_SCRIPT $laz -o $raster`
    for key in ("chunk-size", "resolution", "epsg")
        haskey(opts, key) && append!(cmd.exec, ["--$key", opts[key]])
    end
    # When threads are pinned explicitly, pass the same value to the processor;
    # otherwise it defaults to Threads.nthreads() chosen by --threads=auto.
    haskey(opts, "threads") && append!(cmd.exec, ["--threads", nthreads])
    run(cmd)
    return nothing
end

"""
Merge per-tile GeoTIFFs into one raster. Overlapping pixels take the maximum
vegetation height. Implemented entirely with ArchGDAL (windowed reads/writes),
so only one tile raster is held in memory at a time.
"""
function merge_rasters(paths::Vector{String}, outpath::AbstractString;
                       nodata::Float32=DEFAULT_NODATA)
    isempty(paths) && error("No raster tiles to merge into $outpath")

    x0min, x1max = Inf, -Inf
    ytopmax, y0min = -Inf, Inf
    resolution, wkt = nothing, nothing

    for p in paths
        ds = ArchGDAL.read(p)
        gt = ArchGDAL.getgeotransform(ds)
        x0, rx, _, ytop, _, ry = gt
        x1 = x0 + rx * ArchGDAL.width(ds)
        y0 = ytop + ry * ArchGDAL.height(ds)
        x0min = min(x0min, x0)
        x1max = max(x1max, x1)
        ytopmax = max(ytopmax, ytop)
        y0min = min(y0min, y0)
        if resolution === nothing
            resolution = abs(rx)
            wkt = ArchGDAL.getproj(ds)
        end
        ArchGDAL.destroy(ds)
    end

    ncols = ceil(Int, (x1max - x0min) / resolution)
    nrows = ceil(Int, (ytopmax - y0min) / resolution)
    mkpath(dirname(abspath(outpath)))
    options = ["COMPRESS=DEFLATE", "TILED=YES", "BLOCKXSIZE=256",
               "BLOCKYSIZE=256", "PREDICTOR=3"]

    ArchGDAL.create(outpath; driver=ArchGDAL.getdriver("GTiff"),
                    width=ncols, height=nrows, nbands=1, dtype=Float32,
                    options=options) do ods
        ArchGDAL.setgeotransform!(ods, [x0min, resolution, 0.0, ytopmax, 0.0, -resolution])
        ArchGDAL.setproj!(ods, wkt)
        band = ArchGDAL.getband(ods, 1)
        ArchGDAL.setnodatavalue!(band, nodata)

        for p in paths
            ids = ArchGDAL.read(p)
            gt = ArchGDAL.getgeotransform(ids)
            offx = round(Int, (gt[1] - x0min) / resolution)
            offy = round(Int, (ytopmax - gt[4]) / resolution)
            tile = ArchGDAL.read(ids, 1)          # laid out [x+1, y+1] => (nx, ny)
            nx, ny = size(tile)
            current = ArchGDAL.read(ods, 1, offx, offy, nx, ny)
            merged = similar(tile)
            @inbounds for i in eachindex(tile)
                a, b = tile[i], current[i]
                merged[i] = a == nodata ? b : (b == nodata ? a : max(a, b))
            end
            ArchGDAL.write!(ods, merged, 1, offx, offy, nx, ny)
            ArchGDAL.destroy(ids)
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    opts, flags = parse_pipeline_cli()
    "help" in flags && (pipeline_help(); return)

    island      = get(opts, "island", nothing)
    limit       = haskey(opts, "limit") ? parse(Int, opts["limit"]) : nothing
    tilesdir    = get(opts, "tiles-dir", joinpath("data", "vegetation_height", "tiles"))
    outdir      = get(opts, "outdir", joinpath("data", "vegetation_height"))
    workdir     = get(opts, "workdir", mktempdir(prefix="azores_laz_"))
    resume      = "resume" in flags
    deletetiles = "delete-tiles" in flags
    keeplaz     = "keep-laz" in flags

    mkpath(tilesdir)
    mkpath(outdir)
    mkpath(workdir)

    println("Querying STAC catalogue for $COLLECTION ...")
    items = collect_items(; island=island, limit=limit)
    isempty(items) && (println("No tiles matched."); return)
    print_summary(items)

    println("Authenticating ...")
    authenticate!(opts)

    # Download tokens are bound to the session that searched the catalogue,
    # so re-run the search now that we are logged in.
    items = collect_items(; island=island, limit=limit)

    byisland = Dict{String,Vector{Any}}()
    for item in items
        code = island_code(String(item.id))
        push!(get!(byisland, code, Any[]), item)
    end

    println("\nComputing 1 m vegetation height")
    println("  LAZ staging   : $workdir")
    println("  tile rasters  : $(abspath(tilesdir))")
    println("  island mosaics: $(abspath(outdir))\n")

    allfailed = String[]
    for code in sort(collect(keys(byisland)))
        its = sort(byisland[code]; by=i -> String(i.id))
        println("=== $code ($(get(ISLAND_NAMES, code, code))): $(length(its)) tile(s) ===")
        rasters = String[]
        failed = String[]

        for (n, item) in enumerate(its)
            id = String(item.id)
            raster = joinpath(tilesdir, id * "_vegheight.tif")

            if resume && isfile(raster)
                println("[$n/$(length(its))] skip     $id (raster exists)")
                push!(rasters, raster)
                continue
            end

            laz = joinpath(workdir, id * OUTFILE_EXT)
            try
                println("[$n/$(length(its))] download $id")
                download_item(item; outdir=workdir)

                println("[$n/$(length(its))] compute  $id")
                process_tile(laz, raster, opts)
                isfile(raster) || error("Processor did not create $raster")
                push!(rasters, raster)

                keeplaz || rm(laz; force=true)
            catch e
                push!(failed, id)
                push!(allfailed, id)
                println(stderr, "[$n/$(length(its))] FAIL    $id : $e")
                keeplaz || rm(laz; force=true)
                isfile(raster) && !resume && rm(raster; force=true)
            end
        end

        if isempty(rasters)
            @warn "No rasters were produced for $code; skipping the merge."
            continue
        end

        merged = joinpath(outdir, code * "_vegetation_height.tif")
        println("Merging $(length(rasters)) raster(s) -> $merged")
        try
            merge_rasters(rasters, merged)
            println("  ✓ $merged")
            if deletetiles
                for r in rasters
                    rm(r; force=true)
                end
            end
        catch e
            println(stderr, "Merge failed for $code: $e")
            push!(allfailed, code * " (merge)")
        end

        isempty(failed) || println("  Failed tiles for $code: ", join(failed, ", "))
        println()
    end

    keeplaz || rm(workdir; force=true, recursive=true)
    println("Done.")
    isempty(allfailed) ||
        println("Failures:\n  " * join(allfailed, "\n  ") *
                "\nRe-run with --resume to retry only the missing tiles.")
end

(abspath(PROGRAM_FILE) == @__FILE__) && main()
