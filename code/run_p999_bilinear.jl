#!/usr/bin/env julia
#
# run_p999_bilinear.jl
#
# Convenience runner: clamp every merged island raster to P99.9 and fill the
# pixels removed as outliers by bilinear interpolation. Settings are embedded
# here so they cannot be lost when the process is run in the background.
#
# Optional environment variables:
#   AZORES_ISLAND=CODE   process a single island
#   AZORES_OUTDIR=DIR    write the outputs somewhere else
#
#   julia --project code/run_p999_bilinear.jl

args = [
    "--percentile", "0.999",
    "--interpolate",
    "--raster-dir", get(ENV, "AZORES_RASTER_DIR", "data/vegetation_height"),
    "--outdir", get(ENV, "AZORES_OUTDIR", "data/vegetation_height"),
]
haskey(ENV, "AZORES_ISLAND") && append!(args, ["--island", ENV["AZORES_ISLAND"]])

empty!(ARGS)
append!(ARGS, args)

include(joinpath(@__DIR__, "filter_raster_p99.jl"))
main()
