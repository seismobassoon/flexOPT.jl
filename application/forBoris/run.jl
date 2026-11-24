# application/forBoris/run.jl -- OPT3 / OPT5 vs your own code (SeismicWavePC.jl
# or anything else): point-force, flat free surface, homogeneous half-space.
#
# Thin entry point to dev/freeSurfacePointSource/run.jl (same pattern as
# Elastic_Cartesian2D_Time -> dev/kirishima_oct2026):
#   julia --project=. application/forBoris/run.jl run           # quick (60 km half-width, 9 s)
#   julia --project=. application/forBoris/run.jl run --full    # fuller (90 km half-width, 14 s)
# See application/forBoris/README.md for the physical setup, units, and how
# this maps onto SeismicWavePC.jl's ForceSource / Station / simulate API.

global APP_ARGS = copy(@isdefined(APP_ARGS) ? APP_ARGS : ARGS)
Base.include(Main, joinpath(@__DIR__, "..", "..", "dev", "freeSurfacePointSource", "run.jl"))
