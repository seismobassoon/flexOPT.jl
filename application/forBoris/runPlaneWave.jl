# application/forBoris/runPlaneWave.jl -- OPT3 / OPT5 vs the EXACT analytical
# solution for a normal-incidence plane P wave on a flat free surface.
#
# Thin entry point to dev/freeSurfacePlaneWave/forBorisRun.jl (same pattern
# as application/forBoris/run.jl -> dev/freeSurfacePointSource/run.jl):
#   julia --project=. application/forBoris/runPlaneWave.jl run           # quick
#   julia --project=. application/forBoris/runPlaneWave.jl run --full    # fuller
# Always runs BOTH the actual free-surface test ("normal") and the no-
# free-surface control ("nosurface") in one invocation.
# See application/forBoris/README.md for the physical setup and the known
# open issue with OPT3's free-surface closure that this benchmark exists to
# pin down.

global APP_ARGS = copy(@isdefined(APP_ARGS) ? APP_ARGS : ARGS)
Base.include(Main, joinpath(@__DIR__, "..", "..", "dev", "freeSurfacePlaneWave", "forBorisRun.jl"))
