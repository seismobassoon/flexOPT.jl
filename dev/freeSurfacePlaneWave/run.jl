# dev/freeSurfacePlaneWave/run.jl -- normal-incidence plane P wave on a flat
# free surface, homogeneous medium, exact closed-form analytical reference
# (method of images, dev/freeSurfacePlaneWave/stages/analytical.jl) -- no
# SPECFEM2D, no artificial-boundary ambiguity (the source/model are x-
# uniform by construction, so there is no lateral wave to reflect off
# anything; the only truncation is the bottom one, kept duration-safe).
#
# Built to check a suspicion raised by dev/kirishima_oct2026's Step 3
# (flat free surface, point source, SPECFEM2D reference): OPT3/OPT5 traces
# there sit at ~zero through the whole physical arrival window and only
# move with a late numerical blow-up, while FD3/convFD3 show real wave
# arrivals -- i.e. a likely bug in the :additive elasticTractionFree2D free-
# surface closure itself, not (only) a SPECFEM-boundary-reflection problem.
#
#   julia --project=. dev/freeSurfacePlaneWave/run.jl normal
#   julia --project=. dev/freeSurfacePlaneWave/run.jl normal --redo=opt
# From a REPL (fastest):
#   include("dev/freeSurfacePlaneWave/run.jl")
#   planewave("normal", "--redo=fd")
# Figures and metrics: dev/freeSurfacePlaneWave/outputs/<step>/

include(joinpath(@__DIR__, "..", "common", "devBootstrap.jl"))
PLANEWAVE_STAGES = joinpath(@__DIR__, "stages")
planewave(args...) = devRun(joinpath(@__DIR__, "run.jl"), args...)

DEV = parseDevArgs(@isdefined(DEV_ARGS) ? DEV_ARGS : ARGS;
    steps=["normal"], defaultStep="normal")
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
println("Plane-wave free-surface dev run: ", DEV)

include(joinpath(PLANEWAVE_STAGES, "env.jl"))
include(joinpath(@__DIR__, "codeVersions.jl"))

for step in DEV.steps
    DEV_STEP[] = step
    if step == "normal"
        include(joinpath(PLANEWAVE_STAGES, "config.jl"))
        include(joinpath(PLANEWAVE_STAGES, "opt.jl"))
        include(joinpath(PLANEWAVE_STAGES, "fd.jl"))
        include(joinpath(PLANEWAVE_STAGES, "analytical.jl"))
        include(joinpath(PLANEWAVE_STAGES, "compare.jl"))
    end
end
