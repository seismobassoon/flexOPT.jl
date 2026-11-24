# dev/freeSurfacePlaneWave/run.jl -- normal-incidence plane P wave on a flat
# free surface, homogeneous medium, exact closed-form analytical reference
# (dev/freeSurfacePlaneWave/stages/analytical.jl) -- no SPECFEM2D, no
# artificial-boundary ambiguity (the model/field are x-uniform by
# construction, so there is no lateral wave to reflect off anything; the
# only truncation is the bottom one, kept duration-safe).
#
# NO BODY FORCE / POINT SOURCE / SHEET SOURCE ANYWHERE (changed 2026-11-15,
# on request: avoid any source-calibration question entirely). The plane
# wave is injected as an EXACT initial condition -- config.jl's planeWaveUz,
# two counter-propagating copies of the same pulse -- seeded at t=-dt and
# t=0, then propagated with zero forcing from there on. See config.jl and
# stages/analytical.jl for the closed form and its derivation.
#
# Built to check a suspicion raised by dev/kirishima_oct2026's Step 3
# (flat free surface, point source, SPECFEM2D reference): OPT3/OPT5 traces
# there sit at ~zero through the whole physical arrival window and only
# move with a late numerical blow-up, while FD3/convFD3 show real wave
# arrivals -- i.e. a likely bug in the :additive elasticTractionFree2D free-
# surface closure itself, not (only) a SPECFEM-boundary-reflection problem.
#
#   julia --project=. dev/freeSurfacePlaneWave/run.jl normal
#   julia --project=. dev/freeSurfacePlaneWave/run.jl nosurface   # control: no free surface at all
#   julia --project=. dev/freeSurfacePlaneWave/run.jl normal --redo=opt
# From a REPL (fastest):
#   include("dev/freeSurfacePlaneWave/run.jl")
#   planewave("normal", "--redo=fd")
#   planewave("nosurface")
# Figures and metrics: dev/freeSurfacePlaneWave/outputs/<step>/

include(joinpath(@__DIR__, "..", "common", "devBootstrap.jl"))
PLANEWAVE_STAGES = joinpath(@__DIR__, "stages")
planewave(args...) = devRun(joinpath(@__DIR__, "run.jl"), args...)

# CFL-sweep diagnostic (dt-halving test on the free-surface closure): set
# PLANEWAVE_CFL_FRACTION[] = 0.10 (default 0.20) before calling planewave(...)
# or include-ing this file, then --redo=opt (the hash doesn't see this value,
# only CODE_VERSION["opt"] does -- bump that in codeVersions.jl too, or pass
# --redo=opt, when sweeping it more than once in the same session).
if !@isdefined(PLANEWAVE_CFL_FRACTION)
    const PLANEWAVE_CFL_FRACTION = Ref(0.20)
end
# boundary_weight sweep: set PLANEWAVE_BOUNDARY_WEIGHT[] = :match_volume
# (default 1.0, the current dev/kirishima_oct2026 convention) before calling
# planewave("normal") to test whether a flat additive weight of 1.0 is
# imbalancing the free-surface closure (see opt_build.jl).
if !@isdefined(PLANEWAVE_BOUNDARY_WEIGHT)
    const PLANEWAVE_BOUNDARY_WEIGHT = Ref{Any}(1.0)
end
# overlap-mode sweep: set PLANEWAVE_OVERLAP_MODE[] = :replace (default
# :additive, the current dev/kirishima_oct2026 convention) to test
# overlapBoundaryLinearSystem's own default mode -- fully substituting the
# interior equation at free-surface rows instead of adding the boundary
# equation on top of it.
if !@isdefined(PLANEWAVE_OVERLAP_MODE)
    const PLANEWAVE_OVERLAP_MODE = Ref{Symbol}(:additive)
end
# Lateral edges: pinning only the single outermost column each step (the
# original design) left too thin a buffer for this narrow (nxHalfWidth=10)
# domain -- a persistent edge/interior mismatch leaked all the way to the
# measured center column within the run (xSpreadFD/xSpreadOPT comparable
# to the signal itself: confirmed 2026-11-24, see codeVersions.jl). Pinning
# a wider band fixes it at no extra grid cost (same domain, same dx) --
# width=8 (of nxHalfWidth=10) brought xSpreadFD from ~1.1-1.4 down to
# ~0.047 and OPT3/OPT5 down to near machine precision in the "nosurface"
# control. See stages/opt_propagate.jl / stages/fd_propagate.jl.
if !@isdefined(PLANEWAVE_EDGE_PIN_WIDTH)
    const PLANEWAVE_EDGE_PIN_WIDTH = Ref(8)
end

DEV = parseDevArgs(@isdefined(DEV_ARGS) ? DEV_ARGS : ARGS;
    steps=["normal", "nosurface"], defaultStep="normal")
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
println("Plane-wave free-surface dev run: ", DEV)

include(joinpath(PLANEWAVE_STAGES, "env.jl"))
include(joinpath(@__DIR__, "codeVersions.jl"))

for step in DEV.steps
    DEV_STEP[] = step
    if step in ("normal", "nosurface")
        # "normal": flat free surface at z=0 (the actual closure under test).
        # "nosurface": SAME plane-wave IC, SAME grid style, but no free
        # surface at all -- fully solid medium, Cerjan on both z ends,
        # compared to the pure direct-wave (no image term) analytical
        # solution. Isolates whether the growing, x-localised spurious
        # signal found in "normal" needs a free surface to appear at all.
        include(joinpath(PLANEWAVE_STAGES, "config.jl"))
        include(joinpath(PLANEWAVE_STAGES, "opt.jl"))
        include(joinpath(PLANEWAVE_STAGES, "fd.jl"))
        include(joinpath(PLANEWAVE_STAGES, "specfem.jl"))
        include(joinpath(PLANEWAVE_STAGES, "analytical.jl"))
        include(joinpath(PLANEWAVE_STAGES, "compare.jl"))
    end
end
