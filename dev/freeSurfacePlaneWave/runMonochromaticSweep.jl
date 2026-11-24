# dev/freeSurfacePlaneWave/runMonochromaticSweep.jl -- sweep a tapered,
# single-carrier-frequency ("monochromatic") plane-P-wave-on-a-free-surface
# benchmark across a range of frequencies, to see how OPT3/OPT5/FD3's error
# against the exact analytical solution depends on frequency (equivalently
# points-per-S-wavelength) -- rather than the single broadband Ricker pulse
# dev/freeSurfacePlaneWave/run.jl tests.
#
# Same NO BODY FORCE / NO POINT SOURCE design as run.jl: the field is still
# injected as an exact closed-form initial condition (stages/config.jl's
# planeWaveUz / stages/analytical.jl's planeWaveNormalIncidenceP), just with
# a Gaussian-enveloped cosine ("Gabor wavelet") in place of the Ricker pulse
# -- the SAME free-surface image-method closed form (two counter-
# propagating copies of the same pulse) applies unchanged to ANY pulse
# shape, so no new derivation is needed, only a different `pulse`. See
# stages/config.jl's monochromaticMode branch for how duration/belowMargin
# are widened per frequency to fit the (frequency-dependent) burst length
# safely inside the bottom-truncation margin.
#
# SPECFEM2D is intentionally NOT run here: it would mean one external
# solver run per frequency (this sweep already runs OPT3+OPT5+FD3 at 9
# frequencies x up to 2 steps), and its plane-wave injection (specfem.jl)
# has only been validated at the single run.jl frequency so far. The
# SPECFEM2D cross-check stays in run.jl/stages/specfem.jl, untouched.
#
#   julia --project=. -t auto dev/freeSurfacePlaneWave/runMonochromaticSweep.jl normal
#   julia --project=. -t auto dev/freeSurfacePlaneWave/runMonochromaticSweep.jl nosurface
#   julia --project=. -t auto dev/freeSurfacePlaneWave/runMonochromaticSweep.jl all
# From a REPL:
#   include("dev/freeSurfacePlaneWave/runMonochromaticSweep.jl")
#   monochromaticSweep("normal")
# Figures and metrics: dev/freeSurfacePlaneWave/outputs/monochromatic/
#
# Smoke-test mode (3 frequencies instead of the full 9, set before include
# or before calling monochromaticSweep(...)):
#   const PLANEWAVE_SWEEP_SMOKE_TEST = Ref(true)

include(joinpath(@__DIR__, "..", "common", "devBootstrap.jl"))
PLANEWAVE_STAGES = joinpath(@__DIR__, "stages")
monochromaticSweep(args...) = devRun(joinpath(@__DIR__, "runMonochromaticSweep.jl"), args...)

if !@isdefined(PLANEWAVE_CFL_FRACTION)
    const PLANEWAVE_CFL_FRACTION = Ref(0.20)
end
if !@isdefined(PLANEWAVE_BOUNDARY_WEIGHT)
    const PLANEWAVE_BOUNDARY_WEIGHT = Ref{Any}(1.0)
end
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
if !@isdefined(PLANEWAVE_PULSE_MODE)
    const PLANEWAVE_PULSE_MODE = Ref{Symbol}(:ricker)
end
if !@isdefined(PLANEWAVE_MONO_FREQUENCY)
    const PLANEWAVE_MONO_FREQUENCY = Ref(0.6495)
end
if !@isdefined(PLANEWAVE_MONO_CYCLES)
    const PLANEWAVE_MONO_CYCLES = Ref(3.0)
end
if !@isdefined(PLANEWAVE_SWEEP_SMOKE_TEST)
    const PLANEWAVE_SWEEP_SMOKE_TEST = Ref(false)
end

DEV = parseDevArgs(@isdefined(DEV_ARGS) ? DEV_ARGS : ARGS;
    steps=["normal", "nosurface"], defaultStep="normal")
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
println("Plane-wave MONOCHROMATIC frequency sweep: ", DEV, "  smoke_test=", PLANEWAVE_SWEEP_SMOKE_TEST[])

include(joinpath(PLANEWAVE_STAGES, "env.jl"))
include(joinpath(@__DIR__, "codeVersions.jl"))

# Points-per-S-wavelength targets this sweep walks through -- spans from
# somewhat better than dev/freeSurfacePlaneWave/run.jl's own single Ricker
# test (ppw=21.33 at f=0.6495 Hz, included below as a consistency check
# against that run) down past minimumAcceptedPointsPerSWavelength (10), to
# see where/how fast OPT3/OPT5/FD3's error grows as resolution degrades.
fullSweepPPW = [27.0, 21.33, 18.0, 15.0, 12.5, 10.5, 9.0, 7.5, 6.4]
smokeSweepPPW = [21.33, 12.5, 7.5]
sweepPointsPerWavelength = PLANEWAVE_SWEEP_SMOKE_TEST[] ? smokeSweepPPW : fullSweepPPW
monoCyclesSweep = 3.0   # fixed across the sweep; see config.jl's
    # monochromaticMode branch for how this sets envelope width/duration/
    # belowMargin. 3 cycles is a practical compromise: narrow enough a
    # spectrum to call "monochromatic" (~1/3 relative bandwidth for a
    # Gabor atom), wide enough that the burst+propagation window still
    # fits comfortably under this geometry's bottom-truncation budget even
    # at the lowest swept frequency.

# One reference :ricker inclusion first, purely to read off vs0/dx (fixed,
# frequency-independent constants) without duplicating config.jl's own
# numbers here.
PLANEWAVE_PULSE_MODE[] = :ricker
include(joinpath(PLANEWAVE_STAGES, "config.jl"))
referenceVs0, referenceDxValue = vs0, dx

PLANEWAVE_PULSE_MODE[] = :monochromatic
PLANEWAVE_MONO_CYCLES[] = monoCyclesSweep

using LinearAlgebra: norm

function _depthTrace(uz, zaxis, depth, timesAxis)
    iz = argmin(abs.(zaxis .- depth))
    midColumn = cld(size(uz, 1), 2)
    (time=Float64.(timesAxis), values=Float64.(uz[midColumn, iz, :]))
end

function _traceMetric(trace, depth, applyFS)
    matched = [planeWaveUz(depth, tt; includeImage=applyFS) for tt in trace.time]
    denom = max(norm(matched), eps())
    (relative_error=norm(trace.values .- matched) / denom,
     max_abs_error=maximum(abs.(trace.values .- matched)),
     amplitude_ratio=maximum(abs, trace.values) / max(maximum(abs, matched), eps()))
end

function _aggregateMetric(uz, zaxis, timesAxis, depths, applyFS)
    perDepth = [_traceMetric(_depthTrace(uz, zaxis, depth, timesAxis), depth, applyFS) for depth in depths]
    worst = argmax([m.relative_error for m in perDepth])
    (relative_error=perDepth[worst].relative_error,
     max_abs_error=maximum(m.max_abs_error for m in perDepth),
     amplitude_ratio=perDepth[worst].amplitude_ratio,
     perDepth=perDepth)
end

for step in DEV.steps
    DEV_STEP[] = step
    allMetrics = NamedTuple[]
    for ppw in sweepPointsPerWavelength
        PLANEWAVE_MONO_FREQUENCY[] = referenceVs0 / (referenceDxValue * ppw)
        include(joinpath(PLANEWAVE_STAGES, "config.jl"))
        include(joinpath(PLANEWAVE_STAGES, "opt.jl"))
        include(joinpath(PLANEWAVE_STAGES, "fd.jl"))
        fd3Metric = _aggregateMetric(uzFD, fdZ, fdTimes, observerDepths, applyFreeSurface)
        opt3Metric = _aggregateMetric(uzOPT, zOPT, optTimes, observerDepths, applyFreeSurface)
        opt5Metric = _aggregateMetric(uzOPT5, zOPT, opt5Times, observerDepths, applyFreeSurface)
        push!(allMetrics, (ppw=ppw, frequency=sourceFrequency, duration=duration, belowMargin=belowMargin,
            FD3=fd3Metric, OPT3=opt3Metric, OPT5=opt5Metric))
        @show step ppw sourceFrequency duration fd3Metric.relative_error opt3Metric.relative_error opt5Metric.relative_error
    end

    summaryFigure = Figure(size=(900, 700))
    axErr = Axis(summaryFigure[1,1]; xlabel="points per S-wavelength",
        ylabel="max relative error (worst observer depth)", yscale=log10,
        title="$(step): monochromatic-sweep error vs resolution ($(monoCyclesSweep)-cycle Gabor burst)")
    ppwValues = [m.ppw for m in allMetrics]
    lines!(axErr, ppwValues, [m.FD3.relative_error for m in allMetrics]; label="FD3")
    scatter!(axErr, ppwValues, [m.FD3.relative_error for m in allMetrics])
    lines!(axErr, ppwValues, [m.OPT3.relative_error for m in allMetrics]; label="OPT3")
    scatter!(axErr, ppwValues, [m.OPT3.relative_error for m in allMetrics])
    lines!(axErr, ppwValues, [m.OPT5.relative_error for m in allMetrics]; label="OPT5", linestyle=:dash)
    scatter!(axErr, ppwValues, [m.OPT5.relative_error for m in allMetrics])
    vlines!(axErr, [minimumAcceptedPointsPerSWavelength]; color=:gray, linestyle=:dot, label="design limit (10 pts/λ)")
    axislegend(axErr; position=:rt)

    axAmp = Axis(summaryFigure[2,1]; xlabel="points per S-wavelength",
        ylabel="amplitude ratio (numerical/analytical, worst depth)")
    lines!(axAmp, ppwValues, [m.FD3.amplitude_ratio for m in allMetrics]; label="FD3")
    scatter!(axAmp, ppwValues, [m.FD3.amplitude_ratio for m in allMetrics])
    lines!(axAmp, ppwValues, [m.OPT3.amplitude_ratio for m in allMetrics]; label="OPT3")
    scatter!(axAmp, ppwValues, [m.OPT3.amplitude_ratio for m in allMetrics])
    lines!(axAmp, ppwValues, [m.OPT5.amplitude_ratio for m in allMetrics]; label="OPT5", linestyle=:dash)
    scatter!(axAmp, ppwValues, [m.OPT5.amplitude_ratio for m in allMetrics])
    hlines!(axAmp, [1.0]; color=:gray, linestyle=:dot)
    axislegend(axAmp; position=:rt)

    devSave(summaryFigure, "error_vs_frequency"; subdir=joinpath("monochromatic", step))
    metricsPath = joinpath(DEV_OUTPUT_DIR[], "monochromatic", "$(step)_monochromatic_metrics.jld2")
    mkpath(dirname(metricsPath))
    jldsave(metricsPath; allMetrics, sweepPointsPerWavelength, monoCyclesSweep)
    println("monochromatic sweep metrics -> ", metricsPath)
end
