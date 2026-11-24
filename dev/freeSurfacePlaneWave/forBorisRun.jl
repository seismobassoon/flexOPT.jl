# dev/freeSurfacePlaneWave/forBorisRun.jl -- OPT3 / OPT5 only, normal-incidence
# plane P wave on a flat free surface, vs the EXACT closed-form analytical
# solution (no SPECFEM2D, no FD3, no external dependency beyond flexOPT
# itself) -- the plane-wave half of application/forBoris/, alongside the
# point-force benchmark in dev/freeSurfacePointSource/run.jl.
#
# A trimmed, non-interactive reuse of dev/freeSurfacePlaneWave/run.jl's own
# stages (config.jl, opt.jl -> opt_build.jl/opt_propagate.jl, analytical.jl
# -- read in full before writing this) with FD3 and SPECFEM2D dropped and
# compare.jl replaced by a simpler OPT3/OPT5-vs-analytical-only version
# below, so this needs nothing beyond flexOPT to run. This REUSES the same
# semi-symbolic/opt caches as the full dev/freeSurfacePlaneWave/run.jl (same
# DEV_STEP names "normal"/"nosurface", same physical config), so if those
# have already been computed in this repository, this goes straight to
# plotting; a cache miss recomputes them (same cost as dev/freeSurfacePlaneWave
# itself).
#
# Why this design beats the point-force benchmark for an independent
# cross-check: there is NO body force / point source / Green's-function
# calibration anywhere -- the plane wave is injected as an EXACT closed-form
# initial condition (two counter-propagating copies of the same pulse; see
# stages/analytical.jl), so there is no source-amplitude or source-time-
# function CONVENTION to translate between codes at all, unlike the point-
# force case. Reproducing this in any other code just means: homogeneous
# half-space, same vp/vs/rho, a purely vertical x-independent u_z(z,0)
# initial field, free surface at z=0, propagate with zero forcing.
#
# KNOWN OPEN ISSUE, surfaced here deliberately, not hidden: the "normal"
# (actual free-surface) run's OPT3 closure (the additive elasticTractionFree2D
# construction) is known, as of this writing, to blow up partway through the
# run (flagged in dev/freeSurfacePlaneWave/stages/opt_propagate.jl; its trace
# is truncated/divergent in the figure below, not asserted away). The
# "nosurface" control (no free surface at all) has never shown this. If an
# independent code's own normal-incidence response does NOT show a divergence
# like this, that is useful evidence the issue is specific to this closure,
# not to the physical setup.
#
# Run (both "normal" and "nosurface" always run together, one invocation):
#   julia --project=. dev/freeSurfacePlaneWave/forBorisRun.jl run
#   julia --project=. dev/freeSurfacePlaneWave/forBorisRun.jl run --full
# REPL:
#   julia> include("dev/freeSurfacePlaneWave/forBorisRun.jl")
#   julia> planewaveForBoris("run", "--full")

if !@isdefined(flexopt_root)
    import Pkg
    function find_flexopt_root(start_dir=@__DIR__)
        candidates = String[]
        haskey(ENV, "FLEXOPT_ROOT") && push!(candidates, abspath(expanduser(ENV["FLEXOPT_ROOT"])))
        directory = abspath(start_dir)
        while true
            push!(candidates, directory)
            parent = dirname(directory)
            parent == directory && break
            directory = parent
        end
        for candidate in unique(candidates)
            isfile(joinpath(candidate, "src", "flexOPT.jl")) && return candidate
        end
        error("Cannot locate flexOPT; set ENV[\"FLEXOPT_ROOT\"]")
    end
    global flexopt_root = find_flexopt_root()
    Pkg.activate(flexopt_root)
end
include(joinpath(flexopt_root, "dev", "common", "devBootstrap.jl"))
PLANEWAVE_STAGES = joinpath(@__DIR__, "stages")

planewaveForBoris(args...) = (global APP_ARGS = collect(String, string.(args));
    include(joinpath(@__DIR__, "forBorisRun.jl")))

# Same knobs as run.jl (kept so this reuses the SAME opt/semi-symbolic
# cache keys when the physical config matches): defaults are the currently-
# preferred values (v5, 2026-11-24 wider-edge-pin-band fix), not swept here.
if !@isdefined(PLANEWAVE_CFL_FRACTION)
    const PLANEWAVE_CFL_FRACTION = Ref(0.20)
end
if !@isdefined(PLANEWAVE_BOUNDARY_WEIGHT)
    const PLANEWAVE_BOUNDARY_WEIGHT = Ref{Any}(1.0)
end
if !@isdefined(PLANEWAVE_OVERLAP_MODE)
    const PLANEWAVE_OVERLAP_MODE = Ref{Symbol}(:additive)
end
if !@isdefined(PLANEWAVE_EDGE_PIN_WIDTH)
    const PLANEWAVE_EDGE_PIN_WIDTH = Ref(8)
end

DEV = parseDevArgs(@isdefined(APP_ARGS) ? APP_ARGS : ARGS; steps=["run"], defaultStep="run")
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
println("forBoris plane-wave run (OPT3/OPT5 vs exact analytical; both normal and nosurface): ", DEV)

include(joinpath(PLANEWAVE_STAGES, "env.jl"))
include(joinpath(@__DIR__, "codeVersions.jl"))

forBorisAllMetrics = Dict{String,Any}()
for step in ("normal", "nosurface")   # always both, not user-selectable here
    DEV_STEP[] = step
    include(joinpath(PLANEWAVE_STAGES, "config.jl"))
    include(joinpath(PLANEWAVE_STAGES, "opt.jl"))       # OPT3/OPT5 only (no FD3/SPECFEM2D stages included)
    include(joinpath(PLANEWAVE_STAGES, "analytical.jl"))

    # ------------------------------------------------ forBoris compare/plot
    analyticalFineTimes = collect(0.0:0.01:duration)
    analyticalResult = planeWaveNormalIncidenceP(analyticalFineTimes, observerDepths;
        z0=z0, vp=vp0, amplitude=planewaveAmplitude, pulse=pulseSource, include_image=applyFreeSurface)
    function _depth_trace_forBoris(uz, zaxis, depth, timesAxis)
        iz = argmin(abs.(zaxis .- depth))
        midColumn = cld(size(uz, 1), 2)   # x-independent by construction: any column works
        (time=Float64.(timesAxis), values=Float64.(uz[midColumn, iz, :]))
    end
    _analytical_at_forBoris(depth, t) = planeWaveUz(depth, t; includeImage=applyFreeSurface)

    forBorisCaptionLine1 = applyFreeSurface ?
        "Normal-incidence plane P wave, FLAT FREE SURFACE, u_z: OPT3 & OPT5 vs exact analytical" :
        "Normal-incidence plane P wave, NO free surface (control run), u_z: OPT3 & OPT5 vs exact analytical"
    forBorisCaptionLine2 = "Model: vp=$(round(Int,vp0)) m/s, vs=$(round(Int,vs0)) m/s, rho=$(round(Int,rho0)) kg/m3 " *
        (applyFreeSurface ? "(homogeneous half-space, flat free surface at z=0)" :
            "(homogeneous full-space, no free surface -- isolates the propagation scheme from the closure)")
    forBorisCaptionLine3 = "Incidence: 0 deg, true normal incidence -- a genuine plane wave, purely vertical and " *
        "x-independent by construction (NOT a point source: injected as an exact closed-form initial condition, " *
        "two counter-propagating copies of the same pulse -- see stages/analytical.jl; no force-amplitude or " *
        "source-time-function CONVENTION to match between codes at all)"
    forBorisCaptionLine4 = "Wavelet: Ricker, f0=$(round(sourceFrequency;digits=4)) Hz, delay=$(round(sourceDelay;digits=3)) s " *
        "-- (1-2a^2)exp(-a^2), a=pi*f0*(t-delay)"
    forBorisCaptionLine5 = applyFreeSurface ?
        "KNOWN OPEN ISSUE (not hidden): OPT3's free-surface closure is known to blow up partway through this run " *
        "-- its trace is truncated/divergent where that happens. If your code's own normal-incidence response does " *
        "NOT diverge like this, that is useful evidence the issue is specific to this closure." :
        "Control run: no free surface, so no closure is exercised here at all -- should track analytical closely."

    forBorisFigure = Figure(size=(1100, 220 * length(observerDepths) + 150))
    forBorisMetrics = NamedTuple[]
    for (row, depth) in enumerate(observerDepths)
        axis = Axis(forBorisFigure[row + 1, 1]; ylabel="z=$(round(depth/1e3; digits=1)) km",
            xlabel=row == length(observerDepths) ? "time (s)" : "")
        analyticalColumn = argmin(abs.(analyticalResult.depths .- depth))
        analyticalPeak = maximum(abs, analyticalResult.displacement[:, analyticalColumn])
        ylims!(axis, -3analyticalPeak, 3analyticalPeak)
        optTrace = _depth_trace_forBoris(uzOPT, zOPT, depth, optTimes)
        opt5Trace = _depth_trace_forBoris(uzOPT5, zOPT, depth, opt5Times)
        lines!(axis, analyticalResult.time, analyticalResult.displacement[:, analyticalColumn];
            color=:black, linewidth=2, label="analytical")
        lines!(axis, optTrace.time, optTrace.values; label="OPT3", color=:crimson)
        lines!(axis, opt5Trace.time, opt5Trace.values; label="OPT5", color=:steelblue, linestyle=:dash)
        row == 1 && axislegend(axis; position=:rb)
        function metricAt(trace)
            matched = [_analytical_at_forBoris(depth, tt) for tt in trace.time]
            denom = max(norm(matched), eps())
            (max_abs_error=maximum(abs.(trace.values .- matched)),
                relative_error=norm(trace.values .- matched) / denom,
                max_analytical=maximum(abs, matched), max_numerical=maximum(abs, trace.values))
        end
        push!(forBorisMetrics, (depth=depth, OPT3=metricAt(optTrace), OPT5=metricAt(opt5Trace)))
    end
    Label(forBorisFigure[0, 1], forBorisCaptionLine1; fontsize=14, font=:bold)
    Label(forBorisFigure[1, 1], forBorisCaptionLine2 * "\n" * forBorisCaptionLine3 * "\n" *
        forBorisCaptionLine4 * "\n" * forBorisCaptionLine5; fontsize=10, justification=:left, halign=:left)
    foreach(display, forBorisMetrics)
    devSave(forBorisFigure, "forBoris_vs_analytical_uz")

    outdir = joinpath(DEV_OUTPUT_DIR[], DEV_STEP[])
    mkpath(outdir)
    jldsave(joinpath(outdir, "forBoris_vs_analytical_metrics.jld2");
        vp0, vs0, rho0, z0, applyFreeSurface, sourceFrequency, sourceDelay, planewaveAmplitude,
        observerDepths, forBorisMetrics)
    println("saved -> ", joinpath(outdir, "forBoris_vs_analytical_metrics.jld2"))
    forBorisAllMetrics[step] = forBorisMetrics
end
println("DONE")
