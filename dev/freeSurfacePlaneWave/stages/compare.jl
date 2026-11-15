# dev/freeSurfacePlaneWave/stages/compare.jl
analyticalFineTimes = collect(0.0:0.01:duration)
analyticalResult = planeWaveNormalIncidenceP(analyticalFineTimes, observerDepths;
    z0=z0, vp=vp0, amplitude=planewaveAmplitude, pulse=rickerSource,
    include_image=applyFreeSurface)

function _depth_trace(uz, zaxis, depth, timesAxis)
    iz = argmin(abs.(zaxis .- depth))
    midColumn = cld(size(uz, 1), 2)
    (time=Float64.(timesAxis), values=Float64.(uz[midColumn, iz, :]))
end

# Exact closed form: no interpolation needed, evaluate planeWaveUz directly
# at whatever time the simulation reports.
_analytical_at(depth, t) = planeWaveUz(depth, t; includeImage=applyFreeSurface)

compareFigure = Figure(size=(1100, 220 * length(observerDepths)))
compareMetrics = NamedTuple[]
for (row, depth) in enumerate(observerDepths)
    axis = Axis(compareFigure[row, 1]; ylabel="z=$(round(depth/1e3; digits=1)) km",
        xlabel=row == length(observerDepths) ? "time (s)" : "")
    optTrace = _depth_trace(uzOPT, zOPT, depth, optTimes)
    opt5Trace = _depth_trace(uzOPT5, zOPT, depth, opt5Times)
    fdTrace = _depth_trace(uzFD, fdZ, depth, fdTimes)
    analyticalColumn = argmin(abs.(analyticalResult.depths .- depth))
    lines!(axis, analyticalResult.time, analyticalResult.displacement[:, analyticalColumn];
        color=:black, linewidth=2, label="analytical")
    lines!(axis, fdTrace.time, fdTrace.values; label="FD3")
    lines!(axis, optTrace.time, optTrace.values; label="OPT3")
    lines!(axis, opt5Trace.time, opt5Trace.values; label="OPT5", linestyle=:dash)
    row == 1 && axislegend(axis; position=:rb)
    function metric(trace)
        matched = [_analytical_at(depth, tt) for tt in trace.time]
        denom = max(norm(matched), eps())
        (max_abs_error=maximum(abs.(trace.values .- matched)),
         relative_error=norm(trace.values .- matched) / denom,
         max_analytical=maximum(abs, matched), max_numerical=maximum(abs, trace.values))
    end
    push!(compareMetrics, (depth=depth, FD3=metric(fdTrace), OPT3=metric(optTrace), OPT5=metric(opt5Trace)))
end
devSave(compareFigure, applyFreeSurface ? "vs_analytical_uz" : "vs_analytical_uz_nosurface")
foreach(display, compareMetrics)
metricsPath = joinpath(DEV_OUTPUT_DIR[], DEV_STEP[], "vs_analytical_metrics.jld2")
jldsave(metricsPath; compareMetrics, observerDepths, optCacheFile, fdCacheFile)
println("metrics -> ", metricsPath)
