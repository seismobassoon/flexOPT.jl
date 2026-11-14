# dev/freeSurfacePlaneWave/stages/compare.jl
analyticalFineTimes = collect(0.0:0.01:duration)
analyticalResult = freeSurfaceNormalIncidenceP(analyticalFineTimes, observerDepths;
    z0=z0, vp=vp0, rho=rho0, dx=dx, line_force=sourceForce,
    force_time_function=rickerSource, source_frequency=sourceFrequency)

function _depth_trace(uz, zaxis, depth, timesAxis)
    iz = argmin(abs.(zaxis .- depth))
    midColumn = cld(size(uz, 1), 2)
    (time=Float64.(timesAxis), values=Float64.(uz[midColumn, iz, :]))
end

function _analytical_at(depth, t)
    col = argmin(abs.(analyticalResult.depths .- depth))
    _sample_linear_zero(analyticalResult.time, analyticalResult.displacement[:, col], t)
end

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
devSave(compareFigure, "vs_analytical_uz")
foreach(display, compareMetrics)
metricsPath = joinpath(DEV_OUTPUT_DIR[], DEV_STEP[], "vs_analytical_metrics.jld2")
jldsave(metricsPath; compareMetrics, observerDepths, optCacheFile, fdCacheFile)
println("metrics -> ", metricsPath)
