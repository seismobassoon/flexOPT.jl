function recordSectionFigure(results, names, pole, comp; tmax = SEISMO.duration,
        gain = SEISMO_SECTION.gain, clipAt = SEISMO_SECTION.clip, every = SEISMO_SECTION.every)
    r0 = results[("PREM", pole)]
    t = r0["time"]; dt = t[2] - t[1]
    φ = seismoReceivers.φ
    spacing = every * (φ[2] - φ[1])                 # degrees between plotted traces
    v = Dict(n => seismoVelocity(results[(n, pole)][compKey(comp)], dt) for n in names)
    scale = [max(maximum(abs, @view v["PREM"][:, j]), eps()) for j in eachindex(φ)]
    palette = (:black, :firebrick, :dodgerblue, :darkgreen)
    colors = Dict(n => palette[mod1(i, length(palette))] for (i, n) in enumerate(names))
    groups = vcat([[n] for n in names], [collect(names)])
    fig = Figure(size = (760 * length(groups), 2600))
    for (k, group) in enumerate(groups)
        ax = Axis(fig[1, k]; title = length(group) == 1 ? shortName(group[1]) : "overlay (" *
                join(shortName.(group), ", ") * ")", xlabel = "time (s)",
            ylabel = k == 1 ? "station angle φ (°)" : "", yticks = 0:10:360, titlesize = 16)
        k > 1 && hideydecorations!(ax; grid = false, ticks = false, ticklabels = false)
        overlayCurves!(ax, pole, SEISMO_DEPTH_KM; labels = k == 1, tmax = tmax)
        for n in group
            xs = Float64[]; ys = Float64[]
            for j in 1:every:length(φ)
                w = clamp.(gain .* v[n][:, j] ./ scale[j], -clipAt, clipAt)
                append!(xs, t); push!(xs, NaN)
                append!(ys, φ[j] .+ 0.5spacing .* w); push!(ys, NaN)
            end
            lines!(ax, xs, ys; color = (colors[n], length(group) == 1 ? 1.0 : 0.8), linewidth = 0.55)
        end
        xlims!(ax, 0, tmax); ylims!(ax, -spacing, 360)
    end
    Label(fig[0, :], "$(compName(comp)) ground velocity, source at the $(pole) pole (2000 km). " *
        "Each trace / max|PREM trace| × $(gain), clipped at ±$(clipAt); coloured lines: PREM ray theory " *
        "(black PREM, red sb11, blue sb21)", fontsize = 18)
    return fig
end
