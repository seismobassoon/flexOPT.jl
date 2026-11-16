# dev/new_Paul_models/lib/paulSeismoFigures.jl
# Figures and phase-window measurements for the synthetics (lib/paulSeismo.jl).
# Travel times / ray paths: PREM with dev/raytheory (τ–p + ray tracer).
# All analysis is on ground VELOCITY (∂t of the stored displacement): the
# smooth-step source leaves static offsets in displacement.

SOURCE_COLORS = Dict(:north => :gold, :east => :limegreen, :south => :magenta, :west => :cyan)

seismoVelocity(u, dt) = vcat(zeros(1, size(u, 2)), diff(Float64.(u); dims = 1) ./ dt)

"Last leg P -> vertical, last leg S -> horizontal (radial SV)."
phaseComponent(ph) = uppercase(last(ph)) == 'P' ? :vert : :horiz
"First leg P -> δlnVp background, S -> δlnVs."
phaseFirstMode(ph) = uppercase(first(ph)) == 'P' ? :P : :S
compKey(comp) = comp === :vert ? "uVert" : "uHoriz"
compName(comp) = comp === :vert ? "u_vertical" : "u_horizontal"
shortName(name) = replace(name, "G1m20_" => "")

@isdefined(_TTCURVES) || (_TTCURVES = Dict{Tuple{String,Float64},Any}())
"PREM T(Δ) branches (memoised); empty if the phase does not exist / fails."
function premCurves(phase, depthKm)
    get!(_TTCURVES, (phase, depthKm)) do
        try
            travelTimeCurve(SEISMO_RAY_PREM, phase, depthKm)
        catch err
            @warn "travelTimeCurve failed" phase exception = err
            Any[]
        end
    end
end

"Station angles φ ∈ [0, 360) of epicentral distance Δ, both sides of the source."
stationAngles(Δ, pole) = (mod.(sourceAngle(pole) .+ Δ, 360), mod.(sourceAngle(pole) .- Δ, 360))

"Break a curve where φ wraps through 0/360 (NaN separators) so lines do not jump."
function wrapBreak(T, φ)
    x = Float64[]; y = Float64[]
    for i in eachindex(φ)
        i > 1 && abs(φ[i] - φ[i-1]) > 180 && (push!(x, NaN); push!(y, NaN))
        push!(x, T[i]); push!(y, φ[i])
    end
    return x, y
end

function overlayCurves!(ax, pole, depthKm; phases = SEISMO_SECTION_PHASES, labels = true, tmax = Inf)
    colors = Makie.wong_colors()
    for (k, ph) in enumerate(phases)
        col = colors[mod1(k, length(colors))]
        for (Δs, Ts, _) in premCurves(ph, depthKm)
            for φ in stationAngles(Δs, pole)
                lines!(ax, wrapBreak(Ts, φ)...; color = (col, 0.8), linewidth = 0.8)
            end
        end
        if labels
            br = premCurves(ph, depthKm)
            isempty(br) && continue
            Δs, Ts, _ = br[1]
            j = clamp(length(Δs) ÷ 3, 1, length(Δs))
            Ts[j] < tmax && text!(ax, Ts[j], stationAngles(Δs[j], pole)[1]; text = ph,
                color = col, fontsize = 11, offset = (3, 0))
        end
    end
end

function seismoModelFigure(anoms, disk, sources, planetRadius, cmbRadius)
    names = sort(collect(keys(anoms)))
    fig = Figure(size = (1500, 470 * length(names)))
    θ = range(0, 2π; length = 361)
    for (i, name) in enumerate(names), (j, (key, lab)) in enumerate(((:δvp, "δln Vp"), (:δvs, "δln Vs"), (:δρ, "δln ρ")))
        A = 100 .* getproperty(anoms[name], key)
        A = ifelse.(disk.mantle, A, NaN)
        gl = fig[i, j] = GridLayout()
        ax = Axis(gl[1, 1]; aspect = DataAspect(), title = "$(shortName(name)): $lab (%)", titlesize = 14)
        hidedecorations!(ax); hidespines!(ax)
        lim = robustLimits(A; symmetric = true, hi = 0.99)
        hm = heatmap!(ax, disk.xs ./ 1e3, disk.zs ./ 1e3, A; colormap = Reverse(:vik), colorrange = lim,
            nan_color = :transparent)
        for r in (cmbRadius, planetRadius)
            lines!(ax, r / 1e3 .* cos.(θ), r / 1e3 .* sin.(θ); color = :black, linewidth = 0.6)
        end
        for (pole, s) in sources
            scatter!(ax, [s.position[1] / 1e3], [s.position[2] / 1e3]; marker = :star5,
                markersize = 18, color = SOURCE_COLORS[pole], strokecolor = :black, strokewidth = 1)
        end
        for φ in 0:90:270
            text!(ax, 1.04planetRadius / 1e3 * sind(φ), 1.04planetRadius / 1e3 * cosd(φ);
                text = "$(φ)°", align = (:center, :center), fontsize = 11)
        end
        Colorbar(gl[1, 2], hm; height = Relative(0.7), label = "%")
    end
    Label(fig[0, :], "PREM × (1 + δ): lateral anomalies of Paul's Pvel, Svel, rho (mantle only, clipped). " *
        "Stars: sources at $(round(Int, SEISMO.sourceDepth / 1e3)) km (" *
        join(["$(SOURCE_COLORS[p]) $p" for p in SEISMO_POLES], ", ") * "). Red = slow / light.", fontsize = 15)
    return fig
end

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
    Label(fig[0, :], "$(compName(comp)) ground velocity, $(pole) source (φ = $(round(Int, sourceAngle(pole)))°, $(round(Int, SEISMO.sourceDepth / 1e3)) km). " *
        "Each trace / max|PREM trace| × $(gain), clipped at ±$(clipAt); coloured lines: PREM ray theory " *
        "(black PREM, red sb11, blue sb21)", fontsize = 18)
    return fig
end

function snapshotFigure(results, names, pole, disk)
    times = sort(collect(SEISMO.snapshotTimes))
    fig = Figure(size = (330 * length(times), 330 * length(names) + 50))
    for (i, n) in enumerate(names), (j, ts) in enumerate(times)
        snaps = results[(n, pole)]["snapshots"]
        haskey(snaps, ts) || continue
        ref = results[("PREM", pole)]["snapshots"][ts]
        top = log10(max(maximum(ref), eps(Float32)))
        A = log10.(max.(snaps[ts], eps(Float32)))
        A = ifelse.(disk.material, A, NaN)
        ax = Axis(fig[i, j]; aspect = DataAspect(), title = "$(shortName(n)), t = $(Int(ts)) s", titlesize = 12)
        hidedecorations!(ax); hidespines!(ax)
        heatmap!(ax, disk.xs ./ 1e3, disk.zs ./ 1e3, A; colormap = :inferno,
            colorrange = (top - 3, top), nan_color = :transparent)
    end
    Label(fig[0, :], "log₁₀ |u|, $(pole) source (3 decades below the PREM maximum at each time)", fontsize = 15)
    return fig
end

# ---- phase windows -----------------------------------------------------------

"""
PREM first arrival (time, p, startDir) of each phase on the Δ grid (NaN if none).
"""
function taupTable(m, phases, depthKm, Δgrid)
    tab = Dict{String,Any}()
    for ph in phases
        T = fill(NaN, length(Δgrid)); P = fill(NaN, length(Δgrid)); dirs = fill(:none, length(Δgrid))
        for (i, Δ) in enumerate(Δgrid)
            (0.5 < Δ < 179.5) || continue
            arr = try
                taupArrivals(m, ph, depthKm, Δ)
            catch err
                Arrival[]
            end
            isempty(arr) && continue
            T[i] = arr[1].time; P[i] = arr[1].p; dirs[i] = arr[1].startDir
        end
        tab[ph] = Dict("time" => T, "p" => P, "dir" => String.(dirs))
        @printf("  %-6s found at %3d of %d distances\n", ph, count(isfinite, T), length(Δgrid))
    end
    return tab
end

"Cross-correlation delay (model − PREM, s, >0 = model late), ln amplitude ratio, cc."
function measureWindow(v0, v1, t, tp, w)
    dt = t[2] - t[1]
    idx = findall(x -> tp - w.pre <= x <= tp + w.post, t)
    length(idx) < 5 && return (δt = NaN, dlnA = NaN, cc = NaN)
    nL = round(Int, w.maxLag / dt)
    a = v0[idx]
    na = sqrt(sum(abs2, a))
    na > 0 || return (δt = NaN, dlnA = NaN, cc = NaN)
    ccs = fill(-Inf, 2nL + 1)
    for (k, ℓ) in enumerate(-nL:nL)
        j = idx .+ ℓ
        (first(j) >= 1 && last(j) <= length(v1)) || continue
        b = @view v1[j]
        nb = sqrt(sum(abs2, b))
        nb > 0 && (ccs[k] = sum(a .* b) / (na * nb))
    end
    k = argmax(ccs)
    isfinite(ccs[k]) || return (δt = NaN, dlnA = NaN, cc = NaN)
    shift = Float64(k - nL - 1)
    if 1 < k < length(ccs) && all(isfinite, ccs[k-1:k+1])   # parabolic refinement
        c1, c2, c3 = ccs[k-1], ccs[k], ccs[k+1]
        den = c1 - 2c2 + c3
        den < 0 && (shift += 0.5 * (c1 - c3) / den)
    end
    dlnA = log(sqrt(sum(abs2, v1[idx])) / na)
    return (δt = shift * dt, dlnA = dlnA, cc = ccs[k])
end

"""
All measurements: Vector of NamedTuples (phase, pole, model, φ, Δ, tPred, δt, dlnA, cc, clean).
`clean` = no other phase of SEISMO_PHASES predicted inside the window.
"""
function measurePhases(results, names, table, Δgrid, φs, w)
    out = NamedTuple[]
    for pole in SEISMO_POLES
        t = results[("PREM", pole)]["time"]; dt = t[2] - t[1]
        vel = Dict((n, c) => seismoVelocity(results[(n, pole)][compKey(c)], dt)
            for n in names, c in (:vert, :horiz))
        for ph in keys(table), (j, φ) in enumerate(φs)
            Δ = epicentralDistance(φ, pole)
            i = findfirst(==(round(Δ / step(Δgrid)) * step(Δgrid)), Δgrid)
            i === nothing && continue
            tp = table[ph]["time"][i]
            isfinite(tp) || continue
            clean = !any(q -> q != ph && isfinite(table[q]["time"][i]) &&
                tp - w.pre <= table[q]["time"][i] <= tp + w.post, keys(table))
            c = phaseComponent(ph)
            v0 = @view vel[("PREM", c)][:, j]
            for n in names[2:end]
                m = measureWindow(v0, @view(vel[(n, c)][:, j]), t, tp, w)
                push!(out, (phase = ph, pole = pole, model = n, φ = φ, Δ = Δ, tPred = tp,
                    δt = m.δt, dlnA = m.dlnA, cc = m.cc, clean = clean))
            end
        end
    end
    return out
end

"PREM ray paths of `phase` to Δ = 10°, 20°, ..., as (x, z) in km for a pole source, both sides."
function phaseRayPaths(m, phase, depthKm, pole; Δs = 10.0:10.0:170.0)
    paths = Vector{Tuple{Vector{Float64},Vector{Float64}}}()
    ph = parsePhase(phase)
    φs = deg2rad(sourceAngle(pole))
    for Δ in Δs
        arr = try taupArrivals(m, phase, depthKm, Δ) catch; Arrival[] end
        isempty(arr) && continue
        st = sphericalStartP(m, depthKm, arr[1].p, ph.mode, arr[1].startDir)
        st === nothing && continue
        ray = try traceRay(m, ph, st...; dT = 2.0) catch; nothing end
        ray === nothing && continue
        r = [x[1] for x in ray.path]; θ = [x[2] for x in ray.path]
        for side in (1.0, -1.0)
            ang = φs .+ side .* θ                  # station angle along the ray
            push!(paths, (r .* sin.(ang), r .* cos.(ang)))
        end
    end
    return paths
end

function phaseFigure(phase, results, names, anoms, disk, table, Δgrid, meas, planetRadius, cmbRadius)
    comp = phaseComponent(phase)
    key = phaseFirstMode(phase) === :P ? :δvp : :δvs
    w = SEISMO_WINDOW
    models = names[2:end]
    colors = Dict("PREM" => :black, models[1] => :firebrick)
    length(models) > 1 && (colors[models[2]] = :dodgerblue)
    fig = Figure(size = (2700, 560 * length(SEISMO_POLES) + 60))
    θc = range(0, 2π; length = 361)
    for (row, pole) in enumerate(SEISMO_POLES)
        paths = phaseRayPaths(SEISMO_RAY_PREM, phase, SEISMO_DEPTH_KM, pole)
        # ray maps
        for (col, n) in enumerate(models)
            A = ifelse.(disk.mantle, 100 .* getproperty(anoms[n], key), NaN)
            ax = Axis(fig[row, col]; aspect = DataAspect(),
                title = "$(shortName(n)) $(key === :δvp ? "δlnVp" : "δlnVs"), PREM $phase rays ($(pole) source)",
                titlesize = 12)
            hidedecorations!(ax); hidespines!(ax)
            heatmap!(ax, disk.xs ./ 1e3, disk.zs ./ 1e3, A; colormap = Reverse(:vik),
                colorrange = robustLimits(A; symmetric = true, hi = 0.99), nan_color = :transparent)
            for r in (1221.5, cmbRadius / 1e3, planetRadius / 1e3)
                lines!(ax, r .* cos.(θc), r .* sin.(θc); color = :gray30, linewidth = 0.6)
            end
            for (x, z) in paths
                lines!(ax, x, z; color = (:black, 0.7), linewidth = 0.7)
            end
        end
        # gather, aligned on the PREM time
        t = results[("PREM", pole)]["time"]; dt = t[2] - t[1]
        ax = Axis(fig[row, 3]; title = "$phase, $(compName(comp)) velocity, aligned on PREM time",
            xlabel = "t − t_PREM (s)", ylabel = "station angle φ (°)", yticks = 0:30:360)
        for n in names
            V = seismoVelocity(results[(n, pole)][compKey(comp)], dt)
            V0 = seismoVelocity(results[("PREM", pole)][compKey(comp)], dt)
            xsAll = Float64[]; ysAll = Float64[]
            for (j, φ) in enumerate(seismoReceivers.φ)
                Δ = epicentralDistance(φ, pole)
                i = findfirst(==(round(Δ / step(Δgrid)) * step(Δgrid)), Δgrid)
                i === nothing && continue
                tp = table[phase]["time"][i]
                isfinite(tp) || continue
                idx = findall(x -> tp - w.pre <= x <= tp + w.post, t)
                isempty(idx) && continue
                s = max(maximum(abs, @view V0[idx, j]), eps())
                append!(xsAll, t[idx] .- tp); push!(xsAll, NaN)
                append!(ysAll, φ .+ 3.0 .* clamp.(V[idx, j] ./ s, -1.3, 1.3)); push!(ysAll, NaN)
            end
            lines!(ax, xsAll, ysAll; color = (colors[n], n == "PREM" ? 0.9 : 0.75),
                linewidth = 0.6, label = shortName(n))
        end
        vlines!(ax, [0.0]; color = :gray, linestyle = :dash)
        xlims!(ax, -w.pre, w.post); ylims!(ax, -2, 360)
        row == 1 && axislegend(ax; position = :rt, framevisible = true, labelsize = 10)
        # δt and dlnA versus φ
        for (col, field, lab) in ((4, :δt, "δt = t_model − t_PREM (s)"), (5, :dlnA, "ln(A_model / A_PREM)"))
            ax2 = Axis(fig[row, col]; xlabel = "station angle φ (°)", ylabel = lab, xticks = 0:30:360,
                title = "$phase, $(pole) source")
            for n in models
                sel = [x for x in meas if x.phase == phase && x.pole == pole && x.model == n && isfinite(getproperty(x, field))]
                isempty(sel) && continue
                φc = [x.φ for x in sel if x.clean]; yc = [getproperty(x, field) for x in sel if x.clean]
                φd = [x.φ for x in sel if !x.clean]; yd = [getproperty(x, field) for x in sel if !x.clean]
                isempty(φc) || scatter!(ax2, Float64.(φc), Float64.(yc); color = colors[n],
                    markersize = 6, label = shortName(n))
                isempty(φd) || scatter!(ax2, Float64.(φd), Float64.(yd); color = :transparent,
                    strokecolor = colors[n], strokewidth = 1, markersize = 6)
            end
            hlines!(ax2, [0.0]; color = :gray, linestyle = :dash)
            xlims!(ax2, 0, 360)
            row == 1 && col == 4 && axislegend(ax2; position = :rt, framevisible = true, labelsize = 10)
        end
    end
    Label(fig[0, :], "$phase: sources at $(round(Int, SEISMO.sourceDepth / 1e3)) km depth. Hollow markers: another phase predicted inside the " *
        "window [−$(Int(w.pre)), +$(Int(w.post))] s (measurement contaminated). φ: 0 north, 90 east, 180 south, 270 west.",
        fontsize = 15)
    colsize!(fig.layout, 1, Relative(0.18)); colsize!(fig.layout, 2, Relative(0.18))
    return fig
end

"δt of every phase (rows) vs φ for both sources (columns): the overview."
function phaseSummaryFigure(meas, phases, models)
    fig = Figure(size = (850 * length(SEISMO_POLES), 210 * length(phases) + 80))
    colors = Dict(models[1] => :firebrick)
    length(models) > 1 && (colors[models[2]] = :dodgerblue)
    for (i, ph) in enumerate(phases), (j, pole) in enumerate(SEISMO_POLES)
        ax = Axis(fig[i, j]; ylabel = "$ph δt (s)", xticks = 0:30:360,
            xlabel = i == length(phases) ? "station angle φ (°)" : "",
            title = i == 1 ? "$(pole) source" : "")
        for n in models
            sel = [x for x in meas if x.phase == ph && x.pole == pole && x.model == n && isfinite(x.δt)]
            isempty(sel) && continue
            c = [x for x in sel if x.clean]; d = [x for x in sel if !x.clean]
            isempty(c) || scatter!(ax, Float64[x.φ for x in c], Float64[x.δt for x in c]; markersize = 4,
                color = colors[n], label = shortName(n))
            isempty(d) || scatter!(ax, Float64[x.φ for x in d], Float64[x.δt for x in d]; markersize = 4,
                color = (colors[n], 0.25))
        end
        hlines!(ax, [0.0]; color = :gray, linestyle = :dash)
        xlims!(ax, 0, 360)
        i == 1 && j == length(SEISMO_POLES) && axislegend(ax; position = :rt, labelsize = 10)
    end
    Label(fig[0, :], "Cross-correlation delays model − PREM per phase (faded: window contaminated by another phase)",
        fontsize = 15)
    return fig
end

"""
    anrSectionFigure(results, names, pole, comp)

Record section in the style of the ANR 2027 figure
(data/earth2D_stagyy/prototype_python/make_figure.py): for each heterogeneous
model and each side of the source, grey = model trace, red = δu = model − PREM,
both divided by max|model trace| for t < SEISMO_ANR.tRef and clipped at
±SEISMO_ANR.clip; traces every SEISMO_ANR.everyDeg of epicentral distance.
"""
function anrSectionFigure(results, names, pole, comp; cfg = SEISMO_ANR)
    t = results[("PREM", pole)]["time"]; dt = t[2] - t[1]
    φ = seismoReceivers.φ
    key = compKey(comp)
    v0 = seismoVelocity(results[("PREM", pole)][key], dt)
    early = t .< cfg.tRef
    sp = 0.48 * cfg.everyDeg
    models = names[2:end]
    fig = Figure(size = (620 * 2 * length(models), 1150))
    sides = ((1, "φ = $(round(Int, sourceAngle(pole))) + Δ"), (-1, "φ = $(round(Int, sourceAngle(pole))) − Δ"))
    col = 0
    for m in models
        v = seismoVelocity(results[(m, pole)][key], dt)
        for (sgn, sideLabel) in sides
            col += 1
            ax = Axis(fig[1, col]; xlabel = "time (s)", ylabel = col == 1 ? "epicentral distance Δ (°)" : "",
                title = "$(shortName(m)), $(compName(comp)), $sideLabel", titlesize = 14,
                xgridvisible = false, ygridvisible = false, rightspinevisible = false, topspinevisible = false,
                yticks = 0:30:180)
            xg = Float64[]; yg = Float64[]; xr = Float64[]; yr = Float64[]
            for Δ in 0.0:cfg.everyDeg:180.0
                φt = mod(sourceAngle(pole) + sgn * Δ, 360)
                j = findfirst(x -> isapprox(x, φt; atol = 1e-6), φ)
                j === nothing && continue
                ref = maximum(abs, @view v[early, j])
                ref > 0 || continue
                append!(xg, t); push!(xg, NaN)
                append!(yg, Δ .+ sp .* clamp.(v[:, j] ./ ref, -cfg.clip, cfg.clip)); push!(yg, NaN)
                append!(xr, t); push!(xr, NaN)
                append!(yr, Δ .+ sp .* clamp.((v[:, j] .- v0[:, j]) ./ ref, -cfg.clip, cfg.clip)); push!(yr, NaN)
            end
            lines!(ax, xg, yg; color = (:gray60, 1.0), linewidth = 0.8, label = "$(shortName(m))")
            lines!(ax, xr, yr; color = :firebrick, linewidth = 0.7, label = "δu = $(shortName(m)) − PREM")
            col == 1 && axislegend(ax; position = :rb, labelsize = 11, framevisible = true)
            xlims!(ax, 0, cfg.tmax); ylims!(ax, -4, 184)
        end
    end
    Label(fig[0, :], "$(compName(comp)) ground velocity, $(pole) source ($(round(Int, SEISMO.sourceDepth / 1e3)) km, " *
        "$(SEISMO.dip)° $(SEISMO.slip), T = $(round(Int, SEISMO.period)) s); each trace / max|model trace| for " *
        "t < $(round(Int, cfg.tRef)) s, clipped ±$(cfg.clip)", fontsize = 16)
    return fig
end
