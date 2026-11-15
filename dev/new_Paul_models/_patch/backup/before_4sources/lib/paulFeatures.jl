# dev/new_Paul_models/lib/paulFeatures.jl
#
# Which waveform feature tells G1m20_sb11 from G1m20_sb21 most efficiently?
#
# For every source (north/south pole), component (vertical/horizontal), time
# window (PREM-predicted phase windows and fixed sliding windows) and station,
# we compute features of the ground-velocity segment:
#   lnRMS          ln rms amplitude (= ½ ln variance)
#   lnPeak         ln peak |v|
#   centroidFreq   spectral centroid (mHz)
#   bandwidth      spectral standard deviation around the centroid (mHz)
#   duration       t90 − t10 of the cumulative energy in the window (s)
#   energyCentroid energy-weighted mean time inside the window (s)
#   kurtosis       m4/m2² (impulsive vs emergent)
#   ccDelay        cross-correlation delay against the noise-free PREM trace (s)
#   decorrelation  1 − max normalised cc against the noise-free PREM trace
#   lnHV           ln(rms horizontal / rms vertical)   (component "HV")
#
# Discrimination: each model trace gets `nReal` realisations of band-limited
# Gaussian noise (same band as the source, rms = noiseLevel × rms of the PREM
# trace, the median over stations, for that source and component). Per station
#   z = (mean f_sb11 − mean f_sb21) / sqrt(var f_sb11 + var f_sb21)
# and a (feature, window, component, source) is scored by the fraction of
# stations with |z| > zThreshold ("detectable stations", ranked by their NUMBER,
# so that a phase seen at only a few stations does not win) and by rms(clamp(z, ±10)).

const FEATURE_NAMES = (:lnRMS, :lnPeak, :centroidFreq, :bandwidth, :duration,
    :energyCentroid, :kurtosis, :ccDelay, :decorrelation)

"Unit-rms Gaussian noise with the source spectrum exp(−(f T)²) (Gaussian kernel)."
function bandLimitedNoise(rng, n, dt, period)
    s = period / (π * sqrt(2)) / dt            # kernel σ in samples
    h = ceil(Int, 4s)
    kern = [exp(-0.5 * (k / s)^2) for k in -h:h]
    w = randn(rng, n + 2h)
    x = [sum(@view(w[i:i+2h]) .* kern) for i in 1:n]
    return x ./ sqrt(mean(abs2, x))
end

"Best normalised cc and delay (s, >0 = x later than ref) over ±nL samples."
function ccDelayIdx(ref, x, idx, nL, dt)
    a = @view ref[idx]
    na = sqrt(sum(abs2, a))
    na > 0 || return (NaN, NaN)
    best, kbest = -Inf, 0
    ccs = fill(-Inf, 2nL + 1)
    for (k, ℓ) in enumerate(-nL:nL)
        (first(idx) + ℓ >= 1 && last(idx) + ℓ <= length(x)) || continue
        b = @view x[idx .+ ℓ]
        nb = sqrt(sum(abs2, b))
        nb > 0 || continue
        ccs[k] = sum(a .* b) / (na * nb)
        ccs[k] > best && (best = ccs[k]; kbest = k)
    end
    kbest == 0 && return (NaN, NaN)
    shift = Float64(kbest - nL - 1)
    if 1 < kbest < length(ccs) && isfinite(ccs[kbest-1]) && isfinite(ccs[kbest+1])
        c1, c2, c3 = ccs[kbest-1], ccs[kbest], ccs[kbest+1]
        den = c1 - 2c2 + c3
        den < 0 && (shift += 0.5 * (c1 - c3) / den)
    end
    return (shift * dt, best)
end

"Feature vector (FEATURE_NAMES order) of trace x on sample range idx."
function segmentFeatures(x, ref, idx, dt, nL)
    seg = x[idx]
    n = length(seg)
    e = sum(abs2, seg)
    e > 0 || return fill(NaN, length(FEATURE_NAMES))
    pg = DSP.periodogram(seg .- mean(seg); fs = 1 / dt, nfft = max(nextpow(2, n), 256),
        window = DSP.hanning)
    f = DSP.freq(pg); P = DSP.power(pg)
    sel = f .> 0
    Ps = sum(P[sel])
    fc = sum(f[sel] .* P[sel]) / Ps
    bw = sqrt(sum((f[sel] .- fc) .^ 2 .* P[sel]) / Ps)
    c = cumsum(seg .^ 2) ./ e
    t10 = something(findfirst(>=(0.1), c), 1); t90 = something(findfirst(>=(0.9), c), n)
    tc = sum((0:n-1) .* seg .^ 2) / e * dt
    m2 = e / n
    kurt = mean(seg .^ 4) / m2^2
    δt, cc = ccDelayIdx(ref, x, idx, nL, dt)
    return [log(sqrt(m2)), log(maximum(abs, seg)), 1e3fc, 1e3bw, (t90 - t10) * dt, tc, kurt, δt, 1 - cc]
end

"""
    featureWindows(table, Δgrid, cfg) -> Vector of (label, kind, window(φ, pole) -> (t1, t2) | nothing)
"""
function featureWindows(table, Δgrid, cfg)
    wins = Any[]
    for ph in cfg.phases
        haskey(table, ph) || continue
        push!(wins, (label = ph, kind = :phase, window = (φ, pole) -> begin
            Δ = epicentralDistance(φ, pole)
            i = findfirst(==(round(Δ / step(Δgrid)) * step(Δgrid)), Δgrid)
            i === nothing && return nothing
            tp = table[ph]["time"][i]
            isfinite(tp) ? (tp - cfg.pre, tp + cfg.post) : nothing
        end))
    end
    for t0 in cfg.slideStart:cfg.slideStep:(cfg.duration - cfg.slideLength - cfg.maxLag)  # room for the cc lag
        t1 = t0 + cfg.slideLength
        push!(wins, (label = "$(round(Int, t0))-$(round(Int, t1))s", kind = :slide,
            window = (φ, pole) -> (t0, t1)))
    end
    return wins
end

"""
    computeFeatures(results, models, φs, windows, cfg) -> Dict

Key (pole, comp, model, windowLabel) => Array (nStation, nFeature, nReal+1);
slice [:, :, 1] is noise-free, [:, :, 2:end] the noisy realisations.
comp ∈ (:vert, :horiz); plus (pole, :HV, model, label) => (nStation, 1, nReal+1).
"""
function computeFeatures(results, models, φs, windows, cfg)
    out = Dict{Tuple{Symbol,Symbol,String,String},Array{Float64,3}}()
    nF = length(FEATURE_NAMES)
    for pole in SEISMO_POLES
        t = results[("PREM", pole)]["time"]; dt = t[2] - t[1]
        nL = round(Int, cfg.maxLag / dt)
        for comp in (:vert, :horiz)
            ref = seismoVelocity(results[("PREM", pole)][compKey(comp)], dt)
            # one noise level per source and component (not per trace: nodal stations
            # would get ~no noise and infinite z): noiseLevel × median PREM trace rms
            σ0 = cfg.noiseLevel * median([sqrt(mean(abs2, @view ref[:, j])) for j in eachindex(φs)])
            σ = fill(σ0, length(φs))
            for (mi, m) in enumerate(models)
                v = seismoVelocity(results[(m, pole)][compKey(comp)], dt)
                arrays = Dict(w.label => fill(NaN, length(φs), nF, cfg.nReal + 1) for w in windows)
                Threads.@threads for j in eachindex(φs)
                    # reproducible, independent per (pole, comp, model, station)
                    rng = Random.Xoshiro(hash((cfg.seed, pole, comp, mi, j)))
                    traces = Vector{Vector{Float64}}(undef, cfg.nReal + 1)
                    traces[1] = v[:, j]
                    for r in 1:cfg.nReal
                        traces[r+1] = v[:, j] .+ σ[j] .* bandLimitedNoise(rng, length(t), dt, cfg.period)
                    end
                    for w in windows
                        tw = w.window(φs[j], pole)
                        tw === nothing && continue
                        idx = findall(x -> tw[1] <= x <= tw[2], t)
                        (length(idx) >= 8 && first(idx) > nL && last(idx) + nL <= length(t)) || continue
                        A = arrays[w.label]
                        for r in 1:cfg.nReal+1
                            A[j, :, r] .= segmentFeatures(traces[r], @view(ref[:, j]), idx, dt, nL)
                        end
                    end
                end
                for (lab, A) in arrays
                    out[(pole, comp, m, lab)] = A
                end
                @info "features" pole comp model = m
            end
        end
        for m in models, w in windows
            h = out[(pole, :horiz, m, w.label)][:, 1:1, :]
            vv = out[(pole, :vert, m, w.label)][:, 1:1, :]
            out[(pole, :HV, m, w.label)] = h .- vv          # lnRMS_h − lnRMS_v
        end
    end
    return out
end

featureList(comp) = comp === :HV ? (:lnHV,) : FEATURE_NAMES

"""
    scoreFeatures(feat, models, windows, cfg) -> Vector of NamedTuples, sorted best first
"""
function scoreFeatures(feat, models, windows, cfg)
    m1, m2 = models
    rows = NamedTuple[]
    for pole in SEISMO_POLES, comp in (:vert, :horiz, :HV), w in windows
        A1 = feat[(pole, comp, m1, w.label)]; A2 = feat[(pole, comp, m2, w.label)]
        for (k, fname) in enumerate(featureList(comp))
            z = Float64[]; dfree = Float64[]
            for j in axes(A1, 1)
                a = filter(isfinite, A1[j, k, 2:end]); b = filter(isfinite, A2[j, k, 2:end])
                (length(a) >= 3 && length(b) >= 3) || continue
                s = sqrt(var(a) + var(b))
                s > 0 || continue
                push!(z, (mean(a) - mean(b)) / s)
                push!(dfree, A1[j, k, 1] - A2[j, k, 1])
            end
            isempty(z) && continue
            push!(rows, (source = pole, comp = comp, window = w.label, kind = w.kind, feature = fname,
                nStations = length(z), detectable = mean(abs.(z) .> cfg.zThreshold),
                nDetectable = count(abs.(z) .> cfg.zThreshold),
                rmsZ = sqrt(mean(abs2, clamp.(z, -10, 10))), medianAbsZ = median(abs.(z)),
                rmsDiffNoiseFree = sqrt(mean(abs2, filter(isfinite, dfree)))))
        end
    end
    sort!(rows; by = r -> (-r.nDetectable, -r.detectable, -r.rmsZ))
    return rows
end

compLabel(c) = c === :vert ? "vertical" : c === :horiz ? "horizontal" : "H/V"

"Heat map: detectable-station fraction, rows feature×component, columns windows."
function featureScoreFigure(rows, windows, pole, cfg)
    labs = [w.label for w in windows]
    rowKeys = [(c, f) for c in (:vert, :horiz) for f in FEATURE_NAMES]
    push!(rowKeys, (:HV, :lnHV))
    M = fill(NaN, length(labs), length(rowKeys))
    for r in rows
        r.source == pole || continue
        i = findfirst(==(r.window), labs); k = findfirst(==((r.comp, r.feature)), rowKeys)
        (i === nothing || k === nothing) || (M[i, k] = r.nDetectable)
    end
    fig = Figure(size = (90 * length(labs) + 400, 34 * length(rowKeys) + 250))
    ax = Axis(fig[1, 1]; xticks = (1:length(labs), labs), xticklabelrotation = π / 3,
        yticks = (1:length(rowKeys), ["$(compLabel(c)) · $f" for (c, f) in rowKeys]),
        yreversed = true, title = "$(pole) source: number of stations (of 180) where sb11 ≠ sb21 (|z| > $(cfg.zThreshold)), " *
            "noise $(round(Int, 100cfg.noiseLevel)) % of the PREM trace rms, $(cfg.nReal) realisations")
    hm = heatmap!(ax, 1:length(labs), 1:length(rowKeys), M; colormap = :viridis, colorrange = (0, 180),
        nan_color = :gray92)
    for i in eachindex(labs), k in eachindex(rowKeys)
        isfinite(M[i, k]) && text!(ax, i, k; text = string(round(Int, M[i, k])), align = (:center, :center),
            fontsize = 9, color = M[i, k] > 110 ? :black : :white)
    end
    nPh = count(w -> w.kind === :phase, windows)
    vlines!(ax, [nPh + 0.5]; color = :white, linewidth = 2)
    Colorbar(fig[1, 2], hm; label = "detectable stations")
    return fig
end

"Feature vs station angle for the best `n` combinations: noise-free value ± noise std."
function topFeaturesFigure(rows, feat, models, φs; n = 12, ncol = 4)
    top = rows[1:min(n, length(rows))]
    nrow = cld(length(top), ncol)
    fig = Figure(size = (520 * ncol, 330 * nrow + 60))
    colors = (:firebrick, :dodgerblue)
    for (i, r) in enumerate(top)
        ax = Axis(fig[cld(i, ncol), mod1(i, ncol)]; xlabel = "station angle φ (°)", ylabel = string(r.feature),
            xticks = 0:60:360, titlesize = 12,
            title = @sprintf("#%d %s · %s · %s src · %d/%d st.", i, r.window, compLabel(r.comp), r.source,
                r.nDetectable, r.nStations))
        k = findfirst(==(r.feature), collect(featureList(r.comp)))
        for (mi, m) in enumerate(models)
            A = feat[(r.source, r.comp, m, r.window)]
            ok = findall(j -> isfinite(A[j, k, 1]), axes(A, 1))
            isempty(ok) && continue
            μ = [A[j, k, 1] for j in ok]
            s = [std(filter(isfinite, A[j, k, 2:end])) for j in ok]
            band!(ax, φs[ok], μ .- s, μ .+ s; color = (colors[mi], 0.25))
            scatterlines!(ax, φs[ok], μ; color = colors[mi], markersize = 4, linewidth = 1,
                label = shortName(m))
        end
        i == 1 && axislegend(ax; position = :rt, labelsize = 10)
        xlims!(ax, 0, 360)
    end
    Label(fig[0, :], "Best discriminating features (noise-free value, band = ±1σ over the noise realisations)",
        fontsize = 16)
    return fig
end

"Best score per feature (any window/component/source): which kind of measurement works."
function featureRankingFigure(rows)
    fnames = vcat(collect(FEATURE_NAMES), [:lnHV])
    fig = Figure(size = (1500, 650))
    for (col, kind, title) in ((1, :phase, "phase windows"), (2, :slide, "sliding windows"))
        best = Float64[]; lab = String[]
        for f in fnames
            sel = [r for r in rows if r.feature == f && r.kind == kind]
            isempty(sel) && (push!(best, 0.0); push!(lab, string(f)); continue)
            b = sel[argmax([r.nDetectable + 1e-3r.rmsZ for r in sel])]
            push!(best, b.nDetectable)
            push!(lab, "$(f)\n$(b.window), $(compLabel(b.comp)), $(b.source)")
        end
        ax = Axis(fig[1, col]; yticks = (1:length(fnames), lab), yticklabelsize = 10,
            xlabel = "stations where sb11 ≠ sb21 (best window, of 180)", title = title, yreversed = true)
        barplot!(ax, 1:length(fnames), best; direction = :x, color = best, colormap = :viridis,
            colorrange = (0, 180))
        xlims!(ax, 0, 180)
    end
    return fig
end
