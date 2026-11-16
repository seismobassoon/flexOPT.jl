# dev/new_Paul_models/lib/paulNeutrino.jl
#
# Atmospheric neutrinos through the same 2-D Earth as the synthetics, seen by
# detectors at the surface (2.5 km deep) at station angles SEISMO_SOURCE_ANGLES
# (north, east, south, west).
#
#   nₑ = ρ Yₑ   (electronDensity, mol e⁻ cm⁻³)
#     PREM  : ρ_PREM (planet1D, LID above the Moho, as the seismic runs),
#             Yₑ = NU.yeCore in the core, NU.yeMantle in the mantle;
#     models: ρ_PREM (1 + δlnρ)  ×  Yₑ_mantle (1 + δlnYₑ),
#             δlnYₑ = lateral anomaly of Yₑ(bs) = Yₑ,hz + bs (Yₑ,MORB − Yₑ,hz)
#             (bs = mass fraction, end-members from Paul's NCFMAS compositions),
#     so that, like for the seismic models, PREM is the spherical average.
#
#   GeoPoints.creationPaths only shoots from a detector at the north pole
#   towards −x. Each detector/side is therefore obtained by ROTATING (and, for
#   the other side, mirroring) the nₑ field so that the detector sits at the
#   north pole; src/ is not touched. Side "−": rays towards decreasing φ,
#   side "+": towards increasing φ (φ: 0 north, 90 east = +x).
#
#   Then pathsFromElectronDensity → Pνν (ν and ν̄) → interacting_events
#   (Daemonflux + CC cross sections, perfect detector), as in
#   scripts/earth2D_StagYY_neutrino_seismo.jl.

"Map a point of the rotated frame (detector at the north pole) to the original frame."
function nuFrameToOriginal(x, z, φd, side)
    xm = side * x
    s, c = sind(φd), cosd(φd)
    return (xm * c + z * s, -xm * s + z * c)
end

"nₑ resampled so that a detector at angle φd sits at the north pole (side = ±1)."
function rotatedField(F, xs, zs, φd, side)
    itp = LinearInterpolation((xs, zs), F; extrapolation_bc = 0.0)
    out = similar(F, Float64)
    Threads.@threads for k in eachindex(zs)
        for i in eachindex(xs)
            xo, zo = nuFrameToOriginal(xs[i], zs[k], φd, side)
            out[i, k] = (abs(xo) <= xs[end] && abs(zo) <= zs[end]) ? itp(xo, zo) : 0.0
        end
    end
    return out
end

"Yₑ on the disk: PREM layers (core / mantle) and the models' lateral δlnYₑ."
function nuElectronFraction(disk, cmbRadius; yeCore, yeMantle, δlnYe = nothing)
    Ye = fill(0.0, size(disk.radius))
    Ye[disk.material .& (disk.radius .< cmbRadius)] .= yeCore
    Ye[disk.mantle] .= yeMantle
    δlnYe === nothing || (Ye .*= (1 .+ δlnYe))
    return Ye
end

"Lateral anomaly of Yₑ(bs) on the disk (mantle only)."
function yeLateralAnomaly(bsFile, disk, planetRadius, cmbRadius; Ye_hz, Ye_morb, seamBlendDeg = 0.0)
    s = getCartesianField(bsFile, disk.X, disk.Z;
        interpolation_method = :polar, target_cmb_radius = cmbRadius,
        target_surface_radius = planetRadius, clamp_to_surface = true,
        seam_blend_angle = deg2rad(seamBlendDeg))
    ye = Ye_hz .+ clamp.(s.field, 0, 1) .* (Ye_morb - Ye_hz)
    yeAv = Ye_hz .+ clamp.(s.avField, 0, 1) .* (Ye_morb - Ye_hz)
    return ifelse.(disk.mantle, ye ./ yeAv .- 1, 0.0)
end

"""
    runNeutrinoDetector(ne, disk, φd, side, cfg, osc, fluxTable, xsecTable) -> Dict

Pee, Pμe (ν) and the combined νe / νμ CC events (energy × cosθ) for one
detector side; cfg = NU (binning, Δbaseline, detectorDepth, exposure).
"""
function runNeutrinoDetector(ne, disk, planetRadius, φd, side, cfg, osc, fluxTable, xsecTable)
    nRot = rotatedField(ne, disk.xs, disk.zs, φd, side)
    detectorLocal = SVector(0.0, planetRadius - cfg.detectorDepth)
    sampling = creationPaths(nRot, disk.xs, disk.zs, detectorLocal, cfg.binning;
        Δbaseline = cfg.Δbaseline, material_mask = disk.material, modifiedLongitude2Disk = 0.0,
        center = (0.0, 0.0), boundary_search_step = (disk.xs[2] - disk.xs[1]) / 2, boundary_tolerance = 10.0)
    conversion = pathsFromElectronDensity(sampling; reference_z_over_a = 0.5)
    ν = Pνν(osc, cfg.binning.energies, conversion.paths; zoa = conversion.reference_z_over_a, roundU_and_H = false)
    aν = Pνν(osc, cfg.binning.energies, conversion.paths; zoa = conversion.reference_z_over_a, anti = true,
        roundU_and_H = false)
    det = PerfectDetector(mass_kg = cfg.massKg, livetime_s = cfg.years * 365.25 * 24 * 3600)
    ev = interacting_events(fluxTable, xsecTable, ν, aν; detector = det)
    P = Array(parent(ν))
    # ray end points in the ORIGINAL frame (for maps)
    rays = [(nuFrameToOriginal(p.source..., φd, side), nuFrameToOriginal(p.detector..., φd, side), p.cosθ)
        for p in sampling.profiles]
    return Dict{String,Any}("Pee" => P[:, :, 1, 1], "Pmue" => P[:, :, 2, 1],
        "Ne" => Array(ev.combined.νe), "Nmu" => Array(ev.combined.νμ),
        "energies" => collect(cfg.binning.energies), "cosθ" => collect(cfg.binning.cosθgrid),
        "rays" => [(collect(a), collect(b), c) for (a, b, c) in rays])
end

# ---- zenith classes and energy bands -----------------------------------------
"Deepest region crossed by a straight path of zenith cosine c (surface detector)."
function zenithClass(c, R, cmb, icb, d670)
    rmin = R * sqrt(max(0.0, 1 - c^2))
    rmin < icb && return :innerCore
    rmin < cmb && return :outerCore
    rmin < R - d670 && return :lowerMantle
    return :upperMantle
end
const NU_CLASSES = (:upperMantle, :lowerMantle, :outerCore, :innerCore)

"Sum of events over a zenith class and an energy band."
function binSum(N, E, cosθ, classOf, cls, band)
    s = 0.0
    for (j, c) in enumerate(cosθ), (i, e) in enumerate(E)
        (classOf[j] === cls && band[1] <= e < band[2]) && (s += N[i, j])
    end
    return s
end

"Δχ² between two event tables over the (E, cosθ) bins of a class/band (Poisson, perfect detector)."
function binChi2(N1, N2, E, cosθ, classOf, cls, band)
    s = 0.0
    for (j, c) in enumerate(cosθ), (i, e) in enumerate(E)
        (classOf[j] === cls && band[1] <= e < band[2]) || continue
        m = 0.5 * (N1[i, j] + N2[i, j])
        m > 0 && (s += (N1[i, j] - N2[i, j])^2 / m)
    end
    return s
end

"""
    neutrinoTable(res, models, detectors, sides, cfg, R, cmb, icb) -> Vector of NamedTuples

For every detector side, channel (νe, νμ), zenith class and energy band:
N for each model, δN/N against PREM, and Δχ²(model₁ vs model₂).
"""
function neutrinoTable(res, models, detectors, sides, cfg, R, cmb, icb)
    rows = NamedTuple[]
    r0 = res[("PREM", first(detectors), first(sides))]
    E, cθ = r0["energies"], r0["cosθ"]
    classOf = [zenithClass(c, R, cmb, icb, cfg.d670) for c in cθ]
    m1, m2 = models[2], models[3]
    for det in detectors, side in sides, (ch, key) in (("νe", "Ne"), ("νμ", "Nmu")),
            cls in NU_CLASSES, band in cfg.energyBands
        N = Dict(m => binSum(res[(m, det, side)][key], E, cθ, classOf, cls, band) for m in models)
        N[models[1]] > 0 || continue
        push!(rows, (detector = det, side = side, channel = ch, class = cls, band = band,
            N_PREM = N[models[1]], N_1 = N[m1], N_2 = N[m2],
            dNN_1 = N[m1] / N[models[1]] - 1, dNN_2 = N[m2] / N[models[1]] - 1,
            chi2_12 = binChi2(res[(m1, det, side)][key], res[(m2, det, side)][key], E, cθ, classOf, cls, band),
            chi2_1P = binChi2(res[(m1, det, side)][key], res[(models[1], det, side)][key], E, cθ, classOf, cls, band),
            chi2_2P = binChi2(res[(m2, det, side)][key], res[(models[1], det, side)][key], E, cθ, classOf, cls, band)))
    end
    return rows
end

bandLabel(b) = "$(b[1])–$(b[2]) GeV"
classLabel(c) = c === :upperMantle ? "upper mantle" : c === :lowerMantle ? "lower mantle" :
    c === :outerCore ? "outer core" : "inner core"

"""
Bar plots in the style of the seismic phase summary: x = zenith class × energy
band; top: rms over the detector sides of δN/N (sb11 red, sb21 blue, ×ratio);
bottom: √Δχ²(sb11 vs sb21) per detector (both sides summed), i.e. the number
of σ at which the detector separates the two models for the given exposure.
"""
function neutrinoBarFigure(rows, channel, models, detectors, cfg)
    groups = [(c, b) for c in NU_CLASSES for b in cfg.energyBands]
    sel(c, b) = [r for r in rows if r.channel == channel && r.class === c && r.band == b]
    x = 1:length(groups)
    rms1 = [isempty(sel(g...)) ? NaN : 100 * sqrt(mean(abs2, [r.dNN_1 for r in sel(g...)])) for g in groups]
    rms2 = [isempty(sel(g...)) ? NaN : 100 * sqrt(mean(abs2, [r.dNN_2 for r in sel(g...)])) for g in groups]
    fig = Figure(size = (max(1200, 95 * length(groups)), 950))
    labels = ["$(classLabel(c))\n$(bandLabel(b))" for (c, b) in groups]
    ax1 = Axis(fig[1, 1]; xticks = (x, labels), xticklabelsize = 9, xticklabelrotation = π / 4,
        ylabel = "rms δN/N vs PREM (%)",
        title = "$channel CC events: deviation from PREM over the $(length(detectors)) detectors × 2 sides")
    w = 0.38
    barplot!(ax1, x .- w / 2, rms1; width = w, color = :firebrick, label = shortName(models[2]))
    barplot!(ax1, x .+ w / 2, rms2; width = w, color = :dodgerblue, label = shortName(models[3]))
    for (i, (a, b)) in enumerate(zip(rms1, rms2))
        (isfinite(a) && isfinite(b) && b > 0) && text!(ax1, i, max(a, b) * 1.03; text = @sprintf("×%.1f", a / b),
            align = (:center, :bottom), fontsize = 10)
    end
    axislegend(ax1; position = :lt)
    ax2 = Axis(fig[2, 1]; xticks = (x, labels), xticklabelsize = 9, xticklabelrotation = π / 4,
        ylabel = "√Δχ²  sb11 vs sb21 (σ)", yscale = log10,
        title = @sprintf("separation of the two models per detector (%.0f kt × %.0f yr, perfect detector, both sides)",
            cfg.massKg / 1e6, cfg.years))
    nd = length(detectors); wd = 0.8 / nd
    palette = [SOURCE_COLORS[d] for d in detectors]
    for (k, det) in enumerate(detectors)
        v = [sqrt(sum([r.chi2_12 for r in sel(g...) if r.detector == det]; init = 0.0)) for g in groups]
        barplot!(ax2, x .+ (k - (nd + 1) / 2) * wd, max.(v, 1e-3); width = wd, color = palette[k],
            strokecolor = :black, strokewidth = 0.5, label = string(det))
    end
    hlines!(ax2, [3.0]; color = :gray, linestyle = :dash)
    axislegend(ax2; position = :lt, nbanks = 2)
    return fig
end

"Maps of δlnnₑ for the models with the 4 detectors and a few rays per side."
function neutrinoMapFigure(neFields, disk, res, models, detectors, sides, planetRadius, cmbRadius)
    fig = Figure(size = (650 * (length(models) - 1), 650))
    θ = range(0, 2π; length = 361)
    ne0 = neFields[models[1]]
    for (k, m) in enumerate(models[2:end])
        A = ifelse.(disk.mantle, 100 .* (neFields[m] ./ max.(ne0, eps()) .- 1), NaN)
        ax = Axis(fig[1, k][1, 1]; aspect = DataAspect(), title = "$(shortName(m)): δln nₑ (%) vs PREM, detectors + rays")
        hidedecorations!(ax); hidespines!(ax)
        hm = heatmap!(ax, disk.xs ./ 1e3, disk.zs ./ 1e3, A; colormap = Reverse(:vik),
            colorrange = robustLimits(A; symmetric = true, hi = 0.99), nan_color = :transparent)
        for r in (1221.5, cmbRadius / 1e3, planetRadius / 1e3)
            lines!(ax, r .* cos.(θ), r .* sin.(θ); color = :gray30, linewidth = 0.6)
        end
        for det in detectors, side in sides
            rays = res[(m, det, side)]["rays"]
            for cθ in (-0.95, -0.8, -0.5, -0.25)
                a, b, _ = rays[argmin([abs(r[3] - cθ) for r in rays])]
                lines!(ax, [a[1], b[1]] ./ 1e3, [a[2], b[2]] ./ 1e3; color = (SOURCE_COLORS[det], 0.8), linewidth = 0.8)
            end
            φd = sourceAngle(det)
            scatter!(ax, [planetRadius / 1e3 * sind(φd)], [planetRadius / 1e3 * cosd(φd)]; marker = :rect,
                markersize = 14, color = SOURCE_COLORS[det], strokecolor = :black, strokewidth = 1)
        end
        Colorbar(fig[1, k][1, 2], hm; label = "%", height = Relative(0.7))
    end
    return fig
end

"δN(E, cosθ) for one detector: rows = sides, columns = sb11 − PREM, sb21 − PREM, sb11 − sb21."
function neutrinoDetectorFigure(res, models, det, sides, channel)
    key = channel == "νe" ? "Ne" : "Nmu"
    fig = Figure(size = (1500, 420 * length(sides) + 60))
    for (row, side) in enumerate(sides)
        r0 = res[(models[1], det, side)]
        E, cθ = r0["energies"], r0["cosθ"]
        N0, N1, N2 = r0[key], res[(models[2], det, side)][key], res[(models[3], det, side)][key]
        for (col, (D, t)) in enumerate(((N1 .- N0, "$(shortName(models[2])) − PREM"),
                (N2 .- N0, "$(shortName(models[3])) − PREM"), (N1 .- N2, "$(shortName(models[2])) − $(shortName(models[3]))")))
            ax = Axis(fig[row, col][1, 1]; xscale = log10, xlabel = "E (GeV)", ylabel = col == 1 ? "cos θ (side $(side > 0 ? "+" : "−"))" : "",
                title = "$det detector, $channel CC: $t")
            lim = max(quantile(abs.(vec(D)), 0.995), eps())
            hm = heatmap!(ax, E, cθ, D; colormap = :balance, colorrange = (-lim, lim))
            Colorbar(fig[row, col][1, 2], hm; label = "events / bin", height = Relative(0.7))
        end
    end
    return fig
end
