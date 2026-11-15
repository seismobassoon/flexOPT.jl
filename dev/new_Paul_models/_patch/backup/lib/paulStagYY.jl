# dev/new_Paul_models/lib/paulStagYY.jl
#
# Thin layer over src/batchFiles/batchStagYY.jl (module commonBatchs). Every
# byte is read by batchStagYY itself (read_magic, read_header, readField,
# readStagYYFilesAverage) and every polar -> Cartesian resampling uses its
# sampleStagYYPolar. This file only adds what batchStagYY does not do yet:
#
#   * vector files (`vp`: v_theta, v_phi, v_r, p). StagYY writes one extra word
#     (a scale factor) right after the header and stores each block with one
#     ghost column/row, so they are read here with
#     read_header -> read(scalefac) -> readField(...; xyghost = 1).
#     (readStagYYFilesAverage would mis-read them: no scalefac, xyghost = 0.)
#   * surface files (`cs`, `hf`, `cr`: nz <= 2) are functions of angle only.
#   * names for fields not in batchStagYY's table (Tpotl, cs, hf, cr, vp).
#
# Nothing in src/ is modified.

SECONDS_PER_GYR = 3.15576e16
MPS_TO_CMPYR = 3.15576e9

PAUL_EXTRA_NAMES = Dict(
    "Tpotl" => "Potential temperature",
    "cs" => "Topography",
    "hf" => "Heat flux",
    "cr" => "Crustal thickness",
    "vp" => "Velocity and pressure",
)

"Long name of a StagYY file code (batchStagYY table first, then PAUL_EXTRA_NAMES)."
function stagyyLongName(code::AbstractString)
    haskey(PAUL_EXTRA_NAMES, code) && return PAUL_EXTRA_NAMES[code]
    for (k, v) in commonBatchs.fileExt2FieldNameStagYY
        strip(k) == code && return replace(strip(v), '_' => ' ')
    end
    return code
end

# Display scale and unit (value shown = scale * value in file). Anything not
# listed is shown as stored, without unit.
PAUL_DISPLAY = Dict(
    "t" => (scale = 1.0, unit = "K"),
    "Tpotl" => (scale = 1.0, unit = "K"),
    "rho" => (scale = 1.0, unit = "kg m⁻³"),
    "eta" => (scale = 1.0, unit = "Pa s"),
    "age" => (scale = 1 / SECONDS_PER_GYR, unit = "Gyr"),
    "Pdyn" => (scale = 1e-6, unit = "MPa"),
    "str" => (scale = 1e-6, unit = "MPa"),
    "ed" => (scale = 1.0, unit = "s⁻¹"),
    "k" => (scale = 1.0, unit = "W m⁻¹ K⁻¹"),
    "cs" => (scale = 1e-3, unit = "km"),
    "cr" => (scale = 1e-3, unit = "km"),
    "hf" => (scale = 1e3, unit = "mW m⁻²"),
)
displayScale(code) = haskey(PAUL_DISPLAY, code) ? PAUL_DISPLAY[code].scale : 1.0
displayUnit(code) = haskey(PAUL_DISPLAY, code) ? PAUL_DISPLAY[code].unit : ""

"""
    listStagYYSnapshots(outputDir) -> Dict(snapshot => Dict(code => path))

Binary StagYY files are `<prefix>_<code><5-digit snapshot>` (e.g.
`DeepSpace_phs1103140` -> code `phs11`, snapshot `03140`). `.vts`/`.pvd`
files are ignored.
"""
function listStagYYSnapshots(outputDir::AbstractString)
    snaps = Dict{String,Dict{String,String}}()
    for name in readdir(outputDir)
        m = match(r"^(.+)_([A-Za-z][A-Za-z0-9]*?)(\d{5})$", name)
        m === nothing && continue
        code, snap = String(m.captures[2]), String(m.captures[3])
        get!(snaps, snap, Dict{String,String}())[code] = joinpath(outputDir, name)
    end
    return snaps
end

"""
    readStagYYRaw(file)

Header + data of any StagYY binary file through batchStagYY. `data` has size
`(nval, nx, ny, nz, nb)` (ghost points removed by readField), already
multiplied by the vector scale factor.
"""
function readStagYYRaw(file::AbstractString)
    magic, ilen, byteReversed, io = commonBatchs.read_magic(file)
    if byteReversed
        close(io)
        error("byte-reversed StagYY file not handled by batchStagYY.readField: $file")
    end
    try
        magic, nval = commonBatchs.evaluate_nval_from_magicNumber(magic)
        intT, floatT, dx, dy, nx, ny, nz, nb, nnx, nny, nnz, nnb, rcmb, iStep, time, zc, spherical =
            commonBatchs.read_header(io, ilen, magic)
        scalefac = nval > 1 ? read(io, floatT) : one(floatT)
        raw = commonBatchs.readField(io, floatT, nval, nx, ny, nz, nb, nnx, nny, nnz, nnb;
            xyghost = nval > 1 ? 1 : 0)
        eof(io) || @warn "unread bytes left" file remaining = filesize(file) - position(io)
        return (; nval, nx, ny, nz, dy, rcmb, iStep, time, zc, spherical,
            data = Float64.(raw) .* scalefac)
    finally
        close(io)
    end
end

"Header only (cheap): time step, time, grid sizes."
function stagyyHeader(file::AbstractString)
    magic, ilen, _, io = commonBatchs.read_magic(file)
    try
        magic, nval = commonBatchs.evaluate_nval_from_magicNumber(magic)
        intT, floatT, dx, dy, nx, ny, nz, nb, nnx, nny, nnz, nnb, rcmb, iStep, time, zc, spherical =
            commonBatchs.read_header(io, ilen, magic)
        return (; nval, nx, ny, nz, dy, rcmb, iStep, time, zc, spherical)
    finally
        close(io)
    end
end

"Polar nodes exactly as readStagYYFilesAverage builds them (θshift = 0)."
function stagyyPolarNodes(h)
    nz, ny = h.nz, h.ny
    r = h.rcmb .+ h.zc[1, 1:nz+1]
    rmid = h.rcmb .+ h.zc[2, 1:nz]
    phi = ((1:ny+1) .- (ny + 1)) .* h.dy
    phimid = ((1:ny) .- (ny + 1)) .* h.dy .+ 0.5h.dy
    return (angle = vcat(phi[1], phimid, phi[end]), radius = vcat(r[1], rmid, r[end]),
        phimid = collect(phimid), rmid = collect(rmid))
end

"Pad an (ny, nz) cell field to the (ny+2, nz+2) node matrix of batchStagYY."
function padPolar(field::AbstractMatrix)
    ny, nz = size(field)
    out = zeros(Float64, ny + 2, nz + 2)
    out[2:ny+1, 2:nz+1] .= field
    periodic = 0.5 .* (field[1, :] .+ field[end, :])
    out[1, 2:nz+1] .= periodic
    out[end, 2:nz+1] .= periodic
    out[2:ny+1, 1] .= field[:, 1]
    out[2:ny+1, end] .= field[:, end]
    out[1, 1] = out[end, 1] = periodic[1]
    out[1, end] = out[end, end] = periodic[end]
    return out
end

"Regular Cartesian grid covering the disk of radius `rmax`."
function diskGrid(rmax; n = 701)
    xs = collect(range(-rmax, rmax; length = n))
    X = [x for x in xs, y in xs]
    Y = [y for x in xs, y in xs]
    return (; xs, X, Y, R = hypot.(X, Y))
end

"batchStagYY polar resampling, NaN outside the annulus."
function sampleDisk(M, nodes, g)
    v = commonBatchs.sampleStagYYPolar(M, nodes.angle, nodes.radius, g.X, g.Y)
    v[(g.R .< nodes.radius[1]) .| (g.R .> nodes.radius[end])] .= NaN
    return v
end

"Radial statistics over the interior rows (real cells) of a padded matrix."
function radialStats(M)
    inner = @view M[2:end-1, 2:end-1]
    return (mean = vec(mean(inner; dims = 1)), min = vec(minimum(inner; dims = 1)),
        max = vec(maximum(inner; dims = 1)))
end

"Colour limits from quantiles of the finite values (robust to a few spikes)."
function robustLimits(v; lo = 0.005, hi = 0.995, symmetric = false)
    w = filter(isfinite, vec(v))
    isempty(w) && return (-1.0, 1.0)
    if symmetric
        a = quantile(abs.(w), hi)
        a = a > 0 ? a : 1.0
        return (-a, a)
    end
    a, b = quantile(w, [lo, hi])
    if !(b > a)
        a, b = extrema(w)
        b > a || return (a - 0.5, a + 0.5)
    end
    return (a, b)
end

"log10 display when the field is strictly positive and spans > 3 decades."
useLog(values) = (m = minimum(values); M = maximum(values); m > 0 && M / m > 1e3)

fieldColormap(code) = code in ("t", "Tpotl") ? :inferno : :viridis

function diskPanel!(gp, g, V; title, colormap, colorrange, label, nodes)
    ax = Axis(gp[1, 1]; aspect = DataAspect(), title = title, titlesize = 13)
    hidedecorations!(ax)
    hidespines!(ax)
    hm = heatmap!(ax, g.xs ./ 1e3, g.xs ./ 1e3, V; colormap = colormap,
        colorrange = colorrange, nan_color = :transparent)
    θ = range(0, 2π; length = 361)
    for r in (nodes.radius[1], nodes.radius[end])
        lines!(ax, r / 1e3 .* cos.(θ), r / 1e3 .* sin.(θ); color = :black, linewidth = 0.6)
    end
    Colorbar(gp[1, 2], hm; label = label, height = Relative(0.75))
    return ax
end

# ─────────────────────────────────────────────────────────────────────────────
# Figures
# ─────────────────────────────────────────────────────────────────────────────

timeLabel(h) = @sprintf("t = %.4f Gyr (step %d)", h.time / SECONDS_PER_GYR, h.iStep)

"""
    processScalar(file, code, g; seamBlendDeg)

2-D scalar field through readStagYYFilesAverage (batchStagYY), transformed for
display (scale, optional log10), its lateral anomaly (minus the radial mean of
the displayed quantity) and radial statistics, all resampled on grid `g`.
"""
function processScalar(file, code, g; seamBlendDeg = 0.0)
    h = stagyyHeader(file)
    s = readStagYYFilesAverage(file; seam_blend_angle = deg2rad(seamBlendDeg))
    nodes = (angle = s.angleNodes, radius = s.radiusNodes,
        rmid = stagyyPolarNodes(h).rmid)
    M = displayScale(code) .* Float64.(s.fieldMatrix)
    logScale = useLog(@view M[2:end-1, 2:end-1])
    logScale && (M = log10.(M))
    stats = radialStats(M)
    padMean = vcat(stats.mean[1], stats.mean, stats.mean[end])   # radial endpoints
    A = M .- padMean'
    unit = displayUnit(code)
    label = (logScale ? "log₁₀ " : "") * code * (isempty(unit) ? "" : " ($unit)")
    return (; h, nodes, logScale, label, stats,
        V = sampleDisk(M, nodes, g), Anom = sampleDisk(A, nodes, g))
end

function scalarFigure(p, code, model, snap, g)
    name = stagyyLongName(code)
    fig = Figure(size = (1700, 600))
    Label(fig[0, 1:3], "$model · snapshot $snap · $(timeLabel(p.h)) · $name [$code]";
        fontsize = 17, font = :bold)
    diskPanel!(fig[1, 1], g, p.V; title = name, colormap = fieldColormap(code),
        colorrange = robustLimits(p.V), label = p.label, nodes = p.nodes)
    diskPanel!(fig[1, 2], g, p.Anom; title = "lateral anomaly (minus radial mean)",
        colormap = :vik, colorrange = robustLimits(p.Anom; symmetric = true),
        label = "δ " * p.label, nodes = p.nodes)
    depth = (p.nodes.radius[end] .- p.nodes.rmid) ./ 1e3
    ax = Axis(fig[1, 3]; xlabel = p.label, ylabel = "depth (km)", yreversed = true,
        title = "radial mean, min, max")
    lines!(ax, p.stats.min, depth; color = :gray60, linestyle = :dash, label = "min")
    lines!(ax, p.stats.max, depth; color = :gray30, linestyle = :dash, label = "max")
    lines!(ax, p.stats.mean, depth; color = :black, linewidth = 2, label = "mean")
    axislegend(ax; position = :rb, framevisible = false)
    colsize!(fig.layout, 3, Relative(0.2))
    return fig
end

"vp file: v_φ (tangential), v_r (radial), |v| in cm/yr and p in GPa."
function vectorFigure(file, model, snap, g)
    raw = readStagYYRaw(file)
    raw.nval == 4 || error("expected a 4-component vp file, got nval = $(raw.nval)")
    nodes = stagyyPolarNodes(raw)
    vphi = raw.data[2, 1, :, :, 1] .* MPS_TO_CMPYR
    vr = raw.data[3, 1, :, :, 1] .* MPS_TO_CMPYR
    speed = hypot.(vphi, vr)
    pGPa = raw.data[4, 1, :, :, 1] .* 1e-9
    panels = (("v_φ (tangential)", vphi, "v_φ (cm/yr)", :vik, true),
        ("v_r (radial)", vr, "v_r (cm/yr)", :vik, true),
        ("|v|", speed, "|v| (cm/yr)", :viridis, false),
        ("pressure", pGPa, "p (GPa)", :viridis, false))
    fig = Figure(size = (2000, 560))
    Label(fig[0, 1:4], "$model · snapshot $snap · $(timeLabel(raw)) · $(stagyyLongName("vp")) [vp]";
        fontsize = 17, font = :bold)
    for (k, (title, F, label, cmap, sym)) in enumerate(panels)
        V = sampleDisk(padPolar(F), nodes, g)
        diskPanel!(fig[1, k], g, V; title = title, colormap = cmap,
            colorrange = robustLimits(V; symmetric = sym), label = label, nodes = nodes)
    end
    return fig
end

"Surface files (cs, hf, cr): one line per stored z-layer, versus angle."
function surfaceFigure(file, code, model, snap)
    raw = readStagYYRaw(file)
    nodes = stagyyPolarNodes(raw)
    angleDeg = rad2deg.(nodes.phimid .+ 2π)
    unit = displayUnit(code)
    fig = Figure(size = (1400, 300 + 220 * raw.nz))
    Label(fig[0, 1], "$model · snapshot $snap · $(timeLabel(raw)) · $(stagyyLongName(code)) [$code]";
        fontsize = 17, font = :bold)
    for iz in 1:raw.nz
        ax = Axis(fig[iz, 1]; xlabel = iz == raw.nz ? "angle (°)" : "",
            ylabel = code * (isempty(unit) ? "" : " ($unit)"),
            title = "z-index $iz of $(raw.nz)", xticks = 0:30:360)
        lines!(ax, angleDeg, displayScale(code) .* raw.data[1, 1, :, iz, 1]; color = :black)
        xlims!(ax, 0, 360)
    end
    return fig
end

"Contact sheet: every 2-D scalar field of one snapshot."
function overviewFigure(processed, model, snap, g; ncol = 7)
    codes = sort(collect(keys(processed)))
    nrow = cld(length(codes), ncol)
    fig = Figure(size = (ncol * 270, nrow * 330 + 60))
    h = processed[codes[1]].h
    Label(fig[0, 1:ncol], "$model · snapshot $snap · $(timeLabel(h)) · all 2-D fields";
        fontsize = 20, font = :bold)
    for (k, code) in enumerate(codes)
        p = processed[code]
        gl = fig[cld(k, ncol), mod1(k, ncol)] = GridLayout()
        ax = Axis(gl[1, 1]; aspect = DataAspect(), title = "$code: $(stagyyLongName(code))",
            titlesize = 11)
        hidedecorations!(ax)
        hidespines!(ax)
        hm = heatmap!(ax, g.xs ./ 1e3, g.xs ./ 1e3, p.V; colormap = fieldColormap(code),
            colorrange = robustLimits(p.V), nan_color = :transparent)
        Colorbar(gl[2, 1], hm; vertical = false, height = 7, ticklabelsize = 9,
            label = p.label, labelsize = 10)
        rowgap!(gl, 2)
    end
    return fig
end

# ─────────────────────────────────────────────────────────────────────────────
# Pseudo electron number density  nₑ = N_A ρ Yₑ   (Yₑ = Z/A, electrons per nucleon)
# ─────────────────────────────────────────────────────────────────────────────

"""
    electronYield(bs; Ye_hz, Ye_morb, bsIsVolume, ΔρMorb)

Yₑ of a basalt–harzburgite mechanical mixture. Yₑ mixes by MASS: if `bs` is a
volume fraction it is converted with a constant density excess `ΔρMorb`
(ρ_MORB/ρ_hz − 1). Points outside [0,1] are clamped.
"""
function electronYield(bs; Ye_hz, Ye_morb, bsIsVolume = true, ΔρMorb = 0.02)
    f = clamp.(bs, 0.0, 1.0)
    w = bsIsVolume ? f .* (1 + ΔρMorb) ./ (1 .+ f .* ΔρMorb) : f
    return Ye_hz .+ w .* (Ye_morb - Ye_hz)
end

lateralAnomaly(M) = (s = radialStats(M); M .- vcat(s.mean[1], s.mean, s.mean[end])')

"rms over angle (interior rows) of each radial column."
rmsByRadius(A) = vec(sqrt.(mean(abs2, @view(A[2:end-1, 2:end-1]); dims = 1)))

"""
    electronDensityFigure(files, model, snap, g; Ye_hz, Ye_morb, bsIsVolume, ΔρMorb)

Needs the `rho` and `bs` files of one snapshot. nₑ is shown as ρ·Yₑ in
mol e⁻ cm⁻³ (× N_A = 6.022e23 for electrons cm⁻³), computed with
commonBatchs.electronDensity. Panels: Yₑ, nₑ, δln nₑ, δln Yₑ (lateral, %),
and rms(δln ρ), rms(δln Yₑ), rms(δln nₑ) versus depth.
"""
function electronDensityFigure(files, model, snap, g; Ye_hz, Ye_morb,
        bsIsVolume = true, ΔρMorb = 0.02, seamBlendDeg = 0.0)
    h = stagyyHeader(files["rho"])
    readS(code) = readStagYYFilesAverage(files[code]; seam_blend_angle = deg2rad(seamBlendDeg))
    sρ, sbs = readS("rho"), readS("bs")
    nodes = (angle = sρ.angleNodes, radius = sρ.radiusNodes, rmid = stagyyPolarNodes(h).rmid)
    ρ = Float64.(sρ.fieldMatrix) ./ 1e3                         # g cm⁻³
    Ye = electronYield(Float64.(sbs.fieldMatrix); Ye_hz, Ye_morb, bsIsVolume, ΔρMorb)
    ne = electronDensity(ρ, Ye)                                 # mol e⁻ cm⁻³ (src/batchStagYY)
    dlnρ = 100 .* lateralAnomaly(log.(ρ))
    dlnYe = 100 .* lateralAnomaly(log.(Ye))
    dlnne = 100 .* lateralAnomaly(log.(ne))

    fig = Figure(size = (2250, 600))
    Label(fig[0, 1:5], "$model · snapshot $snap · $(timeLabel(h)) · pseudo electron density nₑ = N_A ρ Yₑ" *
        @sprintf("  (Yₑ hz = %.5f, MORB = %.5f, bs as %s fraction)", Ye_hz, Ye_morb,
            bsIsVolume ? "volume" : "mass"); fontsize = 17, font = :bold)
    YeV = sampleDisk(Ye, nodes, g)
    diskPanel!(fig[1, 1], g, YeV; title = "Yₑ = Z/A", colormap = :viridis,
        colorrange = robustLimits(YeV), label = "Yₑ", nodes = nodes)
    neV = sampleDisk(ne, nodes, g)
    diskPanel!(fig[1, 2], g, neV; title = "nₑ", colormap = :viridis,
        colorrange = robustLimits(neV), label = "nₑ (mol e⁻ cm⁻³)", nodes = nodes)
    for (k, (A, title, label)) in enumerate(((dlnne, "δln nₑ (lateral)", "δln nₑ (%)"),
            (dlnYe, "δln Yₑ (lateral, composition only)", "δln Yₑ (%)")))
        V = sampleDisk(A, nodes, g)
        diskPanel!(fig[1, 2 + k], g, V; title = title, colormap = :vik,
            colorrange = robustLimits(V; symmetric = true), label = label, nodes = nodes)
    end
    depth = (nodes.radius[end] .- nodes.rmid) ./ 1e3
    ax = Axis(fig[1, 5]; xlabel = "rms lateral anomaly (%)", ylabel = "depth (km)",
        yreversed = true, title = "who makes δnₑ?", xscale = log10)
    for (A, lab, col) in ((dlnne, "δln nₑ", :black), (dlnρ, "δln ρ", :steelblue),
            (dlnYe, "δln Yₑ", :firebrick))
        lines!(ax, max.(rmsByRadius(A), 1e-6), depth; color = col, linewidth = 2, label = lab)
    end
    axislegend(ax; position = :rb, framevisible = false)
    colsize!(fig.layout, 5, Relative(0.16))
    return fig
end

# (Z, molar mass g/mol) of the NCFMAS oxides
OXIDE_Z_A = Dict("SiO2" => (30, 60.084), "MgO" => (20, 40.304), "FeO" => (34, 71.844),
    "CaO" => (28, 56.077), "Al2O3" => (50, 101.961), "Na2O" => (30, 61.979))

"Yₑ = Σ xᵢ Zᵢ / Σ xᵢ Aᵢ from a composition given in mol% (Dict oxide => mol%)."
function yeFromMolPercent(comp::AbstractDict)
    Z = sum(x * OXIDE_Z_A[o][1] for (o, x) in comp)
    A = sum(x * OXIDE_Z_A[o][2] for (o, x) in comp)
    return Z / A
end
