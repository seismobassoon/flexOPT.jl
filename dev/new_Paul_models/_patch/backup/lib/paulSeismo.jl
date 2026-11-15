# dev/new_Paul_models/lib/paulSeismo.jl
#
# Whole-Earth 2-D P-SV synthetics through Paul's StagYY models.
#
#   * grid, PREM on the disk (crust + ocean replaced by the LID), closed curved
#     free surface and the FD3 stepper: the same as
#     scripts/earth2D_StagYY_neutrino_seismo.jl (step_closed_surface! is copied
#     here unchanged so that this sub-project does not depend on a script);
#   * heterogeneous models: PREM × (1 + δlnVp, δlnVs, δlnρ), the lateral
#     anomalies of Paul's Pvel, Svel, rho (field / radial mean − 1, so the unit
#     of the files does not matter), mantle only, read with batchStagYY
#     getCartesianField (polar sampling, CMB→surface radial mapping);
#   * source: in-plane moment tensor (moment per unit out-of-plane length,
#     N m / m) as nodal force couples, smooth-step moment history;
#   * receivers every `receiverStepDeg` all around the disk, u_vertical = u·e_r
#     (up) and u_horizontal = u·e_φ (φ increasing, i.e. clockwise from the
#     north pole when x points right and z up).
#
# Station angle φ: 0° = north pole, 90° = +x, 180° = south pole, 270° = −x.

"Regular Cartesian disk centred on the planet."
function seismoDisk(dx, planetRadius, cmbRadius)
    halfCells = ceil(Int, planetRadius / dx) + 3
    xs = collect((-halfCells:halfCells) .* dx)
    X = [x for x in xs, z in xs]
    Z = [z for x in xs, z in xs]
    radius = hypot.(X, Z)
    material = radius .<= planetRadius
    mantle = material .& (radius .>= cmbRadius)
    return (; xs, zs = xs, X, Z, radius, material, mantle)
end

"PREM (planet1D) on the disk; above `mohoRadius` the LID values (no crust, no ocean)."
function premOnDisk(disk, planetRadius, mohoRadius)
    radialTableKm = collect(0.0:1.0:planetRadius / 1e3)
    tableRadii, tableParams = planet1D.compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(
        planet1D.my1DDSMmodel, radialTableKm, "below")
    itp(v) = LinearInterpolation(tableRadii .* 1e3, v; extrapolation_bc = Flat())
    ρI, vpI, vsI = itp(tableParams.ρ), itp(tableParams.Vpv), itp(tableParams.Vsv)
    r = min.(disk.radius, mohoRadius - 1.0)
    m = disk.material
    return (ρ = ifelse.(m, ρI.(r), 0.0), vp = ifelse.(m, vpI.(r), 0.0),
        vs = ifelse.(m, vsI.(r), 0.0))   # g cm⁻³, km s⁻¹
end

"""
    stagyyLateralAnomaly(file, disk, planetRadius, cmbRadius; clip, seamBlendDeg)

field/⟨field⟩_r − 1 of a StagYY scalar file on the disk, zero outside the
mantle, clipped to ±clip. Returns (δ, raw extrema before clipping).
"""
function stagyyLateralAnomaly(file, disk, planetRadius, cmbRadius; clip = 0.10, seamBlendDeg = 0.0)
    s = getCartesianField(file, disk.X, disk.Z;
        interpolation_method = :polar, target_cmb_radius = cmbRadius,
        target_surface_radius = planetRadius, clamp_to_surface = true,
        seam_blend_angle = deg2rad(seamBlendDeg))
    δ = s.field ./ s.avField .- 1
    δ = ifelse.(disk.mantle, δ, 0.0)
    rawExtrema = extrema(δ[disk.mantle])
    return clamp.(δ, -clip, clip), rawExtrema
end

"""
    momentTensorLocal(M0, dip, slip) -> (Mrr, Mtt, Mrt)

In-plane dip-slip on a fault striking out of the plane, dip in degrees,
slip = :thrust or :normal: Mrr = s M0 sin2δ, Mtt = −Mrr, Mrt = −s M0 cos2δ.
"""
function momentTensorLocal(M0, dip, slip)
    s = slip === :thrust ? 1.0 : slip === :normal ? -1.0 :
        throw(ArgumentError("slip must be :thrust or :normal"))
    Mrr = s * M0 * sind(2dip)
    return (Mrr = Mrr, Mtt = -Mrr, Mrt = -s * M0 * cosd(2dip))
end

"Local (e_r, e_t) moment tensor -> Cartesian (Mxx, Mzz, Mxz)."
function momentTensorCartesian(Mloc, er, et)
    M = Mloc.Mrr .* (er * er') .+ Mloc.Mtt .* (et * et') .+ Mloc.Mrt .* (er * et' .+ et * er')
    return (Mxx = M[1, 1], Mzz = M[2, 2], Mxz = M[1, 2])
end

"""
Source at the pole `pole` (:north or :south), `depth` below the surface.
e_r is the outward radial unit vector, e_t = e_φ (direction of increasing φ),
so the two sources are mirror images (rotation by 180°) of each other.
"""
function poleSource(disk, pole, depth, planetRadius, Mloc)
    sgn = pole === :north ? 1.0 : pole === :south ? -1.0 : throw(ArgumentError("pole"))
    zTarget = sgn * (planetRadius - depth)
    I = CartesianIndex(argmin(abs.(disk.xs)), argmin(abs.(disk.zs .- zTarget)))
    disk.mantle[I] || error("source node not in the mantle")
    er = [0.0, sgn]
    et = [sgn, 0.0]          # e_φ at φ = 0 (north) is +x, at φ = 180° (south) −x
    return (index = I, pole = pole, M = momentTensorCartesian(Mloc, er, et), Mloc = Mloc,
        position = (disk.xs[I[1]], disk.zs[I[2]]))
end

"Receivers every `stepDeg` of φ, nearest material node at radius R − depth (or deeper)."
function ringReceivers(disk, planetRadius, stepDeg, depth, dx)
    φs = collect(0.0:stepDeg:360.0-stepDeg)
    idx = map(φs) do φ
        r = planetRadius - depth
        while true
            I = CartesianIndex(argmin(abs.(disk.xs .- r * sind(φ))), argmin(abs.(disk.zs .- r * cosd(φ))))
            disk.material[I] && return I
            r -= dx / 2
        end
    end
    er = [(disk.xs[I[1]], disk.zs[I[2]]) ./ hypot(disk.xs[I[1]], disk.zs[I[2]]) for I in idx]
    eφ = [(e[2], -e[1]) for e in er]
    return (; φ = φs, index = idx, er, eφ)
end

"Epicentral distance (°, 0–180) of station angle φ from the north or south pole."
epicentralDistance(φ, pole) = pole === :north ? min(mod(φ, 360), 360 - mod(φ, 360)) :
    abs(180 - mod(φ, 360))

"Smooth moment step: ∫ Gaussian moment rate (half-width a), sampled on n steps."
function momentHistory(dt, nsteps, f0; t0 = 1.5 / f0)
    a = 1 / (π * f0)                    # moment-rate spectrum exp(−(f/f0)²)
    t = (1:nsteps) .* dt
    rate = @. exp(-((t - t0) / a)^2) / (sqrt(π) * a)
    return cumsum(rate) .* dt           # → 1
end

# ---- FD3 on a closed curved free surface (copied from
#      scripts/earth2D_StagYY_neutrino_seismo.jl, unchanged) -------------------
@inline function seismoDpx(u, m, i, k, h)
    p, q = m[i+1, k], m[i-1, k]
    p && q && return (u[i+1, k] - u[i-1, k]) / (2h)
    p && return (u[i+1, k] - u[i, k]) / h
    q && return (u[i, k] - u[i-1, k]) / h
    return zero(eltype(u))
end
@inline function seismoDpz(u, m, i, k, h)
    p, q = m[i, k+1], m[i, k-1]
    p && q && return (u[i, k+1] - u[i, k-1]) / (2h)
    p && return (u[i, k+1] - u[i, k]) / h
    q && return (u[i, k] - u[i, k-1]) / h
    return zero(eltype(u))
end

function seismoStepClosedSurface!(state::ElasticWaveState2D{T}, edgeFactor) where T
    (; ux, uz, ux_previous, uz_previous, ux_next, uz_next,
       ρ, λ, μ, σxx, σzz, σxz, material) = state
    dx, dz = state.spacing
    nx, nz = size(ux)
    Threads.@threads for k in 2:nz-1
        for i in 2:nx-1
            if !material[i, k]
                σxx[i, k] = σzz[i, k] = σxz[i, k] = zero(T)
                continue
            end
            ux_x = seismoDpx(ux, material, i, k, dx); uz_x = seismoDpx(uz, material, i, k, dx)
            ux_z = seismoDpz(ux, material, i, k, dz); uz_z = seismoDpz(uz, material, i, k, dz)
            σxx[i, k] = (λ[i, k] + 2μ[i, k]) * ux_x + λ[i, k] * uz_z
            σzz[i, k] = λ[i, k] * ux_x + (λ[i, k] + 2μ[i, k]) * uz_z
            σxz[i, k] = μ[i, k] * (ux_z + uz_x)
        end
    end
    dt2 = state.dt^2
    Threads.@threads for k in 2:nz-1
        for i in 2:nx-1
            if !material[i, k]
                ux_next[i, k] = uz_next[i, k] = zero(T)
                continue
            end
            σxx_r = material[i+1, k] ? (σxx[i, k] + σxx[i+1, k]) / 2 : zero(T)
            σxx_l = material[i-1, k] ? (σxx[i, k] + σxx[i-1, k]) / 2 : zero(T)
            σxz_t = material[i, k+1] ? (σxz[i, k] + σxz[i, k+1]) / 2 : zero(T)
            σxz_b = material[i, k-1] ? (σxz[i, k] + σxz[i, k-1]) / 2 : zero(T)
            σxz_r = material[i+1, k] ? (σxz[i, k] + σxz[i+1, k]) / 2 : zero(T)
            σxz_l = material[i-1, k] ? (σxz[i, k] + σxz[i-1, k]) / 2 : zero(T)
            σzz_t = material[i, k+1] ? (σzz[i, k] + σzz[i, k+1]) / 2 : zero(T)
            σzz_b = material[i, k-1] ? (σzz[i, k] + σzz[i, k-1]) / 2 : zero(T)
            ax = ((σxx_r - σxx_l) / dx + (σxz_t - σxz_b) / dz) / ρ[i, k]
            az = ((σxz_r - σxz_l) / dx + (σzz_t - σzz_b) / dz) / ρ[i, k]
            f = edgeFactor[i, k]
            ux_next[i, k] = f * (2ux[i, k] - ux_previous[i, k] + dt2 * ax)
            uz_next[i, k] = f * (2uz[i, k] - uz_previous[i, k] + dt2 * az)
        end
    end
    state.ux_previous, state.ux, state.ux_next = state.ux, state.ux_next, state.ux_previous
    state.uz_previous, state.uz, state.uz_next = state.uz, state.uz_next, state.uz_previous
    state.step += 1
    state.time += state.dt
    state
end

"""
Add the moment-tensor force couples for one step: M(t) = M · m(t) per unit
length; nodal force F = M_ij/(2h) on the four neighbours (body force
f = −∇·(M δ)), acceleration F/(ρ dx dz), added as dt²·a like add_ricker_source!.
"""
function addMomentTensor!(state::ElasticWaveState2D{T}, src, mt) where T
    i, k = Tuple(src.index)
    dx, dz = state.spacing
    dt2 = state.dt^2
    (; Mxx, Mzz, Mxz) = src.M
    function couple!(ii, kk, Fx, Fz)
        state.material[ii, kk] || error("source force-couple node outside the material")
        c = dt2 * mt / (state.ρ[ii, kk] * dx * dz)
        state.ux[ii, kk] += T(c * Fx)
        state.uz[ii, kk] += T(c * Fz)
    end
    couple!(i + 1, k, Mxx / (2dx), Mxz / (2dx))
    couple!(i - 1, k, -Mxx / (2dx), -Mxz / (2dx))
    couple!(i, k + 1, Mxz / (2dz), Mzz / (2dz))
    couple!(i, k - 1, -Mxz / (2dz), -Mzz / (2dz))
    return state
end

"""
    runSeismoFD3(model, disk, src, rec, cfg) -> Dict

model = (ρ g cm⁻³, vp, vs km s⁻¹) on the disk. cfg: dt, duration, period,
outputSampling, surfaceDamping, snapshotTimes. Returns time, uVert, uHoriz
(nt × nstations), snapshots of |u|, and run info.
"""
function runSeismoFD3(model, disk, src, rec, cfg; label = "")
    state = prepare_elastic_wave_2d((ρ = model.ρ, Vpv = model.vp, Vsv = model.vs),
        (cfg.dx, cfg.dx); material_mask = disk.material,
        boundary_conditions = (cerjan = nothing,), dt = cfg.dt)
    m = state.material
    edgeFactor = ones(Float32, size(m))
    for k in 2:size(m, 2)-1, i in 2:size(m, 1)-1
        m[i, k] && !(m[i+1, k] && m[i-1, k] && m[i, k+1] && m[i, k-1]) &&
            (edgeFactor[i, k] = 1 - Float32(cfg.surfaceDamping))
    end
    steps = ceil(Int, cfg.duration / state.dt)
    history = momentHistory(state.dt, steps, 1 / cfg.period)
    stride = max(1, round(Int, cfg.outputSampling / state.dt))
    nstore = steps ÷ stride + 1
    nrec = length(rec.index)
    ux = zeros(Float32, nstore, nrec); uz = similar(ux)
    times = zeros(nstore)
    snapshots = Dict{Float64,Matrix{Float32}}()
    stored = 1
    elapsed = @elapsed for step in 1:steps
        seismoStepClosedSurface!(state, edgeFactor)
        addMomentTensor!(state, src, history[step])
        if step % stride == 0 && stored < nstore
            stored += 1
            times[stored] = state.time
            for (j, I) in enumerate(rec.index)
                ux[stored, j] = state.ux[I]; uz[stored, j] = state.uz[I]
            end
        end
        for ts in cfg.snapshotTimes
            abs(state.time - ts) < state.dt / 2 && (snapshots[ts] = Float32.(hypot.(state.ux, state.uz)))
        end
        if step % 2000 == 0
            amp = maximum(abs, state.uz)
            isfinite(amp) || error("$label FD3 diverged at t = $(state.time) s")
            @info label t = round(state.time) max_uz = amp progress = "$(round(Int, 100step / steps)) %"
        end
    end
    uVert = similar(ux); uHoriz = similar(ux)
    for j in 1:nrec
        (ex, ez), (fx, fz) = rec.er[j], rec.eφ[j]
        uVert[:, j] .= ux[:, j] .* ex .+ uz[:, j] .* ez
        uHoriz[:, j] .= ux[:, j] .* fx .+ uz[:, j] .* fz
    end
    @info "$label done" dt = state.dt steps elapsed
    return Dict{String,Any}("time" => times[1:stored], "uVert" => uVert[1:stored, :],
        "uHoriz" => uHoriz[1:stored, :], "snapshots" => snapshots,
        "dt" => Float64(state.dt), "steps" => steps, "elapsed" => elapsed)
end
