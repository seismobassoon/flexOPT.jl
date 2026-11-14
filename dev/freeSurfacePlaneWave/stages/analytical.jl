# dev/freeSurfacePlaneWave/stages/analytical.jl
# Exact closed-form displacement for a homogeneous half-space (z<=0 solid,
# free surface at z=0) driven by a force PER UNIT AREA s(t) at depth z0,
# uniform in x -- i.e. a genuine normal-incidence plane P wave. This is the
# method-of-images solution of the 1-D wave equation rho*u_tt = kappa*u_zz +
# s(t)*delta(z-z0): the free (traction-free, Neumann) surface at z=0 is
# represented by an EVEN image source at z_image=-z0 (same sign -- a
# Dirichlet/fixed boundary would need the opposite sign). No SPECFEM, no
# out-of-plane numerical integral (c.f. src/elasticGreens2D.jl's full 2-D/
# 3-D point-source Green's function): this is the exact solution, limited
# only by how finely `dt_fine` samples the cumulative force integral below.
function _cumulative_integral(f, t0, dt, t1)
    ts = collect(t0:dt:t1)
    s = f.(ts)
    S = zeros(Float64, length(ts))
    for i in 2:length(ts)
        S[i] = S[i - 1] + dt * (s[i - 1] + s[i]) / 2
    end
    return ts, S
end

function _sample_linear_zero(ts, S, t)
    (t < ts[1] || t > ts[end]) && return 0.0
    position = (t - ts[1]) / (ts[2] - ts[1]) + 1
    lower = clamp(floor(Int, position), 1, length(ts))
    lower == length(ts) && return S[end]
    fraction = position - lower
    return S[lower] + fraction * (S[lower + 1] - S[lower])
end

"""
    freeSurfaceNormalIncidenceP(times, depths; z0, vp, rho, dx, line_force,
        force_time_function, source_frequency, tail_periods=6.0, dt_fine=nothing)

Exact displacement u_z(depth, time) for a vertical force `line_force *
force_time_function(t)` (N/m, the SAME per-column magnitude convention used
by every other benchmark here) applied uniformly at EVERY column at depth
`z0`: pass `line_force = sourceForce` and `dx`, and this function divides by
`dx` itself to get the areal force density (N/m^2) the 1-D image-source
formula needs.
"""
function freeSurfaceNormalIncidenceP(times, depths; z0, vp, rho, dx,
    line_force, force_time_function, source_frequency,
    tail_periods::Real=6.0, dt_fine=nothing)
    t = Float64.(collect(times))
    length(t) >= 2 || throw(ArgumentError("at least two times are required"))
    dtf = isnothing(dt_fine) ? minimum(diff(t)) / 8 : Float64(dt_fine)
    margin = tail_periods / Float64(source_frequency)
    t0 = first(t) - margin
    t1 = last(t) + margin
    arealForceDensity = Float64(line_force) / Float64(dx)
    ts, S = _cumulative_integral(
        tau -> arealForceDensity * force_time_function(tau), t0, dtf, t1)
    zImage = -Float64(z0)
    prefactor = inv(2 * Float64(rho) * Float64(vp))
    depthList = collect(Float64, depths)
    u = zeros(Float64, length(t), length(depthList))
    for (iz, zd) in enumerate(depthList)
        delayDirect = abs(zd - z0) / vp
        delayImage = abs(zd - zImage) / vp
        for (it, time) in enumerate(t)
            u[it, iz] = prefactor * (
                _sample_linear_zero(ts, S, time - delayDirect) +
                _sample_linear_zero(ts, S, time - delayImage))
        end
    end
    return (time=t, depths=depthList, displacement=u,
        convention=:image_source_1d_free_surface_normal_incidence_P)
end
