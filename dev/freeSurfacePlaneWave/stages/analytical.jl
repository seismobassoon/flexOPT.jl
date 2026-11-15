# dev/freeSurfacePlaneWave/stages/analytical.jl
# Exact closed-form displacement for a homogeneous half-space (z<=0 solid,
# free surface at z=0): a normal-incidence plane P wave expressed directly
# as two counter-propagating travelling waves, NOT as a point/sheet force
# convolved against a Green's function. No "areal force density" or
# dx-normalisation convention to get right, no cumulative-force integral
# S(t), no interpolation -- the formula is evaluated in closed form at any
# (z,t) directly (config.jl's planeWaveUz), and is EXACTLY what is used to
# build the simulations' own initial condition (see opt_propagate.jl /
# fd_propagate.jl), so this is as much a definition as a cross-check: what
# is being verified is that FD3/OPT3/OPT5, started from this exact field at
# t=-dt and t=0 with NO forcing at all afterwards, keep reproducing it.
#
# Derivation (method of images, expressed as plane waves instead of a
# static mirror source). Write a rightward-travelling pulse
# w(t - (z-z0)/vp) and a leftward-travelling pulse of the SAME sign and
# SAME shape w(t + (z-z0)/vp). Since u_z = (either term), and
#   sigma_zz = (lambda+2mu) d(u_z)/dz,
# term by term: d/dz[w(t-(z-z0)/vp)] = -(1/vp) w'(t-(z-z0)/vp), and
#               d/dz[w(t+(z-z0)/vp)] = +(1/vp) w'(t+(z-z0)/vp).
# At z=0 these two derivatives are exactly equal in magnitude and opposite
# in sign (both pulses share the same argument magnitude |z-z0|/vp there),
# so their sum is EXACTLY ZERO: sigma_zz(0,t) = 0 for every t, for ANY
# pulse shape w -- this is the free-surface traction condition, satisfied
# by construction, with no reference to w's particular shape at all. This
# is the "same-sign (even) image" of a Neumann boundary, the plane-wave
# analogue of the static method-of-images mirror source. Dropping the
# second term instead gives the exact unbounded (no-boundary) direct-wave-
# only solution used by the "nosurface" control (trivially satisfies the
# free 1-D wave equation everywhere, no boundary involved at all).
#
# A purely vertical, x-independent field like this excites ONLY P-wave
# (dilatational) motion: u_x is identically zero everywhere (no shear is
# excited at all without x-dependence), confirmed directly from the 2-D
# elastic equations of motion with partial_x == 0 throughout.
"""
    planeWaveNormalIncidenceP(times, depths; z0, vp, amplitude, pulse, include_image=true)

Exact normal-incidence plane-P-wave displacement u_z(depth, time), z<=0.
`z0` is a travel-time phase reference only (there is no source location in
this benchmark). `pulse` is the scalar waveform (config.jl's `rickerSource`,
shared with the simulations' own initial condition).
"""
function planeWaveNormalIncidenceP(times, depths; z0, vp, amplitude, pulse,
    include_image::Bool=true)
    t = Float64.(collect(times))
    depthList = collect(Float64, depths)
    u = zeros(Float64, length(t), length(depthList))
    for (iz, zd) in enumerate(depthList)
        for (it, time) in enumerate(t)
            direct = pulse(time - (zd - z0) / vp)
            image = include_image ? pulse(time + (zd - z0) / vp) : 0.0
            u[it, iz] = amplitude * (direct + image)
        end
    end
    return (time=t, depths=depthList, displacement=u,
        convention=include_image ?
            :counterpropagating_plane_waves_free_surface_normal_incidence_P :
            :direct_plane_wave_only_no_surface)
end
