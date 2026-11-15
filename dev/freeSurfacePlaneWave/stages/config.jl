# dev/freeSurfacePlaneWave/stages/config.jl
#
# Normal-incidence plane-P-wave-on-a-free-surface benchmark -- NO body
# force, NO point/sheet source anywhere. The plane wave is injected as an
# EXACT initial condition instead: planeWaveUz(z,t) below is the closed-
# form solution of the free (unforced) 1-D wave equation, built from two
# counter-propagating copies of the same pulse (see stages/analytical.jl
# for the derivation), which satisfies the free-surface traction condition
# exactly for ANY pulse shape. FD3/OPT3/OPT5 are seeded with this field at
# t=-dt and t=0 and then propagated with zero forcing from there on, so
# there is no sheet-source/Green's-function calibration (areal force
# density, dx normalisation, cumulative-force integral) anywhere in this
# test -- only the free-surface closure and the propagation scheme matter.
#
# Because this field never depends on x, the whole problem is exactly
# x-independent -- a true normal-incidence plane wave -- and the lateral
# grid edges carry no wave to reflect, by construction, regardless of what
# boundary condition sits there. So this benchmark needs no periodic-BC
# machinery to be free of lateral-boundary doubt: there is simply nothing
# propagating laterally to doubt. The only edge that matters is the bottom
# one (duration-limited, Cerjan-absorbed, same style as every other step);
# the free surface at z=0 is exactly what is under test, against the EXACT
# closed-form solution (no out-of-plane numerical integral, no SPECFEM).
#
# A purely vertical, x-uniform displacement field excites ONLY the P-wave
# (dilatational) motion -- there is no shear response at all without
# x-dependence -- so this is specifically a normal-incidence P-on-free-
# surface test.
#
# Two DEV steps share this file, told apart by DEV_STEP[]:
#   "normal"     -- the actual free surface at z=0 (additive closure, the
#                    thing under test), void above it, Cerjan at the bottom
#                    only. Compared to the two-counter-propagating-waves
#                    (image) closed form.
#   "nosurface"  -- CONTROL: same grid style, but no free surface at all --
#                    fully solid medium, Cerjan on BOTH z ends, no additive
#                    overlay anywhere. Compared to the direct-wave-only
#                    (no image term) closed form. If the growing, x-
#                    localised spurious signal seen in "normal" still
#                    appears here, it has nothing to do with the free-
#                    surface closure at all.
vp0, vs0, rho0 = 6000.0, 3464.0, 2700.0   # same constants as
    # notebooks/HomogeneousElastic2DBenchmark_FreeSurface.ipynb
quickRun = !DEV.full
referenceDx = 500.0
fdOptSpatialRefinement = 2
optSpatialRefinement = 1
dx = referenceDx / fdOptSpatialRefinement   # 250 m
nxHalfWidth = 10   # lateral extent: physically inert (see above), kept small
x = collect((-nxHalfWidth*dx):dx:(nxHalfWidth*dx))
nx = length(x)
z0 = -20.e3             # travel-time phase reference ONLY -- there is no
    # source location any more, but keeping the same value (and the same
    # belowMargin below) keeps the pulse's envelope exactly where it used
    # to be, so the existing bottom-truncation safety margin still applies.
belowMargin = 50.e3     # z-extent below z0: bottom-truncation safety margin
aboveSurfaceVoid = 5.e3 # void cap above the flat surface (FD/OPT stencil padding)
applyFreeSurface = DEV_STEP[] != "nosurface"
z = collect((z0 - belowMargin):dx:(applyFreeSurface ? aboveSurfaceVoid : (z0 + belowMargin)))
nz = length(z)
material = applyFreeSurface ? repeat(reshape(z .<= 0.0, 1, :), nx, 1) : trues(nx, nz)
cerjanCells = round(Int, 24 * referenceDx / dx)
# No Cerjan in x at all (there is no lateral wave to absorb). "normal": Cerjan
# only at the z-lower (bottom) side, none at z-upper (the physical free
# surface). "nosurface": Cerjan at BOTH z ends, symmetric with belowMargin
# above too (checked by earliestBottomArrival below, same margin either side).
cerjan = applyFreeSurface ?
    CerjanBoundarySpec((0, cerjanCells), (0, 0); damping=0.0053 / fdOptSpatialRefinement^2) :
    CerjanBoundarySpec((0, cerjanCells), (0, cerjanCells); damping=0.0053 / fdOptSpatialRefinement^2)
# "nosurface": bypass boundary_geometry's own free-surface auto-detection
# entirely -- with an all-true material it would otherwise treat the bare
# TOP ARRAY ROW as if it were a free surface (findlast of an all-true column
# is just the array's own top), which is exactly the spurious extra surface
# this control run is meant to avoid. free_surface stays `nothing` instead.
bcFD = applyFreeSurface ?
    boundary_geometry(material, (dx, dx); cerjan=cerjan) :
    BoundaryConditionSet(cerjan=cerjan, material_mask=BitArray(material))
duration = quickRun ? 6.0 : 9.0
outputSampling = 0.05
referencePointsPerSWavelength = 16
rickerSharpness = 1.5
minimumAcceptedPointsPerSWavelength = 10
coarsestSpacing = referenceDx
sourceFrequency = rickerSharpness * vs0 / (referencePointsPerSWavelength * coarsestSpacing)
sourceDelay = 1.5 / sourceFrequency
rickerSource(time) = begin
    a = pi * sourceFrequency * (time - sourceDelay)
    (1 - 2a^2) * exp(-a^2)
end
# No body force anywhere in this benchmark. planeWaveUz(z,t) is the exact,
# closed-form, free (unforced) solution: two copies of the SAME pulse
# travelling in opposite directions, same sign (the "same-sign image" of a
# Neumann/traction-free boundary, written as travelling waves instead of a
# static mirror source -- see analytical.jl). It satisfies
# sigma_zz(0,t) = (lambda+2mu) d/dz[...] = 0 EXACTLY, for every t and for
# ANY pulse shape, by construction -- no calibration, no integral. Dropping
# the second term (includeImage=false) gives the exact unbounded (no-
# boundary) direct-wave-only solution used by the "nosurface" control.
planewaveAmplitude = 1.0
planeWaveUz(zd, time; includeImage::Bool=applyFreeSurface) = planewaveAmplitude * (
    rickerSource(time - (zd - z0) / vp0) +
    (includeImage ? rickerSource(time + (zd - z0) / vp0) : 0.0)
)
significantPulseStart = sourceDelay - 1 / sourceFrequency
earliestBottomArrival = significantPulseStart + belowMargin / vp0
@assert duration < earliestBottomArrival "belowMargin is too small for this duration/vp0"
pointsPerSWavelength = (
    FD3=vs0 / (sourceFrequency * dx),
    OPT3=vs0 / (sourceFrequency * (dx / optSpatialRefinement)),
    OPT5=vs0 / (sourceFrequency * (dx / optSpatialRefinement)),
)
@assert minimum(values(pointsPerSWavelength)) >= minimumAcceptedPointsPerSWavelength
# Observation depths: one BELOW z0 (direct downward P wave only -- no free
# surface involved at all, a baseline sanity check of the solver itself)
# and several BETWEEN z0 and the surface (direct + the reflected wave once
# it arrives -- this is what actually tests the free-surface closure).
observerDepths = [z0 - 8e3, z0 + 5e3, z0 + 10e3, z0 + 15e3, -3e3]
@assert all(d -> first(z) < d < last(z), observerDepths)
@show quickRun (nx, nz) dx z0 belowMargin duration sourceFrequency sourceDelay planewaveAmplitude
@show pointsPerSWavelength earliestBottomArrival observerDepths
