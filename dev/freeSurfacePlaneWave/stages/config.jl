# dev/freeSurfacePlaneWave/stages/config.jl
#
# Normal-incidence plane-P-wave-on-a-free-surface benchmark: a vertical body
# force applied IDENTICALLY at every x column (not Gaussian-weighted in x
# like every other source in this repository) at depth z0, in a homogeneous
# medium. Because neither the model nor the source ever depends on x, the
# whole problem is exactly x-independent -- a true normal-incidence plane
# wave -- and the lateral grid edges carry no wave to reflect, by
# construction, regardless of what boundary condition sits there. So this
# benchmark needs no periodic-BC machinery to be free of lateral-boundary
# doubt: there is simply nothing propagating laterally to doubt. The only
# edge that matters is the bottom one (duration-limited, Cerjan-absorbed,
# same style as every other step); the free surface at z=0 is exactly what
# is under test, against the EXACT closed-form image-source solution in
# stages/analytical.jl (no out-of-plane numerical integral, no SPECFEM).
#
# A vertical, x-uniform force excites ONLY the P-wave (dilatational) motion
# -- there is no shear response at all without x-dependence -- so this is
# specifically a normal-incidence P-on-free-surface test.
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
z0 = -20.e3            # source depth, same convention as the existing notebook
belowMargin = 50.e3     # z-extent below the source: bottom-truncation safety margin
aboveSurfaceVoid = 5.e3 # void cap above the flat surface (FD/OPT stencil padding)
z = collect((z0 - belowMargin):dx:aboveSurfaceVoid)
nz = length(z)
sourcePosition = (x=0.0, z=z0)
applyFreeSurface = true
material = repeat(reshape(z .<= 0.0, 1, :), nx, 1)
cerjanCells = round(Int, 24 * referenceDx / dx)
# No Cerjan in x at all (there is no lateral wave to absorb); Cerjan only at
# the z-lower (bottom) side, none at z-upper (the physical free surface).
cerjan = CerjanBoundarySpec((0, cerjanCells), (0, 0);
    damping=0.0053 / fdOptSpatialRefinement^2)
bcFD = boundary_geometry(material, (dx, dx); cerjan=cerjan)
duration = quickRun ? 6.0 : 9.0
outputSampling = 0.05
referencePointsPerSWavelength = 16
rickerSharpness = 1.5
minimumAcceptedPointsPerSWavelength = 10
coarsestSpacing = referenceDx
sourceFrequency = rickerSharpness * vs0 / (referencePointsPerSWavelength * coarsestSpacing)
sourceDelay = 1.5 / sourceFrequency
sourceForce = 1.0e10   # per-column force magnitude, N/m (same scale as every
    # other benchmark here); because it is repeated at EVERY column rather
    # than Gaussian-weighted across x, the physical source is an areal
    # force density of sourceForce/dx (N/m^2) -- see analytical.jl.
rickerSource(time) = begin
    a = pi * sourceFrequency * (time - sourceDelay)
    (1 - 2a^2) * exp(-a^2)
end
# Every x column carries the SAME unnormalised weight (not summing to 1
# across x like the Gaussian sources elsewhere): that is what makes this a
# sheet/plane-wave source rather than a localised point source.
sourceNodesX = collect(1:nx)
significantSourceStart = sourceDelay - 1 / sourceFrequency
earliestBottomArrival = significantSourceStart + belowMargin / vp0
@assert duration < earliestBottomArrival "belowMargin is too small for this duration/vp0"
pointsPerSWavelength = (
    FD3=vs0 / (sourceFrequency * dx),
    OPT3=vs0 / (sourceFrequency * (dx / optSpatialRefinement)),
    OPT5=vs0 / (sourceFrequency * (dx / optSpatialRefinement)),
)
@assert minimum(values(pointsPerSWavelength)) >= minimumAcceptedPointsPerSWavelength
# Observation depths: one BELOW the source (direct downward P wave only --
# no free surface involved at all, a baseline sanity check of the solver
# itself) and several BETWEEN the source and the surface (direct + the
# reflected wave once it arrives -- this is what actually tests the
# free-surface closure).
observerDepths = [z0 - 8e3, z0 + 5e3, z0 + 10e3, z0 + 15e3, -3e3]
@assert all(d -> first(z) < d < last(z), observerDepths)
@show quickRun (nx, nz) dx z0 belowMargin duration sourceFrequency sourceDelay
@show pointsPerSWavelength earliestBottomArrival observerDepths
