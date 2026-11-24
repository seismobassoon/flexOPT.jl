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
#
# PULSE MODE (added 2026-11-24, for the monochromatic frequency-sweep
# benchmark, dev/freeSurfacePlaneWave/runMonochromaticSweep.jl):
# PLANEWAVE_PULSE_MODE[] selects what `pulseSource` (the thing actually fed
# into planeWaveUz) is.
#   :ricker (default) -- EXACTLY the original single broadband Ricker pulse
#       this file always used, same formulas/values as before -- nothing
#       below changes any existing cached opt/fd/specfem result for
#       run.jl's "normal"/"nosurface" steps.
#   :monochromatic -- a tapered, single-carrier-frequency pulse (a
#       Gaussian-enveloped cosine, i.e. a Gabor wavelet) at
#       PLANEWAVE_MONO_FREQUENCY[], PLANEWAVE_MONO_CYCLES[] carrier periods
#       wide. The free-surface closed form (planeWaveUz, two counter-
#       propagating copies of the SAME pulse) does not care what the pulse
#       shape is, so this needs no new derivation -- only duration/
#       belowMargin need to widen with the sweep driver, to fit the
#       (frequency-dependent) burst length within the bottom-truncation
#       safety margin; see the monochromaticMode branch below.
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
aboveSurfaceVoid = 5.e3 # void cap above the flat surface (FD/OPT stencil padding)
applyFreeSurface = DEV_STEP[] != "nosurface"
# Observation depths: one BELOW z0 (direct downward P wave only -- no free
# surface involved at all, a baseline sanity check of the solver itself)
# and several BETWEEN z0 and the surface (direct + the reflected wave once
# it arrives -- this is what actually tests the free-surface closure).
# Moved up (only depends on z0) so the monochromaticMode branch below can
# size its own duration/belowMargin budget off the farthest observer.
observerDepths = [z0 - 8e3, z0 + 5e3, z0 + 10e3, z0 + 15e3, -3e3]

outputSampling = 0.05
referencePointsPerSWavelength = 16
rickerSharpness = 1.5
minimumAcceptedPointsPerSWavelength = 10
coarsestSpacing = referenceDx
rickerFrequency = rickerSharpness * vs0 / (referencePointsPerSWavelength * coarsestSpacing)
rickerDelay = 1.5 / rickerFrequency
rickerSource(time) = begin
    a = pi * rickerFrequency * (time - rickerDelay)
    (1 - 2a^2) * exp(-a^2)
end

if !@isdefined(PLANEWAVE_PULSE_MODE)
    const PLANEWAVE_PULSE_MODE = Ref{Symbol}(:ricker)
end
if !@isdefined(PLANEWAVE_MONO_FREQUENCY)
    const PLANEWAVE_MONO_FREQUENCY = Ref(0.6495)
end
if !@isdefined(PLANEWAVE_MONO_CYCLES)
    const PLANEWAVE_MONO_CYCLES = Ref(3.0)
end
monochromaticMode = PLANEWAVE_PULSE_MODE[] == :monochromatic

if monochromaticMode
    sourceFrequency = PLANEWAVE_MONO_FREQUENCY[]
    monoCycles = PLANEWAVE_MONO_CYCLES[]
    monoSigma = monoCycles / (2sourceFrequency)   # Gaussian envelope std-dev (s)
    monoLeadInSigmas = 2.5   # exp(-2.5^2/2) ~= 4.4% of peak residual at t=0/duration edges
    sourceDelay = monoLeadInSigmas * monoSigma
    monochromaticSource(time) = begin
        tau = time - sourceDelay
        exp(-tau^2 / (2monoSigma^2)) * cos(2pi * sourceFrequency * tau)
    end
    pulseSource = monochromaticSource
    # Duration/belowMargin budget: must (a) let the whole tapered burst play
    # out (~2*sourceDelay, i.e. monoLeadInSigmas on each side of the peak)
    # AND (b) let the wave (direct +/- reflected) reach the farthest
    # observer with margin, (c) while still finishing strictly before the
    # bottom (Cerjan) truncation would be hit -- same safety idea as the
    # Ricker branch's earliestBottomArrival check below, just solved for
    # belowMargin instead of checked against a fixed one, since the
    # required window length now depends on frequency (lower frequency ->
    # longer burst -> needs a deeper bottom margin to stay safe).
    farthestObserverOffset = maximum(abs.(observerDepths .- z0))
    duration = 2monoLeadInSigmas * monoSigma + farthestObserverOffset / vp0 + 1.0
    belowMargin = max(50.e3, 1.3vp0 * duration)
else
    sourceFrequency = rickerFrequency
    sourceDelay = rickerDelay
    pulseSource = rickerSource
    duration = quickRun ? 6.0 : 9.0
    belowMargin = 50.e3   # z-extent below z0: bottom-truncation safety margin
end

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
    pulseSource(time - (zd - z0) / vp0) +
    (includeImage ? pulseSource(time + (zd - z0) / vp0) : 0.0)
)
significantPulseStart = monochromaticMode ? 0.0 : (sourceDelay - 1 / sourceFrequency)
earliestBottomArrival = significantPulseStart + belowMargin / vp0
@assert duration < earliestBottomArrival "belowMargin is too small for this duration/vp0"
pointsPerSWavelength = (
    FD3=vs0 / (sourceFrequency * dx),
    OPT3=vs0 / (sourceFrequency * (dx / optSpatialRefinement)),
    OPT5=vs0 / (sourceFrequency * (dx / optSpatialRefinement)),
)
if monochromaticMode
    minimum(values(pointsPerSWavelength)) >= minimumAcceptedPointsPerSWavelength || @warn(
        "below the usual minimum points-per-S-wavelength -- expected/intentional for a " *
        "frequency sweep probing past the resolution limit",
        pointsPerSWavelength, minimumAcceptedPointsPerSWavelength)
else
    @assert minimum(values(pointsPerSWavelength)) >= minimumAcceptedPointsPerSWavelength
end
@assert all(d -> first(z) < d < last(z), observerDepths)
@show quickRun (nx, nz) dx z0 belowMargin duration sourceFrequency sourceDelay planewaveAmplitude monochromaticMode
@show pointsPerSWavelength earliestBottomArrival observerDepths
