# dev/freeSurfacePlaneWave/codeVersions.jl
CODE_VERSION = Dict(
    "opt" => "v5-2026-11-24-wider-edge-pin-band",
        # v4->v5: default PLANEWAVE_EDGE_PIN_WIDTH 1->8 (set in
        # run.jl/runMonochromaticSweep.jl). Pinning only the single
        # outermost column each step left too thin a buffer for this
        # narrow domain (nxHalfWidth=10) -- a persistent edge/interior
        # mismatch leaked inward and contaminated the measured center
        # column every step (xSpreadOPT/xSpreadOPT5/xSpreadFD comparable
        # in magnitude to the signal itself). Found 2026-11-24 after the
        # user flagged the "normal" comparison plot as "catastrophic" --
        # that's what led to actually checking xSpreadFD/xSpreadOPT
        # quantitatively instead of just visually; they had been that bad
        # since the very first "nosurface" runs (jobs 005/007), just never
        # looked at numerically. width=8 of nxHalfWidth=10 brought
        # xSpreadOPT/xSpreadOPT5 down to ~1e-4 (essentially exact) and
        # xSpreadFD from ~1.1-1.4 down to ~0.047 in the "nosurface"
        # control (still non-negligible for FD3, but now a normal
        # discretization-level error, not a contamination artifact). Pin
        # width is itself now part of the opt/fd cache key (see opt.jl/
        # fd.jl) so any future width change is a clean cache miss, not a
        # silent collision with a differently-pinned result.
        # v3->v4: downgraded the "normal"-step stopped_early @assert to a
        # @warn (on request), so compare.jl still renders analytical/FD3/
        # OPT3/OPT5/SPECFEM2D together even though OPT3's free-surface
        # closure is known to blow up around step ~880/1018 -- that
        # truncated, divergent trace is itself the point being shown.
        # Numerically identical to v3 up to the stop point; only the
        # control flow after it changed.
        # v2->v3: pinning the two x-edge columns to the exact analytical
        # solution every step (v2 left them frozen at their t=0 IC
        # snapshot, since nothing ever updates them otherwise -- harmless
        # when the field starts/stays at rest, wrong once the IC is
        # nonzero everywhere).
    "fd"  => "v4-2026-11-24-wider-edge-pin-band",  # same v4->v5 OPT fix,
        # FD3 side (same PLANEWAVE_EDGE_PIN_WIDTH default, same cache-key
        # addition). v2->v3 was the same pin-lateral-edges fix described
        # under "opt" above, FD3 side.
    "specfem" => "v8-2026-11-24-read-displacement-directly",
        # v7->v8: src/specfemBenchmark.jl's prepare_specfem2d_case/
        # find_specfem2d_traces gained a `seismotype` keyword (default 2,
        # unchanged for every other caller e.g. dev/kirishima_oct2026).
        # This benchmark now passes seismotype=1 to both and reads
        # SPECFEM2D's native displacement output directly, dropping the
        # trapezoidal velocity->displacement integration entirely -- that
        # integration drifted (peak ~20-40% of the true analytical
        # amplitude, a persistent nonzero post-arrival plateau), even
        # though the time axis was independently verified correct. Found
        # 2026-11-24 when the user asked why SPECFEM2D looked "very bad"
        # once the (unrelated) lateral-edge-leak fix made FD3/OPT3/OPT5
        # look much better by comparison -- SPECFEM2D's error hadn't
        # changed, it just stopped being masked by everything else being
        # equally wrong. See stages/specfem_run.jl.
        # Result (2026-11-24): only PARTIALLY explains the mismatch. The
        # one station farthest from any mesh boundary ("nosurface" z=
        # -28km) improved exactly as predicted (relative_error 1.58->0.12).
        # The other four did not improve in "nosurface" at all (same
        # numbers, bit for bit, as the old integrated-velocity run), and in
        # "normal" several stations now show near-exact PEAK AMPLITUDE
        # (within 1-2% of analytical) yet still 60-98% relative_error --
        # the signature of a timing/phase mismatch, not an amplitude one.
        # So there is at least one more, not-yet-identified SPECFEM2D issue
        # layered on top of the integration-drift one, and the original
        # specfemTimeShift=sourceDelay derivation was only ever
        # independently verified at that single far station, never cross-
        # checked at the other four. User decision (2026-11-24): deprior-
        # itize further SPECFEM2D debugging for now -- FD3/OPT3/OPT5
        # against the exact analytical solution are the trustworthy
        # reference; SPECFEM2D stays a known-rough cross-check until
        # revisited.
        # (shift=1.2/f0) fixed the axis label but left every arrival
        # ~0.46s (=sourceDelay-1.2/f0) too early: SPECFEM2D's analytic
        # plane wave uses a bare, zero-centered Ricker, not one centered
        # at our own sourceDelay. v7 shifts by sourceDelay directly (see
        # the derivation in specfem_run.jl) -- a single global shift,
        # independent of depth, confirmed by matching RickerBare
        # arguments on both sides.
        # v6 (shift=1.2/f0)
        # gave nonzero signal everywhere (confirmed consistent between
        # normal/nosurface at the shared z=-28km receiver, same physics
        # there), but every trace still looked shifted ~1.85s too early.
        # Root cause: SPECFEM2D's own seismogram files start at
        # t=-1.2/f0 (confirmed directly off the raw .semv files), and
        # that shift was only ever wired up for the external-STF path
        # (specfemCase.time_axis_shift), which v3 stopped using. v6
        # applies time_shift=1.2/sourceFrequency unconditionally.
        # v5 (zs=z0)
        # finally ran to completion, but the z=-28km (below z0,
        # direct-wave-only) receiver came back essentially zero: zs=
        # meshTop put almost the whole t=0 field within a few km of the
        # surface. v5 uses zs=z0 instead (our own planeWaveUz reference
        # depth) for both steps.
        # v4 (Stacey not PML, built-in
        # Ricker) still blew up at EXACTLY the same step 5, same max
        # norm (~1.8e12) as v1/v2 -- a dead giveaway the source
        # specification was never the problem. Real cause: the
        # template's PML_BOUNDARY_CONDITIONS=.true. default, which is
        # known-incompatible with initialfield (PML's auxiliary
        # variables start at zero and are never consistent with a
        # globally-nonzero seeded field at the boundary). v4 switches
        # to STACEY_ABSORBING_CONDITIONS instead. v2 (source_factor)
        # fix) still blew up at step 5 even with a correctly-scaled t=0
        # displacement (norm ~1.0) -- the external STF (type 8) seems
        # insufficient for SPECFEM2D to build a consistent initial
        # VELOCITY under initialfield. v3 drops source_time_function and
        # uses the template's own built-in Ricker (time_function_type=1,
        # f0=sourceFrequency) instead -- genuine SPECFEM2D
        # initialfield plane-wave injection (source_type=1/4), homogeneous
        # material, no tomography/sheet-source -- see stages/specfem_run.jl.
        # v1 blew up at time step 5 (max norm ~1.8e12, SIGILL):
        # source_factor defaulted to prepare_specfem2d_case's own point-
        # source amplification (1.0e10), which also multiplies the
        # injected plane-wave field under initialfield. v2 fixes this to
        # source_factor=planewaveAmplitude.
)
PLANEWAVE_CACHE_DIR = joinpath("dev", "freeSurfacePlaneWave")
