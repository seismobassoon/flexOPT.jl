# dev/freeSurfacePlaneWave/codeVersions.jl
CODE_VERSION = Dict(
    "opt" => "v4-2026-11-16-no-hard-assert-on-blowup",  # v3 + downgraded
        # the "normal"-step stopped_early @assert to a @warn (on
        # request), so compare.jl still renders analytical/FD3/OPT3/
        # OPT5/SPECFEM2D together even though OPT3's free-surface
        # closure is known to blow up around step ~880/1018 -- that
        # truncated, divergent trace is itself the point being shown.
        # Numerically identical to v3 up to the stop point; only the
        # control flow after it changed. v2->v3: pinning the two x-
        # edge columns to the exact analytical solution every step (v2
        # left them frozen at their t=0 IC snapshot, since nothing ever
        # updates them otherwise -- harmless when the field starts/stays
        # at rest, wrong once the IC is nonzero everywhere)
    "fd"  => "v3-2026-11-15-pin-lateral-edges",  # same fix, FD3 side
    "specfem" => "v7-2026-11-16-time-shift-sourceDelay",  # v6
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
