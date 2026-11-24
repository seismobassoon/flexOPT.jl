# forBoris — OPT3 / OPT5 reference benchmarks: point force and plane wave, flat free surface

Two clean OPT3/OPT5 reference benchmarks meant to be run against an
independent code (SeismicWavePC.jl, or anything else), both with a flat free
surface at z = 0 in a homogeneous elastic half-space:

1. **Point force** (`run.jl`) — a single vertical point force embedded in the
   medium. Needs source-amplitude/source-time-function conventions to be
   matched between codes (see the caveats below).
2. **Plane wave** (`runPlaneWave.jl`) — a normal-incidence plane P wave,
   injected as an *exact closed-form initial condition*, no body force
   anywhere. There is **no amplitude or source-time-function convention to
   match at all** between codes for this one — just the same model
   (vp/vs/rho), the same initial field, and a free surface at z = 0. This
   also carries a known, currently-open issue in OPT3's free-surface closure
   (see below) that an independent cross-check is specifically useful for.

Both are trimmed, non-interactive ports of this repository's own validated
benchmarks: only OPT3 and OPT5 run here — no SPECFEM2D, no FD3 — so neither
has any extra dependency beyond flexOPT itself.

## 1. Point-force benchmark

### The physical setup

| quantity | value | units |
|---|---|---|
| P velocity `vp0` | 6000 | m/s |
| S velocity `vs0` | 3464 | m/s |
| density `rho0` | 2700 | kg/m³ |
| source position | x = −10, z = −20 | km (z is depth below the free surface, which sits at z = 0) |
| source | vertical force only, Gaussian-smoothed in space, σ = 500 m, cutoff 2σ | — |
| total force | `sourceForce = 1.0e10` | N/m (a 2-D line force: total vertical force per unit length in the invariant out-of-plane direction) |
| source time function | Ricker, peak frequency ≈ 0.65 Hz (`1.5·vs0/(16·500 m)`), delay 1.5/frequency | — |
| receivers | horizontal line at z = −10 km, x = −30, −20, −10, 0, 10 km (quick) | km |
| grid spacing `dx` | 500 | m |
| domain half-width | 60 km (quick) / 90 km (`--full`) around the source | — |
| duration | 9 s (quick) / 14 s (`--full`) | — |

The domain/duration are sized so no artificial (side or bottom) boundary
return reaches a receiver within the recorded window — there is no absorbing
layer, the truncation is simply never seen. The free surface itself is
physical, not absorbing.

Output quantity is **displacement** (u_x, u_z), not velocity.

**On "incidence angle":** this is a point source, not a plane wave, so there
is no single incidence angle -- each receiver sits at its own angle from
the source. Measured from vertical (0 deg = straight down from source to
receiver, i.e. no lateral offset):

| receiver x (km) | horizontal offset from source (km) | vertical offset (km) | angle from vertical |
|---|---|---|---|
| -30 | 20 | 10 | 63.4 deg |
| -20 | 10 | 10 | 45.0 deg |
| -10 | 0  | 10 | 0.0 deg (directly above the source) |
| 0   | 10 | 10 | 45.0 deg |
| 10  | 20 | 10 | 63.4 deg |

(source at x=-10 km, z=-20 km; receivers at z=-10 km, i.e. 10 km above the
source -- angle = atan(|horizontal offset| / 10 km).)

The wavelet is a Ricker pulse, `(1-2a^2)*exp(-a^2)` with
`a = pi*f0*(t-delay)`, `f0 ~= 0.6495 Hz`, `delay = 1.5/f0 ~= 2.31 s`, applied
as a vertical force Gaussian-smoothed in space (sigma = 500 m, cutoff 2 sigma)
around the source position -- all of this, plus the model (vp0/vs0/rho0), is
now also printed directly on `opt3_vs_opt5_waveforms.png` itself, not only
here, after an earlier version of this plot proved hard to reproduce from
the image alone.

### Run it

```sh
julia --project=. application/forBoris/run.jl run           # quick: 60 km half-width, 9 s
julia --project=. application/forBoris/run.jl run --full    # fuller: 90 km half-width, 14 s
```

First run compiles and (if needed) builds the OPT3/OPT5 semi-symbolic recipes
(cached under `data/semiSymbolic/`); on this repository those recipes are
already cached from the notebook's own runs, so it should go straight to
propagation. A few safety checks run before any propagation (coefficient
health gate, free-surface factorization audit for both OPT3 and OPT5) and
will fail loudly, before anything expensive runs, if something is off.

### Output

`application/forBoris/outputs/run/`:
* `forBoris_opt3_opt5_traces.jld2` — `receiverX`, `receiverZ`, `dx`, `dxOPT`,
  `dtOPT`, the source parameters, and `opt3DisplacementX/Z`,
  `opt5DisplacementX/Z`, each `(time=…, values=…)` with `values` a
  `(nt, nreceivers)` matrix, receivers in the order of `receiverX`.
* `forBoris_opt3_opt5_traces.csv` — the same traces, one column per
  receiver/component/solver, for anything that isn't Julia.
* `opt3_vs_opt5_waveforms.png` — u_z wiggle record section, OPT3 vs OPT5, one
  panel per receiver, with the model/source/wavelet printed on the figure
  itself.

### Reproducing this in SeismicWavePC.jl

The mapping below is what we worked out comparing this benchmark's
conventions against `SeismicWavePC.jl`'s API (`ForceSource`, `Station`,
`simulate`). Coordinates there are also (x, y, z) in km with z **up** — same
sign convention as here, so no flip is needed, only a km/m unit change:

```julia
using SeismicWavePC, GeophysicalModelGenerator

# half-space + flat free surface: air above z = 0, solid below
xkm = collect(-70:0.5:50)     # km; matches domainHalfWidth=60e3 around sourcePosition.x=-10e3
zkm = collect(-80:0.5:40)     # km; matches domainHalfWidth=60e3 around sourcePosition.z=-20e3
model = CartData(xyz_grid(xkm, 0, zkm))
X, Y, Z = xyz_grid(xkm, 0, zkm)
isAir = Z .> 0.0
model = addfield(model, (; rho=ifelse.(isAir, 1.2, 2700.0),
                           vp=ifelse.(isAir, 0.0, 6.0),      # km/s
                           vs=ifelse.(isAir, 0.0, 3.464)))   # km/s

source = ForceSource(-10.0, 0.0, -20.0; fz=1.0, t0=0.0, trise=1/0.6495)  # see caveat below
stations = [Station("S$i", xr, 0.0, -10.0) for (i, xr) in enumerate([-30,-20,-10,0,10])]

result = simulate(model; dt=..., nt=..., sources=source, stf=:ricker, stations, backend=CPU())
uz = result.u[:, 2, :]   # displacement, components are (:x, :z) in 2-D P-SV
```

**Two conventions do not translate one-to-one, and matter before comparing
absolute amplitudes or exact waveform shape (timing and relative amplitude
across receivers should already agree):**

1. **Force amplitude.** This benchmark's `sourceForce = 1.0e10 N/m` is the
   literal, physical total vertical line force. SeismicWavePC.jl inherits
   OpenSWPC's own source normalization (see `simulate`'s `M0` bookkeeping in
   its `src/simulation.jl`): a `ForceSource`'s `fz` sets direction, and the
   actual applied magnitude is carried by the source-time-function
   normalization instead, not a direct Newton import. If you want the exact
   same physical force rather than a direction, either reverse-engineer that
   normalization from `simulate`'s source code, or simplest: run both codes
   once, then rescale one trace set by the ratio of peak amplitudes at a
   receiver away from the free surface (so the comparison is about waveform
   shape/timing/relative-amplitude-across-receivers, which is the part that
   is convention-independent). (This is exactly what the flexOPT-side
   comparison script does when folding SeismicWavePC.jl traces into its own
   figures: a per-station peak-amplitude rescale plus a cross-correlation
   time-alignment, documented in dev/freeSurfacePointSource/runWithReferences.jl.)
2. **Source time function.** This benchmark's `rickerSource(t)` is
   `(1-2a²)e^{-a²}` with `a = π·f·(t - delay)`, `delay = 1.5/f`. SeismicWavePC's
   built-in `:ricker` is the same functional form but centred and scaled by
   its own `trise` convention (peak frequency `1/trise`, centred at
   `t0 + trise`) — close but not byte-identical in general. For an exact match,
   pass a `SampledSTF` built directly from this benchmark's `rickerSource`
   (saved time axis is in the JLD2 output above) instead of `stf=:ricker`.

If you'd rather skip reproducing the model in SeismicWavePC.jl yourself and
just get numbers to diff against, the CSV/JLD2 above already has everything
(frequency, delay, force, geometry) needed to build the comparison from the
other side.

## 2. Plane-wave benchmark

### The physical setup

| quantity | value | units |
|---|---|---|
| P velocity `vp0` | 6000 | m/s |
| S velocity `vs0` | 3464 | m/s |
| density `rho0` | 2700 | kg/m³ |
| incidence | **0° (normal incidence)** — a true plane wave, purely vertical and x-independent by construction | — |
| source | **none** — the wave is injected as an exact closed-form initial condition (two counter-propagating copies of the same pulse); see the derivation in `dev/freeSurfacePlaneWave/stages/analytical.jl` | — |
| wavelet | Ricker, same formula as the point-force benchmark: `(1-2a²)e^{-a²}`, `a = π·f0·(t-delay)`, `f0 ≈ 0.6495 Hz`, `delay = 1.5/f0 ≈ 2.31 s` | — |
| free surface | flat, at z = 0 (void above; "normal" step) — compared against a no-free-surface control ("nosurface" step) | — |
| observation depths | z0−8, z0+5, z0+10, z0+15, −3 km, where z0 = −20 km (a travel-time phase reference only — there is no source location) | km |
| grid spacing `dx` | 250 m (`dx/2`, finer than the point-force case) | — |
| duration | 6 s (quick) / 9 s (`--full`) | — |

Because the field never depends on x by construction, there is no lateral
wave to reflect off anything — the only truncation that matters is the
bottom one, kept duration-safe the same way as the point-force benchmark.
Both "normal" (actual free surface) and "nosurface" (control: fully solid
medium, no closure exercised at all) always run together, in one invocation.

**Why no source-amplitude or source-time-function convention to match**:
unlike the point-force case, this benchmark needs no body force, no
Green's-function calibration, and no areal-force-density/dx-normalization
convention at all. Reproducing it in any other code is just: build the same
homogeneous model, seed a purely vertical, x-independent `u_z(z, t=0)` field
using the SAME formula (below), put a flat free surface at z = 0, and
propagate with zero forcing from there on.

```
u_z(z, t) = pulse(t - (z - z0)/vp)  +  pulse(t + (z - z0)/vp)     # "normal" (free surface)
u_z(z, t) = pulse(t - (z - z0)/vp)                                 # "nosurface" (control; direct wave only)
pulse(t) = (1 - 2a²)·exp(-a²),   a = π·f0·(t - delay)
```

(the two-counter-propagating-copies form is the plane-wave analogue of the
method-of-images mirror source: it satisfies the zero-traction free-surface
condition *exactly*, for any pulse shape, by construction — see
`stages/analytical.jl`'s full derivation.)

**KNOWN OPEN ISSUE, surfaced here deliberately, not hidden:** OPT3's
free-surface closure (an additive `elasticTractionFree2D` construction) is
known, as of this writing, to **blow up partway through the "normal" run**
(its trace is truncated/divergent in the output figure, not hidden or
asserted away). The "nosurface" control has never shown this. **This is
exactly the kind of thing an independent code is useful for**: if your own
normal-incidence response does *not* diverge like this, that's evidence the
issue is specific to this closure rather than to the physical setup — please
let us know either way.

### Run it

```sh
julia --project=. application/forBoris/runPlaneWave.jl run           # quick
julia --project=. application/forBoris/runPlaneWave.jl run --full    # fuller
```

Always runs both "normal" and "nosurface" in one invocation (not
user-selectable here). Reuses the same OPT3/OPT5 semi-symbolic/operator
caches as the point-force benchmark and as this repository's own
`dev/freeSurfacePlaneWave/run.jl`.

### Output

`dev/freeSurfacePlaneWave/outputs/normal/` and `.../nosurface/`:
* `forBoris_vs_analytical_uz.png` — u_z, OPT3 and OPT5 vs the exact
  analytical solution, one panel per observation depth, model/incidence/
  wavelet/known-issue printed on the figure itself. The y-axis is clipped to
  3× the analytical peak at each depth, so a divergent OPT3 trace is still
  visibly cut off rather than compressing everything else to a flat line.
* `forBoris_vs_analytical_metrics.jld2` — `observerDepths` and
  `forBorisMetrics` (per-depth `max_abs_error`/`relative_error`/
  `max_analytical`/`max_numerical` for OPT3 and OPT5 against the exact
  closed form).

## Files

* `application/forBoris/run.jl` — thin entry point, forwards to
  `dev/freeSurfacePointSource/run.jl` (point-force benchmark).
* `application/forBoris/runPlaneWave.jl` — thin entry point, forwards to
  `dev/freeSurfacePlaneWave/forBorisRun.jl` (plane-wave benchmark).
* `dev/freeSurfacePointSource/run.jl` — the point-force physics: a trimmed,
  byte-faithful port of cells 3 (config), 8 (OPT3/OPT5 semi-symbolic
  recipes), 9 (operators + free-surface closure + safety audits), 10
  (propagation) and 16 (receiver traces) of
  `notebooks/HomogeneousElastic2DBenchmark_FreeSurface.ipynb`. Only the
  `:clipped_available` free-surface closure and `:direct` time marching are
  kept (the notebook's own defaults); OPT4, the predictor-corrector scheme,
  the `:shifted_full` closure, the Takeuchi & Geller A90 mass-hybrid
  experiment, and the Γ-redistribution ablation are all still in the
  notebook if you want to dig deeper, just not ported here.
* `dev/freeSurfacePlaneWave/forBorisRun.jl` — the plane-wave physics: reuses
  `dev/freeSurfacePlaneWave/stages/{config,opt,opt_build,opt_propagate,
  analytical}.jl` unchanged (FD3 and SPECFEM2D stages dropped), with its own
  compare/plot section (OPT3/OPT5 vs analytical only, captioned for
  reproducibility). `dev/freeSurfacePlaneWave/run.jl` is the fuller version
  with FD3/SPECFEM2D cross-checks and the frequency-sweep diagnostic, kept
  for this repository's own ongoing debugging of the free-surface-closure
  issue above.
