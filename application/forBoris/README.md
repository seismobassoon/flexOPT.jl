# forBoris — OPT3 / OPT5 point-force reference on a flat free surface

A clean OPT3/OPT5 reference benchmark meant to be run against an independent
code (SeismicWavePC.jl, or anything else): a single vertical point force in a
homogeneous elastic half-space with a flat free surface at z = 0.

This is a trimmed, non-interactive port of
`notebooks/HomogeneousElastic2DBenchmark_FreeSurface.ipynb` (the same physical
setup that notebook already validates OPT3/OPT5 against FD3 and SPECFEM2D
with). Only OPT3 and OPT5 run here — no SPECFEM2D, no FD3 — so this has no
extra dependency beyond flexOPT itself.

## The physical setup

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

## Run it

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
  panel per receiver.

## Reproducing this in SeismicWavePC.jl

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
   is convention-independent).
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

## Files

* `application/forBoris/run.jl` — thin entry point, forwards to
  `dev/freeSurfacePointSource/run.jl`.
* `dev/freeSurfacePointSource/run.jl` — the actual physics: a trimmed,
  byte-faithful port of cells 3 (config), 8 (OPT3/OPT5 semi-symbolic
  recipes), 9 (operators + free-surface closure + safety audits), 10
  (propagation) and 16 (receiver traces) of
  `notebooks/HomogeneousElastic2DBenchmark_FreeSurface.ipynb`. Only the
  `:clipped_available` free-surface closure and `:direct` time marching are
  kept (the notebook's own defaults); OPT4, the predictor-corrector scheme,
  the `:shifted_full` closure, the Takeuchi & Geller A90 mass-hybrid
  experiment, and the Γ-redistribution ablation are all still in the
  notebook if you want to dig deeper, just not ported here.
