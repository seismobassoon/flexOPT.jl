# helpers/

Claude-authored utility code for the notebooks, kept separate from `src/`
(which Nobu writes and maintains himself). Nothing here replaces or
duplicates `src/` logic -- it only wraps/combines existing `src/` functions
for a specific notebook step. If something here turns out to be broadly
useful, it should graduate into `src/` deliberately, not by accretion.

## Contents

- `sphericalTI_SH.jl` -- 1-D spherical, frequency-domain TI workflow, SH
  (toroidal) Step 1 (notebooks/SphericalTI_SH_oct2026.ipynb): Cartesian TI
  tensor with radial axis, spherical-component divergence, projection onto
  normalised toroidal harmonics, checks, and two `famousEquation` methods
  added from outside src/: `"3DsismoFreqTIradial"` and
  `"1DsphericalSH_TI_Freq"`.

The Kirishima rebuild (notebooks/Elastic2DtopoBenchmark_oct2026.ipynb) still
only calls existing src/ functions directly.
- `sphericalTI_SH_step2.jl` -- Step 2 of the same notebook: planet1D models
  (in-memory DSM-SH zones, PREM), per-zone grids (uniform h per zone by
  default, or planet1D's GT95 `re` grid), scaled local coordinate so one
  Δ=1 recipe serves all zones/h, conventional DSM SH reference matrices,
  and the convergence / zone-by-zone comparison drivers.
- `sphericalTI_SH_step3.jl` -- Step 3: traction-free homogeneous shell,
  exact toroidal eigenfrequencies (j_l, y_l), end-row closures
  (:dsmHalf Galerkin natural row / one-sided traction collocation of order q),
  OPT3 rows for the near-boundary rows wider stencils cannot reach, and the
  generalized eigenproblem A₀v = -Ω²Mv.
- `sphericalTI_SH_step3b.jl` -- Step 3b: time-domain SH synthetics (DSM
  style) for a point force with |m| ≤ 2 normalised toroidal harmonics,
  PREM outer solid region joined at interfaces by adding half-element rows,
  per-scheme global T/H_L/H_N basis matrices, complex-frequency DSM sum over
  l and ω, inverse transform with e^{ω_I t}.
- `sphericalTI_PSV.jl` -- P-SV Step P1 (notebooks/SphericalTI_PSV_oct2026.ipynb):
  spheroidal projection of the TI operator (zonal jet), expected r²-weighted
  radial forms, Cartesian m ≠ 0 check, fluid equation in Q = P/ω
  (Geller & Ohminato 1994), famousEquation entries "1DsphericalPSV_TI_Freq"
  and "1DsphericalFluidQ_Freq".
