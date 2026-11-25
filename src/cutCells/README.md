# src/cutCells/

Cut-cell volume integrals for flexOPT: a free surface or an internal
discontinuity that does not follow the Cartesian grid, in 1-D, 2-D and 3-D.
Started 2026-11-24 (rough-topography benchmark). Claude-authored, added with
Nobu's OK; the only change elsewhere in `src/` is one `include` line in
`src/flexOPT.jl`.

## Why

The OPT recipes integrate the test function over rectangles: in
`_constructAmatrix_single` (`src/motorsOPT/makeOPTrecette.jl`) every
coefficient is a product over coordinates of 1-D integrals
(`coefWYYKK[iCoord]`, contracted in `windowContraction!` as
`prod_int *= int[iDim, mᶜ, m, k, l]`), and `testIntegrationBounds` can only
shorten each axis (that is how `clipped_available` treats a flat surface).
A surface that crosses a cell obliquely breaks this product structure, so
stair steps were the only option so far.

Inside one cell every integrand (test function × interpolants × Taylor
kernels) is a polynomial. Its integral over the material part of the cell is
therefore exactly

    Σ_a  c_a · m_a,     m_a = ∫_{cell ∩ material} y^a dy

with `c_a` the coefficients of the product of the 1-D factors. The 1-D
machinery stays as it is; only the final product is replaced by a
contraction with the cut moments `m_a`. That is what this folder provides.

## Conventions

* φ <= 0: side "−" (solid / material 0); φ > 0: side "+" (void / material 1)
* n = ∇φ/|∇φ| points from "−" to "+" (outward normal of a free surface)
* cell `c` spans nodes `c … c+1`; moments are in the cell's local
  coordinates y ∈ [0,1]^N, tensor degree `maxDegree` per axis

## Stage 1 (done 2026-11-24): geometry

| file | content |
|---|---|
| `levelSets.jl` | `CutLevelSet`; from a function, a height (`z = h(x)`, `z = h(x,y)`), or points + normals (`levelSetFromPointsNormals`, `interfaceFromLocalCoords` for GeoPoints `localCoord2D/3D`); multilinear value and normal |
| `simplexIntegration.jl` | Grundmann–Möller rules, exact monomial moments on simplices |
| `cellMoments.jl` | Kuhn decomposition of a cell, exact clipping of each simplex by the (piecewise-planar) interface, `cellMoments`, `CutCellGeometry` (fractions, states, moments of every cut cell) |
| `cutNodes.jl` | `cutNodes` (rows that need a cut operator), `nodeSupportCells`, `stencilCrossings` (η0, η1 and normal on each node → neighbour segment, Ovcharenko 2015 Algorithm 1) |
| `aggregation.jl` | small cut cells (default since 2026-11-24): a void node whose support reaches a sliver of solid is merged into its solid neighbour along the inward normal -- trial value extended (`u_I = 2u_J - u_K`, or `u_J`), cut row added to the target's row (`sliverAggregation`, `applyAggregation`) |

Checks: `julia --project=. dev/cutCells/run.jl geometry`.

## Stage 2 (started 2026-11-24): rows of cut nodes -- `cutCellOperators.jl`

Option A (Nobu, 2026-11-24): a numerical per-node assembler, no edit of the
recipe code.

* `CutCellScheme(params)` takes the same Dict as `makeOPTsemiSymbolic` and
  rebuilds, with the recipe's own functions and caches (`taylorCoefInv`,
  `WYYKKIntegralSymbolic`), the α tables, Cˡη/Cˡηᶜ, the time integrals and,
  for every space axis, the 1-D integrands cell by cell (antiderivatives
  differentiated back, fitted in s = 2y - 1 and checked to be polynomials).
* `cutCellRows(scheme, geom, nodes, models, prepared)` contracts them with the
  cut moments (`CutCellGeometry`, converted to the s-basis) instead of
  multiplying 1-D integrals, evaluates the material terms numerically, and
  returns the rows in the layout of `prepared` (A_unknown / L_known /
  R_force, same columns as `_split_operator_columns`).
* `cutCellSelfCheck`: with every support cell uncut the rows must equal the
  volume rows (guard against drifting from `_constructAmatrix_single`).
* `installRows`, then `applyAggregation` (void slivers), then `pinRows`
  (void nodes with no material in their support).
* Materials must be extended across the interface: the Taylor expansion
  still reads the whole stencil, only the integral is restricted.
* Limits for now: weak form, one ν geometry, no hierarchical test functions,
  1-D integrand breakpoints on grid nodes (true for OPT3 volume recipes).

Additions 2026-11-25:

* `CutCellGeometry(ls; zigzag=true)` (`zigzagLevelSet`): the stair-step
  surface (a cell is material iff all its corners are), through the SAME
  operators, so cut cell and zigzag differ only by the integration domain;
* `sliverAggregation(...; reach, extension)`: every void node within the
  stencil half-width (`reach`: 1 for OPT3, 2 for OPT5) of a cut node is
  extended by Lagrange extrapolation of degree 0-3 (`:constant`,
  `:linear`, `:quadratic`, `:cubic`; linear for OPT3, cubic for OPT5);
* volume recipes built on the CPU in Float64 (`FLEXOPT_CPU_RECIPE=1` in
  dev/cutCells): the self-check is then at round-off (1e-15 in 2-D); the
  August 2026 Metal (Float32) 2-D recipes agree only to ~5e-7.

1-D result (`dev/cutCells/paper1D.jl`, raw traces vs exact): zigzag
reflection late/early by exactly -2ηΔx/c and first-order convergence; cut
cells: timing error ≤ 0.1Δx/c (OPT3) and ≤ 0.01Δx/c (OPT5), raw L2 error at
16 points/λ ≤ 3e-2 (OPT3) and ≤ 3e-3 (OPT5), convergence order ≈ 3.6 (OPT3)
and ≈ 3.0 (OPT5) between 16 and 32 points/λ.

2-D tilted plane (OPT3, Δx = 500 m): the cut-cell surface adds nothing
measurable on top of the grid anisotropy (rotated runs differ from θ = 0 by
the same amount with and without the free surface).

Benchmarks: `dev/cutCells/run.jl freeEnd1D` (1-D bar, free end η Δx beyond
the last node, exact solution; self-check, waveform error, reflection timing
vs the 2ηΔx/c of a stair step, reflected amplitude). Next: 2-D tilted plane
at 0/10/20/30° (reference = validated flat run rotated), Gaussian hill vs
SPECFEM2D, 3-D.

## Stage 3: internal discontinuities

Each side integrated with its own ρ, C (stage 2 with both sides), and the
trial interpolation across the interface made consistent with the jump
conditions: continuity of u, of its tangential derivatives, of the traction
n·C:∇u and the jump condition on (C:∇u)/ρ (Ovcharenko 2015 §3.3.1,
eq. 31–36). The transfer U0 = (B⁻A⁻)⁻¹ A⁺ B⁺ U1 (eq. 37–38) maps the values
beyond the interface; η0, η1 and n come from `stencilCrossings`. The B and A
matrices (eq. 40–42, derived by hand in the thesis for 2-D P-SV) will be
generated from the interface conditions so that 1-D, 2-D and 3-D come from
the same code. First check: 1-D two-layer bar with the interface between two
nodes against the analytical solution (thesis §4.1, which only compared with
a node-aligned run).

## References

* H. Mizutani (2002), PhD thesis, University of Tokyo.
* O. Ovcharenko (2015), Master thesis, IPGP (advisors N. Fuji, R. Martin):
  "An accurate methodology of synthetic seismogram calculation for
  complicated models with regular meshing".
* A. Grundmann, H. M. Möller (1978), Invariant integration formulas for the
  n-simplex by combinatorial methods, SIAM J. Numer. Anal. 15, 282–290.
