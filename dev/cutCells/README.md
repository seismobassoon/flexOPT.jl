# Free-surface topography benchmarks (tilted planes, Gaussian hill)

A self-contained description of the 2-D elastic topography benchmarks used to test the
cut-cell free surface of flexOPT (`src/cutCells/`). It is written so that anyone can
reproduce the cases with **another code** without reading flexOPT: every number needed is
here, and `geometryCases.jl` (plain Julia, no packages) gives the same definitions in
code. Please read §1–§5 before running anything; §6 says how to hand results back.

Status: 25 Nov 2026. Reference traces are being regenerated (see §8).

---------------------------------------------------------------------------------------

## 1. Physical problem (common to every case)

| quantity | value |
|---|---|
| dimension | 2-D plane strain (P-SV), coordinates (x, z), **z positive upward** |
| medium | homogeneous, isotropic, perfectly elastic (no attenuation) |
| P velocity vp | 6000 m/s |
| S velocity vs | 3464 m/s |
| density ρ | 2700 kg/m³ |
| free surface | traction-free, shape depends on the case (§2); **vacuum above** |
| other boundaries | must not influence the recorded window (absorb, or make the domain big enough) |
| output | **displacement** u (m), not velocity |
| duration | 9 s, t = 0 is the start of the source time function (below) |
| reference grid spacing | Δx = Δz = 500 m (flexOPT, SeismicWavePC); SPECFEM2D uses 2 km elements with 5 GLL points |

### Source

* **Force**, not a moment tensor. 2-D line force: total force per unit length
  F = 1.0 × 10¹⁰ N/m, in the direction of the unit vector **f** given per case (§2).
* **Spatially distributed** over the nodes of the 500 m grid with Gaussian weights, so that
  every code sees the same smooth source (a true point force is only acceptable if the
  code cannot do otherwise; say so when you report):

  w_i = exp(−|r_i − r_s|² / (2σ²)) for |r_i − r_s| ≤ 2σ, σ = 500 m, then normalised so that
  Σ w_i = 1. Node i receives the force F w_i **f** s(t).
  The nodes are those of the grid x = k·500 m, z = m·500 m (k, m integers). The case source
  r_s is generally NOT on a node for θ ≠ 0, which is why 12–14 nodes carry weight.
  `caseSourceWeights(case, x, z)` returns the list `(i, k, w)`.
* **Source time function**, a Ricker wavelet with delay:

  s(t) = (1 − 2a²) exp(−a²), a = π f₀ (t − t₀),
  f₀ = 1.5 · vs / (16 · 500 m) = 0.649500 Hz, t₀ = 1.5 / f₀ = 2.309469 s.

  It is dimensionless with peak value 1 at t = t₀. Codes with their own Ricker convention
  (e.g. defined by a rise time, or centred at a different time) must either use this formula
  as an external source file or state their convention exactly.

### Recorded quantities

Displacement at 10 receivers per case: 5 on a **deep** line and 5 on a **shallow** line
(coordinates in §2). For the tilted planes the components are reported in the rotated frame:

u_t = u · **t**, u_n = u · **n**, with **t** = (cos θ, sin θ), **n** = (−sin θ, cos θ),

i.e. tangential and normal to the free surface (`caseRotate`). For the hill, θ = 0, so
u_t = u_x and u_n = u_z.

---------------------------------------------------------------------------------------

## 2. The cases

### 2.1 Tilted planes: tilt0, tilt10, tilt20, tilt30 (θ = 0°, 10°, 20°, 30°)

This is the flat point-force benchmark (free surface z = 0, source 20 km deep, receivers 10 km
and 1 km deep) **rigidly rotated by θ about the origin**: surface, source, force direction and
receivers all rotate together, while the computational grid stays Cartesian. So the free
surface z = x tan θ cuts the grid obliquely.

* Free surface: z_s(x) = x tan θ; material where z ≤ z_s(x)
  (level set φ = **n**·r ≤ 0).
* A flat-frame point (x′, z′) maps to grid coordinates (x, z) = x′ **t** + z′ **n**.
* Source: flat-frame (−10 km, −20 km). Force **f** = **n** (normal to the surface, pointing
  out of the material).
* Receivers: flat frame x′ = −30, −20, −10, 0, 10 km at z′ = −10 km (deep) and z′ = −1 km
  (shallow). The receiver at x′ = −10 km is directly above the source along **n**.

Grid coordinates (km):

| case | source (x, z) | force **f** |
|---|---|---|
| tilt0 | (−10.000, −20.000) | (0, 1) |
| tilt10 | (−6.375, −21.433) | (−0.173648, 0.984808) |
| tilt20 | (−2.557, −22.214) | (−0.342020, 0.939693) |
| tilt30 | (1.340, −22.321) | (−0.500000, 0.866025) |

| case | deep receivers (x, z), x′ = −30 … 10 km |
|---|---|
| tilt0 | (−30, −10) (−20, −10) (−10, −10) (0, −10) (10, −10) |
| tilt10 | (−27.808, −15.058) (−17.960, −13.321) (−8.112, −11.585) (1.736, −9.848) (11.585, −8.112) |
| tilt20 | (−24.771, −19.658) (−15.374, −16.237) (−5.977, −12.817) (3.420, −9.397) (12.817, −5.977) |
| tilt30 | (−20.981, −23.660) (−12.321, −18.660) (−3.660, −13.660) (5.000, −8.660) (13.660, −3.660) |

| case | shallow receivers (x, z), x′ = −30 … 10 km |
|---|---|
| tilt0 | (−30, −1) (−20, −1) (−10, −1) (0, −1) (10, −1) |
| tilt10 | (−29.371, −6.194) (−19.523, −4.458) (−9.674, −2.721) (0.174, −0.985) (10.022, 0.752) |
| tilt20 | (−27.849, −11.200) (−18.452, −7.780) (−9.055, −4.360) (0.342, −0.940) (9.739, 2.481) |
| tilt30 | (−25.481, −15.866) (−16.821, −10.866) (−8.160, −5.866) (0.500, −0.866) (9.160, 4.134) |

**Exact answer.** Because everything rotates together, the tilted solution in (u_t, u_n) is
*identical* to the flat (tilt0) solution in (u_x, u_z) at the same receiver. The reference for
tilt10/20/30 is therefore the **flat** solution: no tilted mesh enters the reference.
The only thing the tilt changes is how the surface sits on the grid, which is exactly what is
tested. Two checks need no reference at all:
* above the source (receiver x′ = −10 km), u_t = 0 by symmetry, for every θ;
* u_n there must not depend on θ.

### 2.2 Gaussian hill: hill

* Free surface: z_s(x) = H exp(−(x − x_c)² / (2w²)), H = 3 km, w = 5 km, x_c = 0
  (w is the Gaussian σ; the hill is between the source and the right-hand receivers).
  Flat (z = 0) far from the hill. Material where z ≤ z_s(x).
* Source: (−10 km, −20 km), force **f** = (0, 1) (vertical, upward), same σ, F and s(t) as above.
* Deep receivers: x = −30, −20, −10, 0, 10 km at z = −10 km.
* Shallow receivers: same x, **1 km below the local surface**, z = z_s(x) − 1 km:
  (−30, −1.000) (−20, −0.999) (−10, −0.594) (0, 2.000) (10, −0.594) km.
* Components: u_t = u_x, u_n = u_z.
* No analytic answer. Reference: SPECFEM2D with the mesh following the hill (§4), checked
  by refining it (2 km → 1 km elements).

---------------------------------------------------------------------------------------

## 3. Comparison rules

* **Raw amplitudes.** No normalisation per trace or per station, and no amplitude or
  time-shift fit. A code with a different force or source-time-function convention is
  allowed **one** global factor and **one** global time lag, measured on tilt0 against the
  reference and applied unchanged to every other case. The values must be reported.
* **Window.** 0 to 9 s, shortened at a receiver if an artificial-boundary return can
  arrive earlier in your setup (flexOPT has no absorbing layer: its window ends at the first
  edge return, computed per receiver).
* **Metric.** Relative L2 misfit with both components together, per receiver:

  e = sqrt(Σ_t (Δu_t² + Δu_n²)) / sqrt(Σ_t (u_t,ref² + u_n,ref²)),

  with the traces resampled on a common 0.01 s axis. Report the mean and max over the 5
  receivers of each line. (Both components go in one norm because u_t vanishes by symmetry
  above the source.)
* **Figures.** Wiggle plots (traces offset by receiver position, time horizontal), one
  scale per receiver line shared by all codes, plus residuals (code − reference) ×5.

---------------------------------------------------------------------------------------

## 4. How each code sees the surface

| code | free surface | notes |
|---|---|---|
| SPECFEM2D | mesh follows the surface exactly (top interface = z_s(x), 4-node elements, 2 km) | spectral elements, the reference |
| flexOPT OPT3 / OPT5 **cut-cell** | Cartesian grid; weak form integrated over the exact material part of every cut cell (`src/cutCells`) | traction-free by the natural BC |
| flexOPT OPT3 / OPT5 **zigzag** | same operators; a cell is material only if all its 4 corners are | staircase |
| SeismicWavePC.jl | staircase: vacuum-like air above z_s(x) (vp = vs = 0, ρ = 1.2) | OpenSWPC-type staggered FD; true point force; own Ricker (rise-time convention); calibrated once on tilt0 (§3) |

At θ = 0 cut-cell and zigzag are the same computation (the surface lies on grid nodes).

---------------------------------------------------------------------------------------

## 5. Pitfalls already met (please check them in your own setup)

1. **Force direction in SPECFEM2D.** With `anglesource = a`, our SPECFEM2D setup applies the
   force (sin a, cos a), so **f** = **n** = (−sin θ, cos θ) needs `anglesource = −θ`. With +θ
   the force is 2θ off the normal. The signature: above the source u_n ∝ cos 2θ
   (0.945/0.766/0.50 of the flat value) and u_t/u_n = 0.8/1.7/3.6 at 10/20/30°. Check the
   symmetry test of §2.1 in any code.
2. **Time origin.** SPECFEM2D labels its seismograms from −1.2/f₀ even with an external
   source time function; the traces must be shifted by +1.2/f₀ to put t = 0 at the start
   of s(t).
3. **SeismicWavePC Ricker.** Its built-in Ricker peaks about 0.8 s before s(t) above; a
   calibration lag search limited to ±0.6 s silently stops at the bound. Prefer an external
   s(t) if the code allows it.
4. **Units.** SeismicWavePC/GeophysicalModelGenerator models are in km and km/s.
5. **Absorbing boundaries on tilted meshes.** The SPECFEM2D runs keep the left/right/bottom
   boundaries far away (bottom 80 km below the surface, measured along **n**) so that their
   quality does not matter within 9 s.

---------------------------------------------------------------------------------------

## 6. Handing results back

One CSV file per case and code, named `<code>_<case>.csv`, header line then one row per sample:

```
line,receiver,x_prime_km,x_km,z_km,time_s,ut_m,un_m
deep,1,-30.0,-27.808,-15.058,0.000,0.0,0.0
...
```

* `line` ∈ {deep, shallow}, `receiver` 1–5 in the order of §2, `x_prime_km` the flat-frame
  position (hill: x), (`x_km`, `z_km`) the grid position actually used.
* Raw displacement in metres, time in seconds from the start of s(t), any sampling ≤ 0.02 s.
* Add a few lines to the hand-back note: grid spacing / element size, time step, how the
  surface and the source were represented, boundary treatment, and any global factor or lag.

`compareAll.jl` reads JLD2 (`outputs/references/<code>_<case>.jld2` with `times` and
`traces["<line>_<receiver>"] = (ut=…, un=…)`). A CSV reader can be added on request, or
write the JLD2 directly with these names.

---------------------------------------------------------------------------------------

## 7. Files in this folder

| file | what it does | environment |
|---|---|---|
| `geometryCases.jl` | **the case definitions** (§1–§2): surface, level set, source and weights, force, receivers, grid, rotation | plain Julia |
| `run.jl`, `env.jl`, `tiltedPlane.jl` | flexOPT runs (`run.jl tiltedPlane tilt0 … hill`) | flexOPT |
| `specfemCases.jl` | SPECFEM2D reference runs and trace extraction | flexOPT + SPECFEM2D in `../othersPackages/specfem2d` |
| `seismicWavePCCases.jl` | SeismicWavePC runs | isolated env `dev/freeSurfacePointSource/seismicWavePCEnv` |
| `compareAll.jl` | inventory, error table (`outputs/compareAll/errors.csv`), wiggle figures | flexOPT |
| `freeEnd1D.jl`, `paper1D.jl` | the 1-D cut-cell free-end benchmark (paper figure) | flexOPT |
| `compareVariants.jl` | older self-comparison of the flexOPT variants (each against its own θ = 0) | flexOPT |

Commands (from the flexOPT root):

```bash
# flexOPT, the four variants (Float64 CPU recipes)
FLEXOPT_CPU_RECIPE=1 julia --project=. dev/cutCells/run.jl tiltedPlane tilt0 tilt10 tilt20 tilt30 hill
FLEXOPT_CPU_RECIPE=1 FLEXOPT_CUT_GEOMETRY=zigzag julia --project=. dev/cutCells/run.jl tiltedPlane tilt0 tilt10 tilt20 tilt30 hill
FLEXOPT_CPU_RECIPE=1 FLEXOPT_CUT_SCHEME=OPT5 julia --project=. dev/cutCells/run.jl tiltedPlane tilt0 tilt10 tilt20 tilt30 hill
FLEXOPT_CPU_RECIPE=1 FLEXOPT_CUT_SCHEME=OPT5 FLEXOPT_CUT_GEOMETRY=zigzag julia --project=. dev/cutCells/run.jl tiltedPlane tilt0 tilt10 tilt20 tilt30 hill
# SPECFEM2D reference (+ optional 1 km refinement)
julia --project=. dev/cutCells/specfemCases.jl
FLEXOPT_SPECFEM_ELEMENT=1000 julia --project=. dev/cutCells/specfemCases.jl tilt0 hill
# SeismicWavePC
julia --project=dev/freeSurfacePointSource/seismicWavePCEnv -t auto dev/cutCells/seismicWavePCCases.jl
# comparison
julia --project=. dev/cutCells/compareAll.jl
```

---------------------------------------------------------------------------------------

## 8. Status (25 Nov 2026)

Raw relative L2 (mean over a line) against SPECFEM2D:

| case | line | OPT3 cut | OPT3 zigzag | OPT5 cut | OPT5 zigzag |
|---|---|---|---|---|---|
| tilt0 (flat) | deep | 0.226 | = cut | 0.080 | = cut |
| tilt0 (flat) | shallow | 0.108 | = cut | 0.039 | = cut |
| hill | deep | 0.219 | 0.252 | – | – |
| hill | shallow | 0.096 | 0.271 | – | – |

* The OPT3 error at tilt0 is its own S-wave dispersion at 500 m (16 points per S wavelength
  at f₀), not the surface; OPT5 removes most of it.
* Against each variant's own flat run, OPT3 cut-cell stays at 0.09–0.18 for 10–30°, the
  same as a control run without free surface (grid anisotropy). OPT3 zigzag reaches 0.30 on
  the shallow line.
* The tilted SPECFEM2D runs made before 25 Nov 15:50 used `anglesource = +θ` (pitfall 1)
  and are being redone; tilt10/20/30 against the flat reference will follow.
* **Zigzag blow-ups are reported as results, not fixed.** OPT5 zigzag at 10° stopped at the
  first time step: the LU bug below. Any zigzag case that still blows up after the fix is
  listed as such in the tables and not tuned further. The cut-cell versions are the ones
  that must run.
* **OPT5 blow-ups fixed (25 Nov, 17:00): they were a linear-algebra bug, not an instability.**
  The equations (rows) were numbered point-major and the unknowns (columns) field-major.
  On these matrices UMFPACK's sparse LU then returns garbage (relative residual 1e13–1e15)
  although the matrix is full rank and well conditioned (sparse QR solves it to 1e-15).
  The fix is in `src/numSolvers/linearSolution.jl`: `prepareConstantMatrix` factorises
  multi-field sparse systems with the rows renumbered (`RowAlignedLU`, same equations),
  for every flexOPT user. `tiltedPlane.jl` prints the LU residual of every run (about 1e-15)
  and stops if it is not small. With the fix, OPT5 cut-cell runs at tilt20, tilt30 and the
  hill. **Results computed before this fix (including OPT5 tilt0) must be recomputed**
  (`--redo=tilt`), since a bad factorization can pollute a run without blowing it up.

---------------------------------------------------------------------------------------

## 9. 3-D cases (26 Nov 2026): oblique planes, then Gaussian hill, then Kirishima

Order agreed: linear slope (oblique planes) → 3-D Gaussian hill → Kirishima.
Definitions in `geometryCases3D.jl` (plain Julia). Same medium, Ricker (f₀ = 0.6495 Hz,
t₀ = 1.5/f₀) and 500 m grid as in 2-D; smaller box for memory.

| quantity | value |
|---|---|
| coordinates | (x, y, z), z up |
| source | flat frame (0, 0, −6 km); Gaussian on the 500 m nodes, σ = 500 m, cutoff 2σ, Σw = 1; point force F = 10¹⁰ **N** (3-D), along **n** |
| receivers (9 per line) | x′ = −5, −2.5, 0, 2.5, 5 km at y′ = 0, then y′ = −5, −2.5, 2.5, 5 km at x′ = 0 |
| lines | deep z′ = −3 km, shallow z′ = −1 km |
| components | (u₁, u₂, u_n) = Rᵀu (flat-frame components) |
| duration | 6 s |
| flexOPT box | ±8 km around the source, 8 km below it, + 10 Cerjan cells (sides, bottom), 3 void cells on top: 53 × 53 × 42–56 nodes |

Cases: `flat3D` (z = 0), `oblique10`, `oblique20`: normal **n** = (−sinθ cosφ, −sinθ sinφ, cosθ),
φ = 30°, so the plane cuts the grid along all three axes. R is the rotation about **e_z** × **n** by θ;
an oblique case is the flat one rotated by R, so the **exact reference is the flat solution**
at the same receiver (as for the 2-D tilts).

**SPECFEM3D_Cartesian** (`specfem3DCases.jl`): internal mesher, one layer between a flat bottom
and the free surface given as the top interface (2 km grid, bilinear: exact for a plane), 2 km
elements, NGLL = 5, Stacey on 5 sides, ±25 km laterally, 25 km below the source,
DT = 2.5 ms; one FORCESOLUTION entry per source node (factor F·w), external STF (SPECFEM3D's
time origin is then 0), stations and sources at their z (`USE_SOURCES_RECEIVERS_Z`).
Check on the old 10-km-offset geometry: the oblique meshes reproduce the flat run rotated to
1–2 % (rel. L2), and u₁, u₂ above the source are < 0.5 % of u_n.
Built serially: `./configure FC=gfortran CC=gcc --without-mpi && make xmeshfem3D
xgenerate_databases xspecfem3D`; `SPECFEM3D_ROOT` (default `../othersPackages/specfem3d`).

**flexOPT 3-D** (`run3D.jl`): OPT3 recipe of `3DsismoTimeIsoHeteroForce` (weak form),
operators with Cerjan padding as `dev/greenland_nov2026`; level set, cut cells and cut rows
on the whole padded grid; the Cerjan weight of past-time columns of padding points is
re-applied to the installed cut rows. Void treatment as in 2-D, with the extension search
allowed up to `reach + 3` steps (`sliverAggregation(...; maxSteps)`): in 3-D a diagonal
neighbour of a cut node can sit ~2.8 cells from the material. Time marching: the
predictor–corrector with lumped-future blocks, except at constraint points (extended and
pinned void nodes), where the block is the node's own coefficients (the lumped block of an
extrapolation row is 1 − Σw = 0). Only the receiver nodes are kept.
Self-check (cut rows of an uncut node vs the volume rows): 1e-14 in 3-D.

Cost (OPT3): volume operator ~7 min for 44 k nodes on 2 cloud cores; cut rows 0.2 s per
cut node after the 26 Nov speed-up of `_cutCellBlockCoefficients` (a function barrier
around the loop that fills S: it was 2 s per node, same arithmetic); about 15 kB per node
for A, L and R, plus copies while rows are installed.

Commands (flexOPT root):

```bash
julia --project=. dev/cutCells/specfem3DCases.jl                       # flat3D oblique10 oblique20
FLEXOPT_CPU_RECIPE=1 julia --project=. -t auto dev/cutCells/run3D.jl   # OPT3 cut-cell, three cases
FLEXOPT_CPU_RECIPE=1 FLEXOPT_CUT_GEOMETRY=zigzag julia --project=. -t auto dev/cutCells/run3D.jl
julia --project=. dev/cutCells/compare3D.jl                            # errors.csv + wiggles3D_<case>.png
```
`run3D.jl --small` runs a 5 km box for 2 s (pipeline test); `--redo` recomputes.
