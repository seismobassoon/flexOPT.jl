# dev/japanTrench_nov2026 — Japan Trench box (towards a Vs / Qs inversion)

Region of Yamaya, Kubo, Shiomi & Takemura (2024, *JGR Solid Earth* 129,
e2024JB029944; S-net + F-net CMT with 3-D Green's functions), but a small
box first: centre **37°N 142°E**, **±200 km east and north**, down to
**80 km below sea level**. Claude-authored; `src/` is not modified.

```sh
julia --project=. -t auto dev/japanTrench_nov2026/run.jl step1            # quick, dx = 5 km
julia --project=. -t auto dev/japanTrench_nov2026/run.jl step1 standard   # dx = 2.5 km
julia --project=. -t auto dev/japanTrench_nov2026/run.jl step1 full       # dx = 1 km (model only)
julia --project=. -t auto dev/japanTrench_nov2026/run.jl step2 coarse     # plane wave, OPT3 3-D, dx = 8 km
julia --project=. -t auto dev/japanTrench_nov2026/run.jl step3 coarse     # fences + wiggle record sections
julia --project=. -t auto dev/japanTrench_nov2026/run.jl step4 coarse     # movie (~20 s)
julia --project=. -t auto dev/japanTrench_nov2026/run.jl all quick        # everything at dx = 5 km
```
Grid cases: `coarse` 8 km (Cerjan 10 cells), `quick` 5 km (12), `standard` 2.5 km (16), `full` 1 km (20).

## Model (step1)

* **Box**: `GeoPoints.constructLocalBox`, `axis_angle_deg = 90` → x east,
  y north, z up. True Cartesian frame: sea level is ~6 km below the tangent
  plane at the corners, so `altMin` is lowered by that much and every column
  has ≥ 80 km below sea level. Top = 3 km + one air cell.
* **Material**: `getParamsAndTopo(...; velocity_model=:NIED)`, exactly the
  Kirishima recipe: planet1D PREM by radius (ρ, Vp, Vs, **Qμ, Qκ**), NIED
  Vp/Vs where confidence ≤ 0.8 (otherwise PREM), GMT `earth_relief_15s`
  topography/bathymetry, sea water (1.03 g/cm³, 1.5 km/s) between sea floor
  and sea level, air above. Point sampling (no sub-cell mixing yet).
* Solid Vs raised to **1.0 km/s** as in the paper (§3); Vp kept.
* Cached: `data/dev/japanTrench_nov2026/model_<hash>.jld2` (`--redo=model`).

Figures (`outputs/step1_<case>/`): `model_map` (solid-top elevation, top Vs,
NIED coverage), `model_depth_slices` (Vp, Vs, NIED mask at 5–70 km),
`model_section_EW` / `model_section_NS` (Vp, Vs, ρ, Qμ through the centre).

## Ocean: currently AIR in step1 (src bug), used as a load in step2

`getParamsAndTopo` (src/GeoPoints/getSeismicParamTopo.jl) tests
`topography < alt → air` **before** `topography ≤ alt ≤ 0 → water`, so every
point between the sea floor and sea level is air: no water cell exists, the
step1 sections show air down to the sea floor, and the circle on `model_map`
is sea level stepping through the cells (land < 2 km is hidden at 5 km).
Fix (src, not done): test air only for `alt > max(topography, 0)`.

For step2 this is what "ocean as a burden" needs: the rock top of an oceanic
column is the sea floor and carries the free surface. `ocean_load.jl` reads the
water depth h from the same GMT grid (cache `data/topoViaGMT`) and adds
ρw h / dx to ρ of the top rock cell (λ, μ unchanged). One ρ per node acts on
all three components: the horizontal inertia of the sea floor is overestimated
(an inviscid ocean loads the vertical motion only). `OCEAN_LOAD = :none` turns
it off.

## Plane wave (step2)

* OPT3 3-D, same recipe file as Greenland (`data/semiSymbolic`), Cerjan on the
  four sides and the bottom, `apply_free_surface_3d` (helpers) on the rock top
  (land surface + sea floor, staircase, not benchmarked in 3-D yet),
  lumped-future predictor–corrector.
* Source = smooth horizontal sheet at z = −65 km of moment density
  m = M w(x,y) G(z) r(t − p·x_h): f = −∇m for **P** (curl-free), f = ∇m × a for
  SV/SH (div-free). p = sin 45° / V(sheet layer), Ricker 20 s, Tukey taper
  60 km inside the box, Gaussian of half-width 1.5 dx in z. The sheet also
  sends a down-going wave into the bottom Cerjan layer. Default: P, towards
  azimuth 270° (coming from the trench side). Change in `config.jl`.
* Stored: 3 components on every rock-top node every 0.5 s; E–W (y = 0) and
  N–S (x = 0) sections every 2 s. Cache `opt_<case>_<hash>.jld2`, `--redo=opt`.
* step3: `ocean_load`, `fence_uz`, `fence_ux` (3-D fences at 6 times),
  `record_EW`, `record_NS` (wiggles, x/y/z, red = imposed p). step4:
  `planewave_<case>.mp4`.

## Notes for the next steps

* Grid vs period (OPT3, ~4 points per S wavelength at Vs = 1 km/s):
  quick 5 km → ≳ 20 s (the paper's 20–50 s band), standard 2.5 km → ≳ 10 s.
* Qs starting model = PREM Qμ (600 in crust and lid; 80 below 80 km depth,
  which the bottom cells near the centre reach). NIED has no Q.
* ρ comes from PREM, not from NIED Vp (a Vp–ρ relation could be added).
* The ocean is up to ~7.5 km deep near the trench: at 5 km it is 0–2 cells
  of Vs = 0. μ = 0 cells are untested in the elastic OPT operator.
* No free surface / Cerjan / source yet; NIED offshore coverage is shown
  in `model_map` and the NIED rows of `model_depth_slices`.
