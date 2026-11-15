# dev/new_Paul_models — figures of Paul's StagYY models

```sh
julia --project=. -t auto dev/new_Paul_models/run.jl                      # all models
julia --project=. -t auto dev/new_Paul_models/run.jl figures G1m20_sb21    # one model
```

Data folder: `ENV["PAUL_MODELS_DIR"]`, default
`/Users/nobuaki/IPGP Dropbox/nobuaki fuji/ForNobu` (one `<model>/Output/` per model).

Output: `outputs/figures_<model>/`
* `<snap>_<field>.png` — 2-D scalar fields: full field, lateral anomaly (minus
  radial mean), radial mean/min/max. Fields spanning > 3 decades (eta, ed, str,
  vd, hpe, plast) are shown in log10.
* `<snap>_vp.png` — v_φ, v_r, |v| (cm/yr) and pressure (GPa).
* `<snap>_cs.png`, `_hf.png`, `_cr.png` — surface fields vs angle, one panel per
  stored z-layer (units assumed: m, W m⁻², m).
* `<snap>_00_overview.png` — every 2-D scalar field of the snapshot on one sheet.

Reading: all through `src/batchFiles/batchStagYY.jl`. `lib/paulStagYY.jl` adds
only what batchStagYY lacks: vector files have an extra scale-factor word after
the header and one ghost point per block (`readField(...; xyghost=1)`), and
surface files (nz ≤ 2). `stages/config.jl`: grid size, seam blending (0 = raw), skipped fields.

## Pseudo electron density

`<snap>_01_electron_density.png`: nₑ = N_A ρ Yₑ from `rho` and `bs`, with
Yₑ = Yₑ,hz + w·(Yₑ,MORB − Yₑ,hz) (w = mass fraction of basalt). End-member
compositions (mol% NCFMAS, from Paul) are in `stages/config.jl`; Yₑ is computed
from them (harzburgite 0.495573, basalt 0.495049). `bs` is a mass fraction.
Panels: Yₑ, nₑ, lateral δln nₑ and δln Yₑ (%), rms of δln ρ, δln Yₑ, δln nₑ vs depth.

## Synthetics and phase windows

```sh
julia --project=. -t auto dev/new_Paul_models/run.jl synthetics    # 6 FD3 runs (cached), sections, snapshots
julia --project=. -t auto dev/new_Paul_models/run.jl phases        # PREM τ–p windows, δt, dlnA, ray maps
```

* Models: PREM (planet1D; crust + ocean replaced by the LID) and
  PREM × (1 + δlnVp, δlnVs, δlnρ) for `G1m20_sb11` 03140 and `G1m20_sb21` 03240;
  δ = lateral anomaly of Paul's Pvel, Svel, rho (field / radial mean − 1, unit-free),
  mantle only, clipped ±10 %. The spherical average of both is PREM itself.
* FD3 on the closed curved free surface of scripts/earth2D_StagYY_neutrino_seismo.jl,
  dx = 10 km, T = 30 s, 2400 s, one dt for all runs.
* Sources: 2000 km depth under the north and the south pole, 45° in-plane thrust
  (Mrr = M0, Mtt = −M0), M0 = 1e15 N m per m, smooth moment step; the south source
  is the north one rotated by 180°. 2-D = line source (cylindrical spreading).
* Stations every 2° all around (φ = 0 north pole, 90 = +x, 180 south pole),
  20 km deep: u_vertical = u·e_r, u_horizontal = u·e_φ (radial SV; no SH in P-SV).
  All seismograms: `outputs/synthetics/synthetics.jld2`.
* Phase windows: PREM first arrival (dev/raytheory τ–p) −25/+45 s, ground velocity,
  cross-correlation delay model − PREM (±15 s) and ln amplitude ratio;
  hollow/faded = another phase predicted in the window. Table:
  `outputs/phases/phase_measurements.csv`. Cache: `data/dev/new_Paul_models/`
  (`--redo=fd3`, `--redo=taup`; bump `codeVersions.jl` after code changes).
