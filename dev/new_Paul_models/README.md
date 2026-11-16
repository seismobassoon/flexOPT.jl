# dev/new_Paul_models — figures of Paul's StagYY models

```sh
julia --project=. -t auto dev/new_Paul_models/run.jl                      # all models
julia --project=. -t auto dev/new_Paul_models/run.jl figures G1m20_sb21    # one model
```

Data folder: `ENV["PAUL_MODELS_DIR"]`, default
`/Users/nobuaki/IPGP Dropbox/nobuaki fuji/PhysChemSelfConsistantModels_for_neutrino_and_seismic` (one `<model>/Output/` per model).

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
* Sources: 200 km depth under the north, east, south and west points (φ = 0, 90, 180, 270°;
  `SEISMO_SOURCE_ANGLES`), 45° in-plane thrust
  (Mrr = M0, Mtt = −M0), M0 = 1e15 N m per m, smooth moment step; every source
  is the north one rotated by its angle. 2-D = line source (cylindrical spreading).
* Stations every 2° all around (φ = 0 north pole, 90 = +x, 180 south pole),
  20 km deep: u_vertical = u·e_r, u_horizontal = u·e_φ (radial SV; no SH in P-SV).
  All seismograms: `outputs/synthetics/synthetics.jld2`.
* Phase windows: PREM first arrival (dev/raytheory τ–p) −25/+45 s, ground velocity,
  cross-correlation delay model − PREM (±15 s) and ln amplitude ratio;
  hollow/faded = another phase predicted in the window. Table:
  `outputs/phases/phase_measurements.csv`. Cache: `data/dev/new_Paul_models/`
  (`--redo=fd3`, `--redo=taup`; bump `codeVersions.jl` after code changes).

Every figure is also copied to `<PAUL_MODELS_DIR>/tell_me_if_you_dont_want_these_files_that_could_take_some_space/<step>/`
(`PAUL_FIGURE_MIRROR`, empty string = off).

## Features: which measurement separates sb11 from sb21?

```sh
julia --project=. -t auto dev/new_Paul_models/run.jl features
```

Windows: the PREM phase windows (−25/+45 s) and fixed 200-s windows over the record.
Features of the ground-velocity segment: ln rms (variance), ln peak, spectral centroid
and bandwidth (mHz), energy duration t90−t10, energy-centroid time, kurtosis,
cross-correlation delay and decorrelation against noise-free PREM, ln H/V.
Coda after each phase window (150 s): number of distinct later arrivals (envelope maxima
> 20 % of the phase peak, ≥ 15 s apart), ln(coda / phase energy), envelope decay rate.
Phase pairs at the same station (`FEATURES.pairs`, e.g. ScS/S, PcP/P, SKS/S, SS/S):
ln amplitude ratio and differential cross-correlation delay (source and site cancel).
Each model trace gets 20 realisations of band-limited Gaussian noise (5 % of the
median PREM trace rms, one level per source and component); per station z = Δ(mean)/sqrt(var₁ + var₂); score = number of
stations with |z| > 2. Outputs `outputs/features/`: `feature_scores.csv` (ranked),
`features_scores_<pole>.png` (feature × window map), `features_00_ranking.png`,
`features_01_top.png` (best 12 vs φ), `feature_values.jld2`. Settings: `FEATURES` in
`stages/config.jl`.

## Neutrinos: 4 surface detectors

```sh
julia --project=. -t auto dev/new_Paul_models/run.jl neutrinos
```

Detectors 2.5 km deep at the source angles (north, east, south, west), atmospheric
νe/νμ/ν̄ (Daemonflux, CC cross sections, perfect detector 20 Mt × 1 yr), both sides
of each detector (rays towards decreasing / increasing φ). nₑ = ρ Yₑ: PREM (Z/A 0.466
core, 0.4955 mantle) and PREM × (1 + δlnρ)(1 + δlnYₑ) for sb11/sb21, δlnYₑ from bs with
Paul's end-members. creationPaths shoots only from the north pole, so the nₑ field is
rotated / mirrored for each detector side (src/ untouched). Outputs `outputs/neutrinos/`:
`nu_00_maps` (δln nₑ, detectors, rays), `nu_01_bars_e|mu` (rms δN/N vs PREM and
√Δχ²(sb11 vs sb21) per detector, per zenith class × energy band), `nu_02_<det>_e|mu`
(δN(E, cosθ) maps), `neutrino_table.csv`. Cache: `nu_*` in data/dev/new_Paul_models.
