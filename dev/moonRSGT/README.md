# dev/moonRSGT — RSGT catalogue for a family of 1-D Moon models

Receiver-side strain Green's tensors (RSGT) for each `.nd` Moon model, so that
students can make synthetics for **any station, any source position down to
the top of the fluid core, any moment tensor or force**, for every model,
without running a solver.

By reciprocity the strain at a source point Q due to a unit force at the
station R gives the station displacement for a moment tensor at Q
(u_n(R) = ε^{R,n}(Q) : M), and the displacement at Q gives it for a force at
Q (impacts).  In a 1-D model all stations are equivalent, so one pole disk
per model is enough.  The solver and the disk are those of
`dev/spherical1Delastic/lib` and `dev/spherical1DelasticKernel/lib`
(`computeSGTDisk(...; family = :force)`, solved with Kᵀ so reciprocity is
exact for OPT3).  Nothing in `src/` is touched.

**The `.nd` models are confidential.** They are read from a folder outside the
repository (`--nd-dir`, default `/gpfs/scratch/nobuaki/moonRSGT/OUTPUT_ND`)
and the catalogue is written outside it too (`--out-dir`, default
`/gpfs/scratch/nobuaki/moonRSGT/catalogue`).  `dev/moonRSGT/.gitignore`
blocks `*.nd`, `*.zip`, `*.jld2` and `outputs/` (the check figures).

## Catalogue (defaults, `moonDefaultParams` in `lib/moonCatalogue.jl`)

| | |
|---|---|
| frequencies | 0 – 0.05 Hz, Δf = 1/3276.8 s (165 frequencies), elastic (no Q), ω_I from omegaiTlen = 1e-2 |
| depths | 0, 20, 40, … km down to the top of the fluid core, + both sides of every discontinuity |
| distances | 0 – 180°, every 0.5° |
| content | 10 strain + 5 displacement pairs, ComplexF32, ≈ 0.5 GB per model (≈ 0.5 TB for 1000 models) |
| solver | OPT3, radial step h = 5 km, lmax(ω) = ωR/2 + 100 |
| interpolation | local Lagrange, order 3 in depth (inside each physical layer, never across the Moho), order 5 in distance |

The grid fixes the band: 0.5° (15 km) × 20 km needs wavelengths of ≳ 4 grid
steps, i.e. f ≲ 0.05 Hz in the crust (vs ≈ 2.9 km/s).  Going higher in
frequency means a finer grid (`--dz`, `--dtheta`) and much more storage.

Models: each physical layer of the `.nd` file is fitted by one cubic in r/R
(the planet1D format), split in halves only where the fit misses a `.nd` value
by more than `fitTol` = 0.005 km/s or g/cm³ (3–13 zones per model).  The TauP
reader of `src/planet1D` is not used: it makes one zone per 10-km line.

## Building it on stella

```sh
# 1. bring the models by hand (not through git)
scp OUTPUT_ND.zip stella:/gpfs/scratch/nobuaki/moonRSGT/       # from your Mac
cd /gpfs/scratch/nobuaki/moonRSGT && unzip OUTPUT_ND.zip && rm -rf __MACOSX
mkdir -p /gpfs/scratch/nobuaki/moonRSGT/logs

# 2. prepare + check on a compute node (S-CAPAD/DANTE rule: never run julia on stella*)
cd /gpfs/users/nobuaki/github/flexOPT.jl
sbatch dev/moonRSGT/slurm/prepare_check.sbatch     # logs in /gpfs/scratch/nobuaki/moonRSGT/logs/
#   or interactively:  srun -p ncpushort -n1 -c1 --mem=8G -t 2:00:00 --pty bash
#                      julia --project=. dev/moonRSGT/run.jl prepare ; julia --project=. dev/moonRSGT/run.jl check

# 3. everything (Slurm array on ncpu, 10 models per task; re-submitting skips finished models)
sbatch dev/moonRSGT/slurm/build_catalogue.sbatch
srun -p ncpushort -n1 -c1 --mem=4G -t 0:30:00 julia --project=. dev/moonRSGT/run.jl info   # what is missing
```

`check` compares, for 3 sources (on the depth grid at 300 km, off-grid at
913 km, off-grid at 9 km in the crust), the catalogue synthetics with an
independent source-side computation at the exact position (no interpolation,
solved with K), and prints the relative L2 differences of spectra and of
traces (Gaussian half-duration 10 s).  The wiggle overlay is saved in
`dev/moonRSGT/outputs/check/` (git-ignored).  Run it before the array.

## Quick calls (SAC/GCMT names)

Julia (on a compute node: `srun -p ncpushort -n1 -c1 --mem=4G -t 1:00:00 --pty bash`, or Jupyterhub):

```julia
include("/gpfs/users/nobuaki/github/flexOPT.jl/dev/moonRSGT/lib/moonSynthetics.jl")
tr = moonSeismogram("00817"; stla=-3.0, stlo=-23.4, evla=20.0, evlo=15.0, evdp=900.0,
                    Mrr=1e15, Mtt=-5e14, Mpp=-5e14, Mrt=0.0, Mrp=0.0, Mtp=0.0)   # GCMT, N·m
tr.t, tr.data[:Z], tr.data[:R], tr.data[:T], tr.Δdeg
# options: comps=(:Z,:N,:E), dt=0.5, halfDuration=15.0, quantity=:velocity
for id in moonModelIds()                                   # all models (catalogue reloaded per model)
    s = moonSeismogram(id; stla=-3.0, stlo=-23.4, evla=20.0, evlo=15.0, evdp=900.0,
                       Mrr=1e15, Mtt=-5e14, Mpp=-5e14, Mrt=0.0, Mrp=0.0, Mtp=0.0)
end
```

Command line / batch (`synth.jl`, one JLD2 = HDF5 output for many models × source-station pairs):

```sh
julia --project=. dev/moonRSGT/synth.jl --models=all --requests=req.csv --out=synth.jld2
# req.csv header: stla,stlo,evla,evlo,evdp,mrr,mtt,mpp,mrt,mrp,mtp  (one line per synthetic)
```

Python (students; numpy + h5py; calls synth.jl, wrapped in `srun` outside a Slurm job):

```python
import sys; sys.path.append("/gpfs/users/nobuaki/github/flexOPT.jl/dev/moonRSGT/python")
from moon_rsgt import moon_synthetics
t, data, models = moon_synthetics("all", stla=-3.0, stlo=-23.4, evla=20.0, evlo=15.0, evdp=900.0,
                                  mrr=1e15, mtt=-5e14, mpp=-5e14, mrt=0, mrp=0, mtp=0)
# data[ireq, imodel, icomp, it], comps "ZRT"
```

## Using it (students)

The reader only needs JLD2 (it includes the geometry and interpolation files of
`dev/spherical1DelasticKernel/lib`, not the solver):

```julia
include("dev/moonRSGT/lib/moonSynthetics.jl")
ENV["MOON_RSGT_DIR"] = "/gpfs/scratch/nobuaki/moonRSGT/catalogue"   # if not the default

cat = loadMoonRSGT(817)                       # rsgt_00817.jld2, ≈ 0.5 GB
cat.model                                     # the .nd table (depth, vp, vs, ρ)
cat.maxDepth                                  # deepest source depth (top of the fluid core)

M  = doubleCouple(strike=30, dip=60, rake=90, M0=1e15)   # N·m
# or M = fromGCMT(Mrr, Mtt, Mpp, Mrt, Mrp, Mtp)
tr = moonSynthetics(cat; station=(-3.0, -23.4), source=(20.0, 15.0, 900.0), M=M,
                    comps=(:Z, :R, :T), dt=1.0, halfDuration=10.0)
tr.t, tr.data[:Z]                              # s, m

# impact: vertical force (N) at the surface, components up, south, east
tr = moonSynthetics(cat; station=(-3.0, -23.4), source=(5.0, 40.0, 0.0), force=(-1e10, 0.0, 0.0))

# the same source and station for every model
trs = moonSyntheticsAllModels(; station=(-3.0, -23.4), source=(20.0, 15.0, 900.0), M=M)
```

Options of `moonSynthetics`: `comps` among `:Z, :N, :E, :R, :T`; `quantity =
:displacement | :velocity | :acceleration`; `stf = :step` (default) or
`:impulse`; `halfDuration` (Gaussian, s; keep ≳ 10 s with fmax = 0.05 Hz).
Conventions: local basis (r, θ, φ) = (up, south, east), i.e. GCMT (r, t, p);
moment tensor in this code's order (Mrr, Mrθ, Mrφ, Mθθ, Mθφ, Mφφ).

## Files

* `run.jl` steps `prepare`, `check`, `build`, `info`
* `env.jl` solver environment without Metal (CPU backend)
* `lib/moonModel.jl` `.nd` → planet1D model (cubic zones)
* `lib/moonCatalogue.jl` radii, systems, `buildMoonRSGT`
* `lib/moonSynthetics.jl` loader, sources, spectra, time series (student side)
* `lib/moonChecks.jl` validation against source-side synthetics
* `slurm/build_catalogue.sbatch` array job
