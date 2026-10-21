# dev/ — script versions of the notebook sub-projects

Notebooks stay the place to *look* at each step; `dev/<sub_project>/` is the
place to *iterate*: plain `.jl` scripts, run from a terminal or a REPL, with
every expensive stage cached explicitly by `myProduceOrLoad`
(src/batchFiles/batchDrWatson.jl). Claude-authored; `src/` and the
notebooks are not touched.

| sub-project | from notebook | steps |
|---|---|---|
| `kirishima_oct2026/` | `notebooks/Elastic2DtopoBenchmark_oct2026.ipynb` | `step1`, `step2` (cases `homog`, `heter`) |
| `spherical1Delastic/` | `notebooks/SphericalTI_SH_oct2026.ipynb` | `step1`, `step2`, `step3a`, `step3a1`, `step3a2`, `step3b` |

## Running

From the repository root:

```sh
julia --project=. dev/kirishima_oct2026/run.jl step2 heter
julia --project=. dev/kirishima_oct2026/run.jl step2 --redo=opt          # recompute the OPT family only
julia --project=. dev/spherical1Delastic/run.jl step3b --redo=spectra
julia --project=. dev/spherical1Delastic/run.jl all 2>&1 | tee dev/spherical1Delastic/outputs/all.log
```

Fastest for back-and-forth: one REPL **per sub-project** (packages, `src/`
modules and loaded results stay in memory between runs):

```julia
julia> include("dev/kirishima_oct2026/run.jl")        # first run (default step)
julia> kirishima("step2", "heter", "--redo=fd")        # afterwards
julia> spherical("step3a2")                            # (in the spherical REPL)
```

The stage files are re-read on every call, so edits there need nothing.
After editing `helpers/`, add `--reload`; after editing `src/`, restart
Julia (re-including the src/ modules in a live session makes their exported
names ambiguous). All variables of the last run stay available in the REPL
(`uzOPT`, `stationErrorMetrics`, `eigL10p`, `seisTraces`, ...).

Options: `<step>|all`, optional cases, `--redo=<stage>[,<stage>]|all`,
`--full` (Kirishima: `quickRun = false`), `--reload`.

## Caching rules

* Every cached stage is written as
  ```julia
  params = Dict{String,Any}("stage" => ..., "codeVersion" => CODE_VERSION[...], <every parameter>...)
  devRedo!("<stage>", params, CACHE_DIR, prefix)      # only acts with --redo=<stage>
  result = myProduceOrLoad(params, CACHE_DIR, prefix) do cfg ... end
  ```
  so the call site shows what is cached, where and under which parameters.
* Files: `data/dev/<sub_project>/<prefix>_<hash>.jld2` (DrWatson `datadir()`,
  git-ignored). Changing any parameter gives a new hash, i.e. a new file
  (old ones are kept; delete them by hand when you want).
* Changing *code* is not seen by the hash: bump that stage's string in
  `<sub_project>/codeVersions.jl`, or use `--redo=<stage>` once.
* Already-established caches are reused as they are: OPT recipes in
  `data/semiSymbolic` (myProduceOrLoad inside the stages / `shRecipe`), the
  Kirishima velocity models in `notebooks/tmp` (lazyProduceOrLoad, same
  names as the notebook), SPECFEM2D cases in `data/specfem2d_benchmarks`.

## Outputs

Figures (PNG) and small result files go to `dev/<sub_project>/outputs/<step>[_<case>]/`
(git-ignored by `dev/.gitignore`; remove that line to version them).
