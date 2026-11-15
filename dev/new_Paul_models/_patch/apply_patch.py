# run from dev/new_Paul_models:  python3 apply_patch.py <dir with new files>
import os, sys, re, shutil
src = sys.argv[1]
def rd(p): return open(p).read()
def wr(p, s): open(p, 'w').write(s)
def rep(p, old, new, count=1):
    s = rd(p)
    assert old in s, (p, old[:60])
    wr(p, s.replace(old, new, count))

# new files
shutil.copy(os.path.join(src, 'lib/paulFeatures.jl'), 'lib/paulFeatures.jl')
shutil.copy(os.path.join(src, 'stages/features.jl'), 'stages/features.jl')

# A. uninitialised uz
rep('lib/paulSeismo.jl', 'ux = zeros(Float32, nstore, nrec); uz = similar(ux)',
    'ux = zeros(Float32, nstore, nrec); uz = zeros(Float32, nstore, nrec)')

# B. synthetics: sanitise first sample of cached runs; paulSave
rep('stages/synthetics.jl', '''            label = "$name $pole")
    end
end''', '''            label = "$name $pole")
    end
    # runs cached before the 2026-11-15 fix have an uninitialised first sample
    # (uz = similar(ux)); u(t = 0) = 0 exactly, so set it instead of recomputing.
    for key in ("uVert", "uHoriz")
        seismoResults[(name, pole)][key][1, :] .= 0
    end
end''')
for p in ('stages/synthetics.jl', 'stages/phases.jl', 'stages/figures.jl'):
    s = rd(p); wr(p, s.replace('devSave(', 'paulSave('))

# F. wiggle record sections + stronger gathers
s = rd('lib/paulSeismoFigures.jl')
a = s.index('function recordSectionFigure(')
b = s.index('function snapshotFigure(')
s = s[:a] + rd(os.path.join(src, 'lib/wiggle_section.jl')) + '\n' + s[b:]
old = 'append!(ysAll, φ .+ 1.6 .* V[idx, j] ./ s); push!(ysAll, NaN)'
assert old in s
s = s.replace(old, 'append!(ysAll, φ .+ 3.0 .* clamp.(V[idx, j] ./ s, -1.3, 1.3)); push!(ysAll, NaN)')
wr('lib/paulSeismoFigures.jl', s)

# C. config
s = rd('stages/config.jl')
s += '''
# record sections (wiggles): every `every`-th station, trace / max|PREM trace| × gain, clipped
SEISMO_SECTION = (gain = 3.0, clip = 1.2, every = 1)

# ---- features (step `features`): which measurement separates sb11 from sb21?
FEATURES = (
    phases = SEISMO_PHASES,     # PREM-predicted phase windows [t − pre, t + post]
    pre = 25.0, post = 45.0,
    slideStart = 200.0, slideStep = 200.0, slideLength = 200.0,  # fixed windows over the record
    maxLag = 15.0,              # s, cross-correlation search
    noiseLevel = 0.05,          # noise rms / rms of the whole PREM trace (station, component)
    nReal = 20,                 # noise realisations per trace
    zThreshold = 2.0,           # |z| > 2: station distinguishes the two models
    seed = 20261115,
)
'''
wr('stages/config.jl', s)

# D. run.jl
s = rd('run.jl')
s = s.replace('get(ENV, "PAUL_MODELS_DIR", "/Users/nobuaki/IPGP Dropbox/nobuaki fuji/ForNobu")',
    'get(ENV, "PAUL_MODELS_DIR",\n    "/Users/nobuaki/IPGP Dropbox/nobuaki fuji/PhysChemSelfConsistantModels_for_neutrino_and_seismic")')
s = s.replace('steps = ["figures", "synthetics", "phases"]', 'steps = ["figures", "synthetics", "phases", "features"]')
s = s.replace('''include(joinpath(PAUL_STAGES, "config.jl"))
''', '''include(joinpath(PAUL_STAGES, "config.jl"))

# every figure also goes to the shared Dropbox folder (empty string: off)
PAUL_FIGURE_MIRROR = get(ENV, "PAUL_FIGURE_MIRROR",
    joinpath(PAUL_MODELS_DIR, "tell_me_if_you_dont_want_these_files_that_could_take_some_space"))
"devSave + copy to PAUL_FIGURE_MIRROR/<step>/"
function paulSave(fig, name)
    path = devSave(fig, name)
    if !isempty(PAUL_FIGURE_MIRROR)
        dir = joinpath(PAUL_FIGURE_MIRROR, DEV_STEP[])
        mkpath(dir)
        cp(path, joinpath(dir, basename(path)); force = true)
    end
    return path
end
''', 1)
s = s.replace('''        include(joinpath(@__DIR__, "lib", "paulSeismoFigures.jl"))
''', '''        include(joinpath(@__DIR__, "lib", "paulSeismoFigures.jl"))
        include(joinpath(@__DIR__, "lib", "paulFeatures.jl"))
''')
s = s.replace('''        step == "phases" && @time "phases" include(joinpath(PAUL_STAGES, "phases.jl"))''',
'''        step == "phases" && @time "phases" include(joinpath(PAUL_STAGES, "phases.jl"))
        step == "features" && @time "features" include(joinpath(PAUL_STAGES, "features.jl"))''')
s = s.replace('#   julia --project=. -t auto dev/new_Paul_models/run.jl phases             # phase windows (after synthetics)',
 '#   julia --project=. -t auto dev/new_Paul_models/run.jl phases             # phase windows (after synthetics)\n#   julia --project=. -t auto dev/new_Paul_models/run.jl features           # discriminating features sb11 vs sb21')
wr('run.jl', s)
for k in ('paulSave', '"features"', 'PhysChem'):
    assert k in s, k

# env: Random, DSP for features
rep('stages/seismo_env.jl', 'using Interpolations, JLD2, LinearAlgebra, DelimitedFiles',
    'using Interpolations, JLD2, LinearAlgebra, DelimitedFiles, Random, DSP')

# README
s = rd('README.md')
s = s.replace('`/Users/nobuaki/IPGP Dropbox/nobuaki fuji/ForNobu`', '`/Users/nobuaki/IPGP Dropbox/nobuaki fuji/PhysChemSelfConsistantModels_for_neutrino_and_seismic`')
s += '''
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
Each model trace gets 20 realisations of band-limited Gaussian noise (5 % of the
median PREM trace rms, one level per source and component); per station z = Δ(mean)/sqrt(var₁ + var₂); score = number of
stations with |z| > 2. Outputs `outputs/features/`: `feature_scores.csv` (ranked),
`features_scores_<pole>.png` (feature × window map), `features_00_ranking.png`,
`features_01_top.png` (best 12 vs φ), `feature_values.jld2`. Settings: `FEATURES` in
`stages/config.jl`.
'''
wr('README.md', s)
print('patch applied')
