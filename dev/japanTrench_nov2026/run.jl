# dev/japanTrench_nov2026/run.jl -- Japan Trench box for a Vs / Qs inversion
#
# Box of Yamaya, Kubo, Shiomi & Takemura (2024, JGR 129, e2024JB029944) style
# region, but small: centre 37°N 142°E, ±200 km east and north, down to 80 km
# below sea level; NIED 3-D velocity model + GMT earth_relief topography /
# bathymetry + sea water + air (getParamsAndTopo, as the Kirishima box).
#
# From a terminal, in the repository root:
#   julia --project=. -t auto dev/japanTrench_nov2026/run.jl step1            # quick (5 km)
#   julia --project=. -t auto dev/japanTrench_nov2026/run.jl step1 standard   # 2.5 km
#   julia --project=. -t auto dev/japanTrench_nov2026/run.jl step1 full --redo=model
#   julia --project=. -t auto dev/japanTrench_nov2026/run.jl step2 coarse     # plane wave, OPT3 3-D (8 km, fast)
#   julia --project=. -t auto dev/japanTrench_nov2026/run.jl step3 coarse     # fences + record sections
#   julia --project=. -t auto dev/japanTrench_nov2026/run.jl step4 coarse     # movie (~20 s)
#   julia --project=. -t auto dev/japanTrench_nov2026/run.jl all quick
# From a REPL (fastest):
#   include("dev/japanTrench_nov2026/run.jl")      # default: step1 quick
#   japanTrench("step1", "standard")               # afterwards
# Outputs: dev/japanTrench_nov2026/outputs/step1_<case>/
# Cached stages (--redo=...): model, opt  ->  data/dev/japanTrench_nov2026/{model,opt_<case>}_<hash>.jld2

include(joinpath(@__DIR__, "..", "common", "devBootstrap.jl"))
JAPANTRENCH_STAGES = joinpath(@__DIR__, "stages")
japanTrench(args...) = devRun(joinpath(@__DIR__, "run.jl"), args...)

# zsh passes "# comment" text as arguments: drop everything from the first '#'.
japanTrenchArgs = let a = collect(String, @isdefined(DEV_ARGS) ? DEV_ARGS : ARGS)
    k = findfirst(x -> startswith(x, "#"), a)
    isnothing(k) ? a : a[1:k-1]
end
DEV = parseDevArgs(japanTrenchArgs;
    steps=["step1", "step2", "step3", "step4"], defaultStep="step1",
    cases=["coarse", "quick", "standard", "full"], defaultCases=["quick"])
length(DEV.cases) == 1 || error("give one grid case at a time (coarse | quick | standard | full)")
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
println("Japan Trench dev run: ", DEV)

include(joinpath(JAPANTRENCH_STAGES, "env.jl"))
include(joinpath(@__DIR__, "codeVersions.jl"))
japanTrenchCase = only(DEV.cases)
include(joinpath(JAPANTRENCH_STAGES, "config.jl"))

for stage in ("step1_box.jl", "step1_model.jl")
    println("---- ", stage)
    @time stage include(joinpath(JAPANTRENCH_STAGES, stage))
end
japanTrenchNeedOPT && include(joinpath(JAPANTRENCH_STAGES, "ocean_load.jl"))
japanTrenchOptLoaded = false        # per run: never reuse another case's optResult from the REPL
for step in DEV.steps
    DEV_STEP[] = step * "_" * japanTrenchCase
    println("\n================ Japan Trench ", step, " -- ", japanTrenchCase, " ================")
    if step == "step1"
        include(joinpath(JAPANTRENCH_STAGES, "step1_figures.jl"))
    elseif step == "step2"
        @time "step2" include(joinpath(JAPANTRENCH_STAGES, "step2_opt.jl"))
        global japanTrenchOptLoaded = true
    elseif step == "step3"
        japanTrenchOptLoaded || include(joinpath(JAPANTRENCH_STAGES, "step2_opt.jl"))   # cache hit
        global japanTrenchOptLoaded = true
        include(joinpath(JAPANTRENCH_STAGES, "step3_common.jl"))
        include(joinpath(JAPANTRENCH_STAGES, "step3_figures.jl"))
    elseif step == "step4"
        japanTrenchOptLoaded || include(joinpath(JAPANTRENCH_STAGES, "step2_opt.jl"))   # cache hit
        global japanTrenchOptLoaded = true
        include(joinpath(JAPANTRENCH_STAGES, "step3_common.jl"))
        include(joinpath(JAPANTRENCH_STAGES, "step4_movie.jl"))
    end
end
