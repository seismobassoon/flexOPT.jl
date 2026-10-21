# dev/kirishima_oct2026/run.jl  -- Kirishima elastic 2-D benchmark, script version
#
# Same code as notebooks/Elastic2DtopoBenchmark_oct2026.ipynb (Steps 1-2),
# split into stages/*.jl; the expensive stages are cached with
# myProduceOrLoad (see stages/step2_opt.jl, step2_fd.jl, step2_specfem.jl).
#
# From a terminal, in the repository root:
#   julia --project=. dev/kirishima_oct2026/run.jl step2              # both cases
#   julia --project=. dev/kirishima_oct2026/run.jl step2 heter        # one case
#   julia --project=. dev/kirishima_oct2026/run.jl step2 heter --redo=opt
#   julia --project=. dev/kirishima_oct2026/run.jl step1
# From a REPL (fastest: everything stays compiled and loaded):
#   include("dev/kirishima_oct2026/run.jl")      # first time (runs the default: step2, both)
#   kirishima("step2", "heter", "--redo=opt")    # afterwards
# Figures and metrics: dev/kirishima_oct2026/outputs/<step>[_<case>]/

include(joinpath(@__DIR__, "..", "common", "devBootstrap.jl"))
KIRISHIMA_STAGES = joinpath(@__DIR__, "stages")
kirishima(args...) = devRun(joinpath(@__DIR__, "run.jl"), args...)

DEV = parseDevArgs(@isdefined(DEV_ARGS) ? DEV_ARGS : ARGS;
    steps=["step1", "step2"], defaultStep="step2",
    cases=["homog", "heter"])
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
println("Kirishima dev run: ", DEV)

include(joinpath(KIRISHIMA_STAGES, "env_step1.jl"))
include(joinpath(@__DIR__, "codeVersions.jl"))
include(joinpath(KIRISHIMA_STAGES, "shared_site.jl"))

STEP2_CASE_STAGES = (
    "step2_model.jl",     # homog / heter fields, safety checks, modelHash
    "step2_opt.jl",       # cached: recipes + operators + convFD3/OPT3/(OPT4)/OPT5 propagation
    "step2_fd.jl",        # cached: explicit FD3
    "step2_traces.jl",    # receiver traces, pre-SPECFEM figure
    "step2_specfem.jl",   # cached: SPECFEM2D reference
    "step2_compare.jl",   # metrics vs SPECFEM2D, figure + metrics file
)

for step in DEV.steps
    DEV_STEP[] = step
    if step == "step1"
        @time "step1" include(joinpath(KIRISHIMA_STAGES, "step1_box.jl"))
    elseif step == "step2"
        include(joinpath(KIRISHIMA_STAGES, "step2_env.jl"))
        include(joinpath(KIRISHIMA_STAGES, "step2_box.jl"))
        include(joinpath(KIRISHIMA_STAGES, "step2_config.jl"))
        for case in DEV.cases
            global step2Case = Symbol(case)
            DEV_STEP[] = "step2_" * case
            println("\n================ Step 2 -- ", case, " ================")
            for stage in STEP2_CASE_STAGES
                println("---- ", stage)
                @time stage include(joinpath(KIRISHIMA_STAGES, stage))
            end
        end
    end
end
