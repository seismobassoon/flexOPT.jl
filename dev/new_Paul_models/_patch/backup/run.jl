# dev/new_Paul_models/run.jl -- figures of every field of Paul's StagYY models
#
# Data: <PAUL_MODELS_DIR>/<model>/Output/<prefix>_<field><snapshot>, read with
# src/batchFiles/batchStagYY.jl (see lib/paulStagYY.jl for the two additions:
# vector `vp` files and surface `cs`/`hf`/`cr` files).
#
# From a terminal, in the repository root:
#   julia --project=. -t auto dev/new_Paul_models/run.jl                    # figures, all models
#   julia --project=. -t auto dev/new_Paul_models/run.jl figures G1m20_sb21  # one model
#   julia --project=. -t auto dev/new_Paul_models/run.jl synthetics         # FD3: PREM, sb11, sb21 x N/S pole
#   julia --project=. -t auto dev/new_Paul_models/run.jl phases             # phase windows (after synthetics)
#   julia --project=. -t auto dev/new_Paul_models/run.jl synthetics --redo=fd3   # recompute the runs
# From a REPL (fastest):
#   include("dev/new_Paul_models/run.jl")
#   paul("figures", "G1m20_sb11")
# Another data folder: ENV["PAUL_MODELS_DIR"] = "/path/to/folder" (one sub-folder per model).
# Outputs: dev/new_Paul_models/outputs/figures_<model>/<snapshot>_<field>.png
#          and <snapshot>_00_overview.png (all 2-D fields on one sheet).

include(joinpath(@__DIR__, "..", "common", "devBootstrap.jl"))
PAUL_STAGES = joinpath(@__DIR__, "stages")
paul(args...) = devRun(joinpath(@__DIR__, "run.jl"), args...)

PAUL_MODELS_DIR = get(ENV, "PAUL_MODELS_DIR", "/Users/nobuaki/IPGP Dropbox/nobuaki fuji/PhysChemSelfConsistantModels_for_neutrino_and_seismic")
isdir(PAUL_MODELS_DIR) || error("PAUL_MODELS_DIR not found: $PAUL_MODELS_DIR")
paulModels = sort(filter(d -> isdir(joinpath(PAUL_MODELS_DIR, d, "Output")), readdir(PAUL_MODELS_DIR)))
isempty(paulModels) && error("no <model>/Output folder in $PAUL_MODELS_DIR")

paulArgs = let a = collect(String, @isdefined(DEV_ARGS) ? DEV_ARGS : ARGS)
    k = findfirst(x -> startswith(x, "#"), a)
    isnothing(k) ? a : a[1:k-1]
end
DEV = parseDevArgs(paulArgs; steps = ["figures", "synthetics", "phases"], defaultStep = "figures",
    cases = paulModels, defaultCases = paulModels)
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
println("Paul models dev run: ", DEV)

include(joinpath(PAUL_STAGES, "env.jl"))
include(joinpath(@__DIR__, "lib", "paulStagYY.jl"))
include(joinpath(@__DIR__, "codeVersions.jl"))
include(joinpath(PAUL_STAGES, "config.jl"))

for step in DEV.steps
    if step == "figures"
        for model in DEV.cases
            global paulModel = model
            DEV_STEP[] = "figures_" * model
            println("\n================ Paul models: ", model, " ================")
            @time "figures $model" include(joinpath(PAUL_STAGES, "figures.jl"))
        end
    else
        include(joinpath(PAUL_STAGES, "seismo_env.jl"))
        include(joinpath(@__DIR__, "lib", "paulSeismo.jl"))
        include(joinpath(@__DIR__, "lib", "paulSeismoFigures.jl"))
        println("\n================ Paul models: ", step, " ================")
        step == "synthetics" && @time "synthetics" include(joinpath(PAUL_STAGES, "synthetics.jl"))
        step == "phases" && @time "phases" include(joinpath(PAUL_STAGES, "phases.jl"))
    end
end
