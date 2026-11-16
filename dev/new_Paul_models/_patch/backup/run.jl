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
#   julia --project=. -t auto dev/new_Paul_models/run.jl features           # discriminating features sb11 vs sb21
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
DEV = parseDevArgs(paulArgs; steps = ["figures", "synthetics", "phases", "features"], defaultStep = "figures",
    cases = paulModels, defaultCases = paulModels)
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
println("Paul models dev run: ", DEV)

include(joinpath(PAUL_STAGES, "env.jl"))
include(joinpath(@__DIR__, "lib", "paulStagYY.jl"))
include(joinpath(@__DIR__, "codeVersions.jl"))
include(joinpath(PAUL_STAGES, "config.jl"))

# every figure also goes to the shared Dropbox folder (empty string: off)
PAUL_FIGURE_MIRROR = get(ENV, "PAUL_FIGURE_MIRROR",
    joinpath(PAUL_MODELS_DIR, "tell_me_if_you_dont_want_these_files_that_could_take_some_space"))
"devSave + copy to PAUL_FIGURE_MIRROR/<step>/"
function paulSave(fig, name; attempts = 3)
    path = nothing
    for k in 1:attempts
        try
            path = devSave(fig, name)
            break
        catch err
            # e.g. "Operation timed out" on close() when the folder is synced / busy
            @warn "saving $name failed (attempt $k/$attempts)" exception = err
            k < attempts && sleep(5k)
        end
    end
    path === nothing && (@warn "figure $name NOT saved"; return nothing)
    if !isempty(PAUL_FIGURE_MIRROR)
        try
            dir = joinpath(PAUL_FIGURE_MIRROR, DEV_STEP[])
            mkpath(dir)
            cp(path, joinpath(dir, basename(path)); force = true)
        catch err
            @warn "copy of $name to the mirror failed" exception = err
        end
    end
    return path
end

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
        include(joinpath(@__DIR__, "lib", "paulFeatures.jl"))
        println("\n================ Paul models: ", step, " ================")
        global PAUL_SYNTH_FIGURES = step == "synthetics"   # phases/features reuse the runs without redrawing
        step == "synthetics" && @time "synthetics" include(joinpath(PAUL_STAGES, "synthetics.jl"))
        step == "phases" && @time "phases" include(joinpath(PAUL_STAGES, "phases.jl"))
        step == "features" && @time "features" include(joinpath(PAUL_STAGES, "features.jl"))
    end
end
