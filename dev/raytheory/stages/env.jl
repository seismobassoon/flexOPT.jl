# dev/raytheory/stages/env.jl
# The ray-theory core (lib/) uses only the Julia standard library: src/ is not
# loaded at all. CairoMakie is loaded only by the plotting steps (plots.jl).
# lib/*.jl and checks.jl are re-read with --reload.
rtLoad = !@isdefined(RAYTHEORY_ENV) || DEV.reload
@show VERSION Threads.nthreads() Base.active_project()
using LinearAlgebra, Printf
RT_DIR = joinpath(flexopt_root, "dev", "raytheory")
rtLoad && include(joinpath(RT_DIR, "lib", "rayTheory.jl"))
rtLoad && include(joinpath(RT_DIR, "stages", "checks.jl"))

# 1-D models from dataInput/model (DSM .psv.model format), by case name
RT_MODELS = Dict(
    "PREM"   => joinpath("earth", "PREM.psv.model"),
    "ak135"  => joinpath("earth", "ak135.psv.model"),
    "iasp"   => joinpath("earth", "iasp.psv.model"),
    "DWThot" => joinpath("mars", "DWThot.psv.model"),
    "TAYAKr" => joinpath("mars", "TAYAKr.psv.model"),
    "vprem"  => joinpath("moon", "vprem.psv.model"),
)
rtModel(name; isotropic=:voigt) = readDSMModel(joinpath(flexopt_root, "dataInput", "model", RT_MODELS[name]); isotropic=isotropic, name=name)
RT_TAUP_CSV = joinpath(RT_DIR, "reference", "taup_premDSM.csv")
RAYTHEORY_ENV = true
