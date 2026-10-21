# dev/spherical1Delastic/run.jl -- 1-D spherical TI SH (DSM-style), script version
#
# Same code as notebooks/SphericalTI_SH_oct2026.ipynb, split into stages/*.jl;
# the expensive parts are cached with myProduceOrLoad (explicitly, in each
# stage file). Steps:
#   step1    symbolic derivation and checks           (stages/step1_symbolic.jl)
#   step2    recipes, convergence, zone comparisons   (stages/step2_operators.jl)
#   step3a   shell eigenfrequencies                   (stages/step3a_eigen.jl)
#   step3a1  3a′: clipped closures                    (stages/step3a1_clipped.jl)
#   step3a2  3a″: fallback vs moving ν                (stages/step3a2_movingnu.jl)
#   step3b   PREM SH seismograms                      (stages/step3b_seismograms.jl)
#
# From a terminal, in the repository root:
#   julia --project=. dev/spherical1Delastic/run.jl step3b
#   julia --project=. dev/spherical1Delastic/run.jl step3b --redo=spectra
#   julia --project=. dev/spherical1Delastic/run.jl all
# From a REPL (fastest):
#   include("dev/spherical1Delastic/run.jl")      # first time (runs the default: step3b)
#   spherical("step3a2", "--redo=eigen")          # afterwards
# Figures: dev/spherical1Delastic/outputs/<step>/
# --redo= stage names: cartesian, convergence, exactness, zones, eigen, systems, spectra, all

include(joinpath(@__DIR__, "..", "common", "devBootstrap.jl"))
SPH_STAGES = joinpath(@__DIR__, "stages")
spherical(args...) = devRun(joinpath(@__DIR__, "run.jl"), args...)

SPH_STEP_FILES = Dict(
    "step1"   => "step1_symbolic.jl",
    "step2"   => "step2_operators.jl",
    "step3a"  => "step3a_eigen.jl",
    "step3a1" => "step3a1_clipped.jl",
    "step3a2" => "step3a2_movingnu.jl",
    "step3b"  => "step3b_seismograms.jl",
)
DEV = parseDevArgs(@isdefined(DEV_ARGS) ? DEV_ARGS : ARGS;
    steps=["step1", "step2", "step3a", "step3a1", "step3a2", "step3b"], defaultStep="step3b")
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
println("spherical1Delastic dev run: ", DEV)

include(joinpath(SPH_STAGES, "env.jl"))
include(joinpath(@__DIR__, "codeVersions.jl"))
include(joinpath(SPH_STAGES, "configs.jl"))
include(joinpath(SPH_STAGES, "plots.jl"))

for step in DEV.steps
    DEV_STEP[] = step
    println("\n================ ", step, " ================")
    @time step include(joinpath(SPH_STAGES, SPH_STEP_FILES[step]))
end
