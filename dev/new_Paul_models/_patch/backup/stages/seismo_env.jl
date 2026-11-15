# dev/new_Paul_models/stages/seismo_env.jl -- extra modules for synthetics/phases
seismoLoadEnv = !@isdefined(PAUL_SEISMO_ENV)
seismoLoadEnv && !@isdefined(planet1D) && include(joinpath(flexopt_root, "src", "planet1D.jl"))
seismoLoadEnv && planet1D.configure_input!()
seismoLoadEnv && !@isdefined(elasticWave2D) && include(joinpath(flexopt_root, "src", "elasticWave2D.jl"))
using .planet1D, .elasticWave2D
using Interpolations, JLD2, LinearAlgebra, DelimitedFiles
# ray theory (stdlib only, plain include into Main): PREM travel times and paths
seismoLoadEnv && include(joinpath(flexopt_root, "dev", "raytheory", "lib", "rayTheory.jl"))
PAUL_SEISMO_ENV = true
