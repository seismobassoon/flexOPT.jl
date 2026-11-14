# dev/freeSurfacePlaneWave/stages/env.jl
# Minimal environment: no GeoPoints/planet1D/NIED -- this is a synthetic
# homogeneous box, not tied to any real location.
planewaveLoadEnv = !@isdefined(PLANEWAVE_ENV)
planewaveLoadEnv && !@isdefined(commonBatchs) && include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
using .commonBatchs
using Metal
planewaveLoadEnv && !@isdefined(backend) && include(joinpath(flexopt_root, "src", "batchFiles", "batchGPU.jl"))
backend isa KernelAbstractions.CPU && @warn("No GPU backend is available; OPT recipes will be built on CPU")
planewaveLoadEnv && !@isdefined(flexOPT) && include(joinpath(flexopt_root, "src", "flexOPT.jl"))
planewaveLoadEnv && !@isdefined(elasticWave2D) && include(joinpath(flexopt_root, "src", "elasticWave2D.jl"))
using .flexOPT, .elasticWave2D
using CairoMakie, Statistics, LinearAlgebra, SparseArrays, JLD2
CairoMakie.activate!(type="png")
PLANEWAVE_ENV = true
@show VERSION Threads.nthreads() backend
