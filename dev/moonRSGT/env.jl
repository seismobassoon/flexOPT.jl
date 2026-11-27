# dev/moonRSGT/env.jl -- solver environment WITHOUT Metal (Linux clusters: stella).
# Same loading order as dev/spherical1Delastic/stages/env.jl and
# dev/spherical1DelasticKernel/stages/env.jl, except `using Metal`:
# batchGPU.jl then selects the CPU backend (the 1-D spherical code is CPU
# anyway; recipes are always built on the CPU).
# One top-level statement per line on purpose (bindings created by include
# must be visible to the next statement).
moonLoadEnv = !@isdefined(MOON_ENV)
@show VERSION Threads.nthreads() Base.active_project()
moonLoadEnv && !@isdefined(SPHERICAL_ENV) && include(joinpath(flexopt_root, "src", "batchFiles", "batchGPU.jl"))
moonLoadEnv && !@isdefined(SPHERICAL_ENV) && include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
moonLoadEnv && !@isdefined(SPHERICAL_ENV) && include(joinpath(flexopt_root, "src", "planet1D.jl"))
moonLoadEnv && !@isdefined(SPHERICAL_ENV) && planet1D.configure_input!()
moonLoadEnv && !@isdefined(SPHERICAL_ENV) && include(joinpath(flexopt_root, "src", "flexOPT.jl"))
using .commonBatchs, .planet1D, .flexOPT
using Symbolics, LinearAlgebra, Random, SparseArrays, Printf
MOON_SPH_LIB = joinpath(flexopt_root, "dev", "spherical1Delastic", "lib")
MOON_KER_LIB = joinpath(flexopt_root, "dev", "spherical1DelasticKernel", "lib")
moonLoadEnv && include(joinpath(MOON_SPH_LIB, "sphericalTI_SH.jl"))
moonLoadEnv && include(joinpath(MOON_SPH_LIB, "sphericalTI_SH_step2.jl"))
moonLoadEnv && include(joinpath(MOON_SPH_LIB, "sphericalTI_SH_step3.jl"))
moonLoadEnv && include(joinpath(MOON_SPH_LIB, "sphericalTI_SH_step3b.jl"))
moonLoadEnv && include(joinpath(MOON_SPH_LIB, "sphericalTI_PSV.jl"))
moonLoadEnv && include(joinpath(MOON_SPH_LIB, "elastic1D.jl"))
moonLoadEnv && include(joinpath(MOON_SPH_LIB, "sources1D.jl"))
moonLoadEnv && include(joinpath(MOON_SPH_LIB, "nearBoundary1D.jl"))
moonLoadEnv && include(joinpath(flexopt_root, "dev", "spherical1Delastic", "stages", "configs.jl"))   # shConfigs
moonLoadEnv && include(joinpath(MOON_KER_LIB, "sgtHarmonics1D.jl"))
moonLoadEnv && include(joinpath(MOON_KER_LIB, "sgtDisk1D.jl"))
moonLoadEnv && include(joinpath(MOON_KER_LIB, "sgtGeometry1D.jl"))
moonLoadEnv && include(joinpath(MOON_KER_LIB, "sgtInterp1D.jl"))
MOON_LIB = joinpath(flexopt_root, "dev", "moonRSGT", "lib")
include(joinpath(MOON_LIB, "moonModel.jl"))
include(joinpath(MOON_LIB, "moonCatalogue.jl"))
include(joinpath(MOON_LIB, "moonSynthetics.jl"))
include(joinpath(MOON_LIB, "moonChecks.jl"))
include(joinpath(MOON_LIB, "moonExportH5.jl"))
moonRecipes = @isdefined(moonRecipes) ? moonRecipes : Dict{Any,Any}()
MOON_ENV = true
