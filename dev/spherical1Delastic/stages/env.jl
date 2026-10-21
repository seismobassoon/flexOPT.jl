# dev/spherical1Delastic/stages/env.jl
# Ported from notebooks/SphericalTI_SH_oct2026.ipynb, cell sh-step1-01 (+ the
# helper includes of cells sh-step2-01, sh-step3-01, sh-step3b-01).
# src/ modules: loaded once per Julia session (after editing src/, restart
# Julia). helpers/sphericalTI_SH*.jl: re-read with --reload. One top-level
# statement per line on purpose (bindings created by include must be
# visible to the next statement).
sphLoadEnv = !@isdefined(SPHERICAL_ENV)
sphLoadHelpers = sphLoadEnv || DEV.reload
@show VERSION Threads.nthreads() Base.active_project()
# Metal must be loaded before batchGPU.jl selects the backend.
using Metal
Metal.functional() || error("Metal.jl is loaded, but cannot access the Apple GPU")
sphLoadEnv && include(joinpath(flexopt_root, "src", "batchFiles", "batchGPU.jl"))
sphLoadEnv && include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
sphLoadEnv && include(joinpath(flexopt_root, "src", "planet1D.jl"))
sphLoadEnv && planet1D.configure_input!()
sphLoadEnv && include(joinpath(flexopt_root, "src", "flexOPT.jl"))
using .commonBatchs, .planet1D, .flexOPT
using Symbolics, LinearAlgebra, Random, SparseArrays, Printf
sphLoadHelpers && include(joinpath(flexopt_root, "helpers", "sphericalTI_SH.jl"))
sphLoadHelpers && include(joinpath(flexopt_root, "helpers", "sphericalTI_SH_step2.jl"))
sphLoadHelpers && include(joinpath(flexopt_root, "helpers", "sphericalTI_SH_step3.jl"))
sphLoadHelpers && include(joinpath(flexopt_root, "helpers", "sphericalTI_SH_step3b.jl"))
using CairoMakie
CairoMakie.activate!(type="png")
SPHERICAL_ENV = true
