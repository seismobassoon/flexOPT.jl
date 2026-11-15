# dev/new_Paul_models/stages/env.jl
# Load src/ modules once per Julia session (after editing src/: restart Julia).
paulLoadEnv = !@isdefined(PAUL_ENV)
paulLoadEnv && !@isdefined(commonBatchs) && include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
using .commonBatchs
using CairoMakie, Statistics, Printf
CairoMakie.activate!(type = "png")
PAUL_ENV = true
@show VERSION Threads.nthreads()
