# dev/japanTrench_nov2026/stages/env.jl
# Load src/ modules once per Julia session (after editing src/: restart Julia).
# step1 needs commonBatchs + planet1D + GeoPoints; step2-4 also flexOPT, the
# GPU/KernelAbstractions plumbing and helpers/greenlandIceWater.jl (3-D free
# surface closure, memory estimate). Every line is its own top-level statement.
japanTrenchLoadEnv = !@isdefined(JAPANTRENCH_ENV)
japanTrenchNeedOPT = any(!=("step1"), DEV.steps)
japanTrenchLoadEnv && !@isdefined(commonBatchs) && include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
japanTrenchLoadEnv && !@isdefined(planet1D) && include(joinpath(flexopt_root, "src", "planet1D.jl"))
japanTrenchLoadEnv && planet1D.configure_input!()
japanTrenchLoadEnv && !@isdefined(GeoPoints) && include(joinpath(flexopt_root, "src", "GeoPoints.jl"))
using .commonBatchs, .planet1D, .GeoPoints
using CairoMakie, Statistics, LinearAlgebra, Printf, Interpolations
CairoMakie.activate!(type="png")
JAPANTRENCH_ENV = true
@show VERSION Threads.nthreads()

japanTrenchNeedOPT && !@isdefined(JAPANTRENCH_OPT_ENV) && include(joinpath(@__DIR__, "env_opt.jl"))
