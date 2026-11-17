# dev/japanTrench_nov2026/stages/env_opt.jl -- step2-4 environment (as Greenland's env.jl)
using Metal
!@isdefined(backend) && include(joinpath(flexopt_root, "src", "batchFiles", "batchGPU.jl"))
!@isdefined(flexOPT) && include(joinpath(flexopt_root, "src", "flexOPT.jl"))
using .flexOPT
!@isdefined(GreenlandIceWater) && include(joinpath(flexopt_root, "helpers", "greenlandIceWater.jl"))
using .GreenlandIceWater
using SparseArrays, StaticArrays
JAPANTRENCH_OPT_ENV = true
@show backend
