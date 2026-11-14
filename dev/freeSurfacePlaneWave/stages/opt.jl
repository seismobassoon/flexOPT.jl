# dev/freeSurfacePlaneWave/stages/opt.jl
optCacheParams = Dict{String,Any}(
    "stage" => "opt", "codeVersion" => CODE_VERSION["opt"],
    "grid" => (x0=first(x), z0=first(z), nx=nx, nz=nz, dx=dx),
    "duration" => duration, "outputSampling" => outputSampling,
    "source" => (position=sourcePosition, frequency=sourceFrequency, delay=sourceDelay,
        force=sourceForce, nColumns=length(sourceNodesX)),
    "boundaries" => (cerjanCells=cerjanCells, damping=cerjan.damping),
    "vp0" => vp0, "vs0" => vs0, "rho0" => rho0,
)
optCachePrefix = "opt_planewave"
devRedo!("opt", optCacheParams, PLANEWAVE_CACHE_DIR, optCachePrefix)
optCacheFile = devCacheFile(optCacheParams, PLANEWAVE_CACHE_DIR, optCachePrefix)
println("OPT family cache: ", optCacheFile, isfile(optCacheFile) ? "  (hit)" : "  (miss -> computing)")
optResult = myProduceOrLoad(optCacheParams, PLANEWAVE_CACHE_DIR, optCachePrefix) do config
    include(joinpath(PLANEWAVE_STAGES, "opt_build.jl"))
    include(joinpath(PLANEWAVE_STAGES, "opt_propagate.jl"))
end
xOPT, zOPT, dxOPT, dtOPT = optResult["xOPT"], optResult["zOPT"], optResult["dxOPT"], optResult["dtOPT"]
uxOPT, uzOPT, optTimes = optResult["uxOPT"], optResult["uzOPT"], optResult["optTimes"]
uxOPT5, uzOPT5, opt5Times = optResult["uxOPT5"], optResult["uzOPT5"], optResult["opt5Times"]
