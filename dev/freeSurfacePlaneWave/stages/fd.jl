# dev/freeSurfacePlaneWave/stages/fd.jl
fdCacheParams = Dict{String,Any}(
    "stage" => "fd", "codeVersion" => CODE_VERSION["fd"],
    "grid" => (x0=first(x), z0=first(z), nx=nx, nz=nz, dx=dx),
    "duration" => duration, "outputSampling" => outputSampling,
    "source" => (position=sourcePosition, frequency=sourceFrequency, delay=sourceDelay,
        force=sourceForce, nColumns=nx),
    "boundaries" => (cerjanCells=cerjanCells, damping=cerjan.damping),
    "vp0" => vp0, "vs0" => vs0, "rho0" => rho0,
)
fdCachePrefix = "fd_planewave"
devRedo!("fd", fdCacheParams, PLANEWAVE_CACHE_DIR, fdCachePrefix)
fdCacheFile = devCacheFile(fdCacheParams, PLANEWAVE_CACHE_DIR, fdCachePrefix)
println("FD3 cache: ", fdCacheFile, isfile(fdCacheFile) ? "  (hit)" : "  (miss -> computing)")
fdResult = myProduceOrLoad(fdCacheParams, PLANEWAVE_CACHE_DIR, fdCachePrefix) do config
    include(joinpath(PLANEWAVE_STAGES, "fd_propagate.jl"))
end
uxFD, uzFD, fdTimes = fdResult["uxFD"], fdResult["uzFD"], fdResult["fdTimes"]
fdX, fdZ, fdPadding, fdDt = fdResult["fdX"], fdResult["fdZ"], fdResult["fdPadding"], fdResult["fdDt"]
@show fdDt size(uzFD) maximum(abs, uzFD)
