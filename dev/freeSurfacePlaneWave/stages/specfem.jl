# dev/freeSurfacePlaneWave/stages/specfem.jl
#
# Genuine SPECFEM2D plane-wave reference: initialfield=.true. with
# source_type=1 ("normal", P wave including the free-surface reflection) or
# source_type=4 ("nosurface", incident P wave only -- no reflection/
# conversion at all) -- SPECFEM2D's own analytic plane-wave seeding, NOT a
# point source or a grid of sources. Same no-source-calibration spirit as
# the FD3/OPT3/OPT5 side of this benchmark (config.jl's planeWaveUz IC):
# the wave enters SPECFEM2D's domain as an exact analytic field too, built
# and verified by SPECFEM2D itself rather than by our own closed form.
# See stages/specfem_run.jl for the actual case construction; this file is
# only the myProduceOrLoad cache wrapper (same pattern as opt.jl/fd.jl).
#
# Toggle PLANEWAVE_RUN_SPECFEM2D[] = false to reuse a previous run's
# OUTPUT_FILES without re-meshing/re-solving (e.g. while only tweaking the
# trace post-processing or compare.jl's plotting).
if !@isdefined(PLANEWAVE_RUN_SPECFEM2D)
    const PLANEWAVE_RUN_SPECFEM2D = Ref(true)
end

specfemCacheParams = Dict{String,Any}(
    "stage" => "specfem", "codeVersion" => CODE_VERSION["specfem"],
    "grid" => (x0=first(x), z0=first(z), nx=nx, dx=dx),
    "duration" => duration, "dt" => dtOPT,
    "planewave" => (z0=z0, frequency=sourceFrequency, delay=sourceDelay,
        amplitude=planewaveAmplitude),
    "vp0" => vp0, "vs0" => vs0, "rho0" => rho0,
    "applyFreeSurface" => applyFreeSurface,
    "observerDepths" => collect(observerDepths),
    "runSPECFEM2D" => PLANEWAVE_RUN_SPECFEM2D[],
)
specfemCachePrefix = "specfem_planewave"
devRedo!("specfem", specfemCacheParams, PLANEWAVE_CACHE_DIR, specfemCachePrefix)
specfemCacheFile = devCacheFile(specfemCacheParams, PLANEWAVE_CACHE_DIR, specfemCachePrefix)
println("SPECFEM2D cache: ", specfemCacheFile, isfile(specfemCacheFile) ? "  (hit)" : "  (miss -> computing)")
specfemResult = myProduceOrLoad(specfemCacheParams, PLANEWAVE_CACHE_DIR, specfemCachePrefix) do config
    include(joinpath(PLANEWAVE_STAGES, "specfem_run.jl"))
end
specfemTracesZ = specfemResult["specfemTracesZ"]
specfemCase = specfemResult["specfemCase"]
specfemObserverDepths = specfemResult["specfemObserverDepths"]
specfemCFL = specfemResult["specfemCFL"]
@show specfemCase.case_directory specfemCase.time_axis_shift specfemCFL length(specfemTracesZ)
