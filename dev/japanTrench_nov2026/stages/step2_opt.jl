# dev/japanTrench_nov2026/stages/step2_opt.jl
#
# OPT3 3-D elastic propagation of a plane wave, cached with myProduceOrLoad.
# On a cache hit nothing is assembled (surface wavefield + sections loaded).
# On a miss:
#   step2_build.jl -- OPT3 recipe (same cache as Greenland, data/semiSymbolic),
#                     operators, Cerjan (sides + bottom), 3-D free surface;
#   step2_march.jl -- plane-wave sheet source + predictor-corrector marching.

dtOPT = courant * dx / (sqrt(3) * vpMaximum * 1e3)
optSteps = ceil(Int, duration / dtOPT)
traceStride = max(1, round(Int, traceSampling / dtOPT))
snapStride = traceStride * max(1, round(Int, snapshotEvery / traceSampling))
nWholePoints = (Nx + 2cerjanCells) * (Ny + 2cerjanCells) * (Nz + cerjanCells)
@printf("dt = %.4f s, %d steps for %.0f s; whole grid %.2f M points (%d Cerjan cells), operators ≈ %.1f GB\n",
    dtOPT, optSteps, duration, nWholePoints / 1e6, cerjanCells, estimate_operator_memory_gb(nWholePoints))

optCacheParams = Dict{String,Any}(
    "stage" => "opt",
    "codeVersion" => CODE_VERSION["opt"],
    "modelHash" => modelHash,
    "case" => japanTrenchCase,
    "vsMinSolid" => VS_MIN_SOLID,
    "oceanLoad" => string(OCEAN_LOAD),
    "dt" => dtOPT, "duration" => duration,
    "cerjan" => (cells=cerjanCells, damping=cerjanDamping),
    "planeWave" => merge(PLANE_WAVE, (type=string(PLANE_WAVE.type),)),
    "recipe" => (supplementaryOrder=optSupplementaryOrder, gpu=recipeOnGPU),
    "marching" => (scheme=string(optTimeMarchingScheme), predictor=string(optPredictorMode),
        iterations=optCorrectorIterations, tolerance=optCorrectorTolerance,
        relaxation=optCorrectorRelaxation),
    "output" => (traceStride=traceStride, snapStride=snapStride),
)
optPrefix = "opt_$(japanTrenchCase)"
devRedo!("opt", optCacheParams, JAPANTRENCH_CACHE_DIR, optPrefix)
optCacheFile = devCacheFile(optCacheParams, JAPANTRENCH_CACHE_DIR, optPrefix)
println("OPT cache: ", optCacheFile, isfile(optCacheFile) ? "  (hit)" : "  (miss -> computing)")

optResult = myProduceOrLoad(optCacheParams, JAPANTRENCH_CACHE_DIR, optPrefix) do cfg
    include(joinpath(JAPANTRENCH_STAGES, "step2_build.jl"))
    include(joinpath(JAPANTRENCH_STAGES, "step2_march.jl"))   # last expression = Dict
end

surfaceTimes = optResult["surfaceTimes"]
surfaceU = optResult["surfaceU"]          # (x|y|z of the box, i, j, time), m, on the top rock node
snapTimes = optResult["snapTimes"]
secEW = optResult["secEW"]                # (component, i, k, snapshot) along y = 0
secNS = optResult["secNS"]                # (component, j, k, snapshot) along x = 0
sheet = optResult["sheet"]
pc = optResult["pcDiagnostics"]
@printf("plane wave %s: sheet at z = %.1f km, V_ref = %.2f km/s, p = %.4f s/km (apparent velocity %.2f km/s)\n",
    PLANE_WAVE.type, sheet.z0 / 1e3, sheet.vRef / 1e3, sheet.p * 1e3, 1e-3 / sheet.p)
@printf("P-C: iterations per step %s, accepted relative residual in [%.1e, %.1e]\n",
    string(extrema(pc.iterations)), extrema(pc.residuals)...)
@show optResult["timings"]
optResult["stoppedEarly"] && @warn "the run stopped early (blow-up guard); output is partial"
