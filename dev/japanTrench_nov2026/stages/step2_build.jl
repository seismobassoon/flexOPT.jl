# dev/japanTrench_nov2026/stages/step2_build.jl
# Runs ONLY on a cache miss of step2_opt.jl. Same construction as
# dev/greenland_nov2026/stages/step2_build.jl:
#   1. OPT3 recipe of famousEquations("3DsismoTimeIsoHeteroForce") in scaled
#      coordinates (Δ = 1): identical parameters -> the same cache file as Greenland;
#   2. operators on this model, Cerjan on the four sides and the bottom;
#   3. apply_free_surface_3d (helpers): air pinned, Dietrich mirror on the
#      rock nodes next to air (land surface and sea floor).

# ---- 1. recipe
recipeBackend = recipeOnGPU ? backend : KernelAbstractions.CPU()
optParameters = Dict{String,Any}(
    "famousEquationType" => "3DsismoTimeIsoHeteroForce",
    "Δ" => (1.0, 1.0, 1.0, 1.0),
    "orderBtime" => 1, "orderBspace" => 1,
    "pointsInSpace" => 3, "pointsInTime" => 3,
    "supplementaryOrder" => optSupplementaryOrder,
    "taylorInverseMode" => :hierarchical_constrained,
    "fieldItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0,
        offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "materItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0,
        offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "recipe_backend" => recipeBackend,
)
recipeCacheParameters = Dict{String,Any}(k => v for (k, v) in optParameters if k != "recipe_backend")
recipeCacheParameters["recipe_cache_version"] = 3
recipeCacheParameters["recipe_backend_name"] = string(nameof(typeof(recipeBackend)))
recipeFile = devCacheFile(recipeCacheParameters, "semiSymbolic", "elastic3DForce_OPT3")
println("OPT3 3-D recipe: ", recipeFile, isfile(recipeFile) ? "  (hit)" : "  (miss -> makeOPTsemiSymbolic, slow)")
timedRecipe = @timed myProduceOrLoad(recipeCacheParameters, "semiSymbolic", "elastic3DForce_OPT3") do cfg
    runtime = Dict{String,Any}(cfg)
    for key in ("hash_id", "recipe_cache_version", "recipe_backend_name")
        pop!(runtime, key, nothing)
    end
    runtime["recipe_backend"] = recipeBackend
    makeOPTsemiSymbolic(runtime)
end
optRecipe = timedRecipe.value

# ---- 2. operators: SI, scaled as Kirishima/Greenland. λ, μ from the rock ρ;
# ρ' from the loaded ρ (ocean mass on the top rock cells).
ρSI = ρModel .* 1e3
vpSI = VpModel .* 1e3
vsSI = VsModel .* 1e3
rho0 = median(ρSI[materialFS])
scale2 = (dtOPT / dx)^2
rhoOPT = max.(ρLoaded .* 1e3 ./ rho0, 1e-3)
muOPT = (ρSI .* vsSI .^ 2 ./ rho0) .* scale2
lambdaOPT = (ρSI .* (vpSI .^ 2 .- 2 .* vsSI .^ 2) ./ rho0) .* scale2
muOPT[.!materialFS] .= 0.0                 # air (and ocean) rows are pinned anyway
lambdaOPT[.!materialFS] .= 0.0
@assert all(isfinite, lambdaOPT) && all(isfinite, muOPT)
@assert minimum(lambdaOPT[materialFS] .+ (2 / 3) .* muOPT[materialFS]) > 0 "negative bulk modulus"

cerjanSpec = CerjanBoundarySpec((cerjanCells, cerjanCells, cerjanCells),
    (cerjanCells, cerjanCells, 0); damping=cerjanDamping)
modelsOPT = [rhoOPT, lambdaOPT, muOPT]
pointsOPT = getModelPoints(modelsOPT[1], 3,
    optRecipe["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
familyOPT = (models=modelsOPT, modelPoints=pointsOPT, Δ=(1.0, 1.0, 1.0, 1.0),
    modelName="japanTrench_$(japanTrenchCase)_$(modelHash)_$(OCEAN_LOAD)")
timedAssembly = @timed numericalOperatorConstruction(Dict{String,Any}(
    "optRec" => optRecipe, "modelFam" => familyOPT,
    "absorbingBoundaries" => nothing,
    "maskedRegionInSpace" => nothing,
    "boundaryConditions" => BoundaryConditionSet(free_surface=nothing, cerjan=cerjanSpec,
        material_mask=nothing),
    "representation" => "matrixfree",
))["numOperators"]
timedPreparation = @timed prepareLinearSystem(timedAssembly.value)
preparedVolume = timedPreparation.value
optPadding = cerjan_padding(cerjanSpec)                    # [lower; upper] × (x, y, z)
@assert preparedVolume.spaceShape == (Nx, Ny, Nz) .+ Tuple(optPadding[1, :]) .+ Tuple(optPadding[2, :])
@assert preparedVolume.NField == 3 && preparedVolume.NForceField == 3
@assert preparedVolume.timePointsUsedForOneStep == 3

# ---- 3. free surface (land surface + sea floor)
timedClosure = @timed apply_free_surface_3d(preparedVolume, materialFS;
    prepared_with_matrices=flexOPT._prepared_with_matrices)
preparedOPT = timedClosure.value
operatorTimings = (recipe=timedRecipe.time, assembly=timedAssembly.time,
    preparation=timedPreparation.time, freeSurface=timedClosure.time)
preparedVolume = nothing; timedAssembly = nothing; timedPreparation = nothing; timedClosure = nothing
GC.gc()
@show operatorTimings
@printf("OPT3 3-D: whole grid %s, %d unknowns, nnz(A) = %.1f M, nnz(L) = %.1f M, nnz(R) = %.1f M\n",
    string(preparedOPT.spaceShape), size(preparedOPT.A_unknown, 1),
    nnz(preparedOPT.A_unknown) / 1e6, nnz(preparedOPT.L_known) / 1e6, nnz(preparedOPT.R_force) / 1e6)
