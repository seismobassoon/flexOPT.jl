# dev/freeSurfacePlaneWave/stages/opt_build.jl
# Recipes + operators for OPT3 and OPT5, closed at the free surface with the
# SAME additive elasticTractionFree2D closure used (but never yet confirmed
# working) in dev/kirishima_oct2026/stages/step2_opt_build.jl -- the thing
# actually under test here. Runs only on a cache miss of opt.jl.
xOPTPhysical = x
zOPTPhysical = z
solidOPT = repeat(reshape(zOPTPhysical .<= 0.0, 1, :), length(xOPTPhysical), 1)
dxOPT = dx
dtOPT = 0.20 * dxOPT / (sqrt(2) * vp0)   # kept simple: no temporalRefinement here
vpFieldOPT = fill(vp0, size(solidOPT))
vsFieldOPT = fill(vs0, size(solidOPT))
rhoFieldOPT = fill(rho0, size(solidOPT))
recipeSpacing = (1.0, 1.0, 1.0)   # scaled-coordinates convention, as elsewhere
rhoOPT = rhoFieldOPT ./ rho0
muOPT = (rhoFieldOPT .* vsFieldOPT .^ 2 ./ rho0) .* (dtOPT / dxOPT)^2
lambdaOPT = (rhoFieldOPT .* (vpFieldOPT .^ 2 .- 2 .* vsFieldOPT .^ 2) ./ rho0) .* (dtOPT / dxOPT)^2
muOPT[.!solidOPT] .= 0.0
lambdaOPT[.!solidOPT] .= 0.0
optParameters = Dict{String,Any}(
    "famousEquationType" => "2DsismoTimeIsoHeteroSingleForce",
    "Δ" => recipeSpacing,
    "orderBtime" => 1, "orderBspace" => 1,
    "pointsInSpace" => 3, "pointsInTime" => 3,
    "supplementaryOrder" => 2,
    "taylorInverseMode" => :hierarchical_constrained,
    "fieldItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0,
        offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "materItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0,
        offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "recipe_backend" => backend,
)
recipeCacheDirectory = "semiSymbolic"
recipeCacheVersion = 3
function cachedOPTRecipe(parameters, prefix)
    cacheParameters = Dict{String,Any}(
        key => value for (key, value) in parameters if key != "recipe_backend")
    cacheParameters["recipe_cache_version"] = recipeCacheVersion
    function produceRecipe(config)
        runtimeParameters = Dict{String,Any}(config)
        pop!(runtimeParameters, "hash_id", nothing)
        pop!(runtimeParameters, "recipe_cache_version", nothing)
        runtimeParameters["recipe_backend"] = backend
        makeOPTsemiSymbolic(runtimeParameters)
    end
    myProduceOrLoad(produceRecipe, cacheParameters, recipeCacheDirectory, prefix)
end
optRecipe = cachedOPTRecipe(optParameters, "elastic2D_OPT3")
opt5Parameters = copy(optParameters)
opt5Parameters["orderBspace"] = 1
opt5Parameters["pointsInSpace"] = 5
opt5Parameters["fieldItpl"] = merge(optParameters["fieldItpl"], (offsetSpace=2.0,))
opt5Parameters["materItpl"] = merge(optParameters["materItpl"], (offsetSpace=2.0,))
opt5Recipe = cachedOPTRecipe(opt5Parameters, "elastic2D_OPT5")
coefficientGatePath = joinpath(flexopt_root, "data", "elastic_lhs_coefficient_gate.jld2")
isfile(coefficientGatePath) || error("Run ElasticLHSOperatorAudit.ipynb before OPT5 propagation")
coefficientGate = load(coefficientGatePath)
@assert all(result.passed for result in values(coefficientGate["coefficient_health"]))

optCerjan = CerjanBoundarySpec(
    (cerjan.lower[1] * optSpatialRefinement, cerjan.lower[2] * optSpatialRefinement),
    (cerjan.upper[1] * optSpatialRefinement, 0);
    damping=cerjan.damping / optSpatialRefinement^2,
)
bcOPT = boundary_geometry(solidOPT, (dxOPT, dxOPT); free_surface_mode=:pinned_void, cerjan=optCerjan)
modelsOPT = [rhoOPT, lambdaOPT, muOPT]
pointsOPT = getModelPoints(modelsOPT[1], 3,
    optRecipe["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
familyOPT = (models=modelsOPT, modelPoints=pointsOPT, Δ=recipeSpacing, modelName="planewave_OPT3")
numericalVolume = numericalOperatorConstruction(Dict{String,Any}(
    "optRec" => optRecipe, "modelFam" => familyOPT,
    "absorbingBoundaries" => nothing, "maskedRegionInSpace" => nothing,
    "boundaryConditions" => bcOPT, "representation" => "matrixfree",
))["numOperators"]
preparedVolume = prepareLinearSystem(numericalVolume; free_surface_spacing=recipeSpacing[1:2])
pointsOPT5 = getModelPoints(modelsOPT[1], 5,
    opt5Recipe["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
familyOPT5 = (models=modelsOPT, modelPoints=pointsOPT5, Δ=recipeSpacing, modelName="planewave_OPT5")
numericalOPT5 = numericalOperatorConstruction(Dict{String,Any}(
    "optRec" => opt5Recipe, "modelFam" => familyOPT5,
    "absorbingBoundaries" => nothing, "maskedRegionInSpace" => nothing,
    "boundaryConditions" => bcOPT, "representation" => "matrixfree",
))["numOperators"]
preparedOPT5Volume = prepareLinearSystem(numericalOPT5; free_surface_spacing=recipeSpacing[1:2])

# Additive free-surface traction closure -- the SAME construction as
# dev/kirishima_oct2026/stages/step2_opt_build.jl's applyFreeSurface branch,
# reused here unchanged on purpose: this is the exact thing under test.
boundaryParameters = copy(optParameters)
boundaryParameters["famousEquationType"] = "elasticTractionFree2D"
boundaryRecipe = cachedOPTRecipe(boundaryParameters, "elasticTractionFree2D_OPT3")
normalX = zeros(size(solidOPT)); normalZ = zeros(size(solidOPT))
for (point, normal) in zip(bcOPT.free_surface.points, bcOPT.free_surface.normals)
    normalX[point], normalZ[point] = normal
end
boundaryModels = [lambdaOPT, muOPT, normalX, normalZ]
boundaryPoints = getModelPoints(boundaryModels[1], 3,
    boundaryRecipe["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
boundaryFamily = (models=boundaryModels, modelPoints=boundaryPoints,
    Δ=recipeSpacing, modelName="planewave_free_surface")
numericalBoundary = numericalOperatorConstruction(Dict{String,Any}(
    "optRec" => boundaryRecipe, "modelFam" => boundaryFamily,
    "absorbingBoundaries" => (isnothing(bcOPT.cerjan) ? nothing : cerjan_padding(bcOPT.cerjan)),
    "maskedRegionInSpace" => bcOPT.free_surface.points,
    "representation" => "matrixfree",
))["numOperators"]
preparedBoundary = prepareLinearSystem(numericalBoundary)
surfaceWhole = numericalVolume.numericalOperators.left.geometry.freeSurfaceBoundary.points
preparedOPT = overlapBoundaryLinearSystem(preparedVolume, preparedBoundary, surfaceWhole;
    mode=:additive, boundary_weight=1.0)
preparedOPT5 = overlapBoundaryLinearSystem(preparedOPT5Volume, preparedBoundary, surfaceWhole;
    mode=:additive, boundary_weight=1.0)
@assert preparedOPT.NForceField == 2
@assert preparedOPT5.NForceField == 2
@show preparedOPT.spaceShape preparedOPT5.spaceShape
