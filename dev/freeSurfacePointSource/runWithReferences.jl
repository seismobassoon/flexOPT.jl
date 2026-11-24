# dev/freeSurfacePointSource/runWithReferences.jl -- the FULL point-force,
# flat-free-surface comparison: FD3, convFD3, OPT3, OPT5, SPECFEM2D.
#
# This extends run.jl (OPT3/OPT5 only, written for application/forBoris/,
# which deliberately skips FD3/SPECFEM2D so a collaborator doesn't need
# SPECFEM2D installed) with the reference solvers, to answer: how well do
# all of flexOPT's point-force/free-surface solvers actually agree with
# each other and with SPECFEM2D, right now, on the SAME physical problem?
#
# A faithful, trimmed port of notebooks/HomogeneousElastic2DBenchmark_FreeSurface.ipynb
# cells 3 (config), 6 (FD3), 8/9/10 (OPT3/OPT5 recipes/operators/propagation,
# reused verbatim from run.jl, PLUS the convFD3 portions of the same cells),
# 16 (receiver traces), 18 (SPECFEM2D -- here, reusing the already-cached
# SPECFEM2D run at data/specfem2d_benchmarks/homogeneous_flat_free_surface_dx500m_nrec23
# rather than re-running the mesher/solver, since that case already has the
# exact 5 stations this script needs: S0016/S0022/S0017/S0023/S0018 =
# x=-30,-20,-10,0,10 km @ z=-10 km -- confirmed against its DATA/STATIONS),
# and 20/21 (per-station waveform-fit metrics + the normalized multi-solver
# waveformFigure). OPT4, the predictor-corrector scheme, the shifted_full
# closure, the Takeuchi A90 experiment, the Gamma ablation, and the full
# 2-D receiverGrid diagnostics (parityAudit, phaseWindowAudit,
# absolute_component_figure) are all dropped, same as run.jl -- this script
# answers the "waveform fit between all of them" question, nothing more.
#
# Run:
#   julia --project=. dev/freeSurfacePointSource/runWithReferences.jl run
#   julia --project=. dev/freeSurfacePointSource/runWithReferences.jl run --full
# REPL:
#   julia> include("dev/freeSurfacePointSource/runWithReferences.jl")
#   julia> pointSourceWithReferences("run", "--full")

if !@isdefined(flexopt_root)
    import Pkg
    function find_flexopt_root(start_dir=@__DIR__)
        candidates = String[]
        haskey(ENV, "FLEXOPT_ROOT") && push!(candidates, abspath(expanduser(ENV["FLEXOPT_ROOT"])))
        directory = abspath(start_dir)
        while true
            push!(candidates, directory)
            parent = dirname(directory)
            parent == directory && break
            directory = parent
        end
        for candidate in unique(candidates)
            isfile(joinpath(candidate, "src", "flexOPT.jl")) && return candidate
        end
        error("Cannot locate flexOPT; set ENV[\"FLEXOPT_ROOT\"]")
    end
    global flexopt_root = find_flexopt_root()
    Pkg.activate(flexopt_root)
end
include(joinpath(flexopt_root, "dev", "common", "devBootstrap.jl"))

pointSourceWithReferences(args...) = (global APP_ARGS = collect(String, string.(args));
    include(joinpath(@__DIR__, "runWithReferences.jl")))

DEV = parseDevArgs(@isdefined(APP_ARGS) ? APP_ARGS : ARGS; steps=["run"], defaultStep="run")
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
DEV_STEP[] = "run_with_references"

# ---------------------------------------------------------------- bootstrap
pointSourceRefLoadEnv = !@isdefined(POINTSOURCEREF_ENV)
pointSourceRefLoadEnv && !@isdefined(commonBatchs) && include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
using .commonBatchs
pointSourceRefLoadEnv && !@isdefined(backend) && include(joinpath(flexopt_root, "src", "batchFiles", "batchGPU.jl"))
backend isa KernelAbstractions.CPU && @warn("No GPU backend is available; OPT recipes will be built on CPU")
pointSourceRefLoadEnv && !@isdefined(elasticWave2D) && include(joinpath(flexopt_root, "src", "elasticWave2D.jl"))
pointSourceRefLoadEnv && !@isdefined(flexOPT) && include(joinpath(flexopt_root, "src", "flexOPT.jl"))
pointSourceRefLoadEnv && !@isdefined(specfemBenchmark) && include(joinpath(flexopt_root, "src", "specfemBenchmark.jl"))
using .elasticWave2D, .flexOPT, .specfemBenchmark
set_default_form!(:weak)
@assert DEFAULT_FORM[] === :weak
using CairoMakie, Statistics, LinearAlgebra, SparseArrays, JLD2, DelimitedFiles
CairoMakie.activate!(type="png")
POINTSOURCEREF_ENV = true

quickRun = !DEV.full

# =========================================================== cell 3 (config)
referenceDx = 500.0
benchmarkDx = parse(Float64, get(ENV, "FLEXOPT_BENCHMARK_DX", string(referenceDx)))
fdOptSpatialRefinement = quickRun ? 1 : 2
dx = benchmarkDx / fdOptSpatialRefinement
optSpatialRefinement = 1
temporalRefinement = 2
sourcePosition = (x=-10e3, z=-20e3)
domainHalfWidth = parse(Float64, get(ENV, "FLEXOPT_DOMAIN_HALF_WIDTH", string(quickRun ? 60e3 : 90e3)))
x = collect((sourcePosition.x-domainHalfWidth):dx:(sourcePosition.x+domainHalfWidth))
applyFreeSurface = true
applyCerjan = false
z = collect((sourcePosition.z-domainHalfWidth):dx:(sourcePosition.z+domainHalfWidth))
nx, nz = length(x), length(z)
vp0, vs0, rho0 = 6000.0, 3464.0, 2700.0
material = applyFreeSurface ? repeat(reshape(z .<= 0.0, 1, :), nx, 1) : trues(nx, nz)
vpField = fill(vp0, nx, nz)
vsField = fill(vs0, nx, nz)
rhoField = fill(rho0, nx, nz)
model = (ρ=rhoField ./ 1e3, Vpv=vpField ./ 1e3, Vsv=vsField ./ 1e3)
vpMaximum, vsMinimum = vp0, vs0
bcFD = applyFreeSurface ? boundary_geometry(material, (dx, dx); cerjan=nothing) :
    BoundaryConditionSet(free_surface=nothing, cerjan=nothing,
        material_mask=BitArray(material), free_surface_mode=:pinned_void)
receiverZ = -10e3
receiverX = quickRun ? [-30e3, -20e3, -10e3, 0.0, 10e3] : [-20e3, -10e3, 0.0, 10e3, 20e3]
receiverGrid = [(x=xr, z=receiverZ) for xr in receiverX]   # trimmed: the line only, not the full 2-D grid
duration = parse(Float64, get(ENV, "FLEXOPT_BENCHMARK_DURATION", string(quickRun ? 9.0 : 14.0)))
outputSampling = 0.10
referencePointsPerSWavelength = 16
rickerSharpness = 1.5
minimumAcceptedPointsPerSWavelength = 10
coarsestSpacing = referenceDx
sourceFrequency = rickerSharpness * vs0 / (referencePointsPerSWavelength * coarsestSpacing)
sourceDelay = 1.5 / sourceFrequency
sourceForce = 1.0e10
sourceMode = :gaussian
sourceSpatialSigma = referenceDx
sourceCutoffSigmas = 2.0
function spatial_source_specs(xaxis, zaxis; center=sourcePosition,
    mode=sourceMode, sigma=sourceSpatialSigma, cutoff=sourceCutoffSigmas)
    if mode === :point
        return [(x=xaxis[argmin(abs.(xaxis .- center.x))],
            z=zaxis[argmin(abs.(zaxis .- center.z))], weight=1.0)]
    end
    specs = [(x=Float64(xv), z=Float64(zv),
        weight=exp(-((xv-center.x)^2 + (zv-center.z)^2) / (2sigma^2)))
        for xv in xaxis for zv in zaxis
        if hypot(xv-center.x, zv-center.z) <= cutoff*sigma]
    total = sum(s.weight for s in specs)
    [(x=s.x, z=s.z, weight=s.weight / total) for s in specs]
end
rickerSource(time) = begin
    a = π * sourceFrequency * (time - sourceDelay)
    (1 - 2a^2) * exp(-a^2)
end
pointsPerSWavelength = (FD3=vsMinimum / (sourceFrequency * dx),
    OPT3=vsMinimum / (sourceFrequency * (dx / optSpatialRefinement)),
    OPT5=vsMinimum / (sourceFrequency * (dx / optSpatialRefinement)))
distanceToNearestArtificialEdge = minimum((
    sourcePosition.x - first(x), last(x) - sourcePosition.x, sourcePosition.z - first(z)))
significantSourceStart = sourceDelay - 1 / sourceFrequency
earliestPEdgeArrival = significantSourceStart + distanceToNearestArtificialEdge / vpMaximum
probe10kmMinus5km = (x=10e3, z=-5e3)
function artificial_boundary_returns(station)
    imageRight = (x=2last(x)-sourcePosition.x, z=sourcePosition.z)
    imageLeft = (x=2first(x)-sourcePosition.x, z=sourcePosition.z)
    imageBottom = (x=sourcePosition.x, z=2first(z)-sourcePosition.z)
    arrival(image) = significantSourceStart + hypot(station.x-image.x, station.z-image.z) / vpMaximum
    (right=arrival(imageRight), left=arrival(imageLeft), bottom=arrival(imageBottom))
end
probeArtificialReturns = artificial_boundary_returns(probe10kmMinus5km)
@assert duration < earliestPEdgeArrival
@assert 8.0 < minimum(values(probeArtificialReturns))
@assert minimum(values(pointsPerSWavelength)) >= minimumAcceptedPointsPerSWavelength
@show quickRun applyFreeSurface (nx, nz) dx sourcePosition receiverZ receiverX
@show sourceFrequency sourceDelay sourceForce sourceMode sourceSpatialSigma

# ============================================================= cell 6 (FD3)
fdConfig = ElasticThreePointConfig2D(pointsInSpace=3, pointsInTime=3, supplementaryOrder=2, cfl=0.38)
fd = prepare_elastic_wave_2d(model, (dx, dx); material_mask=material, boundary_conditions=bcFD, config=fdConfig)
fdCoordinates = elastic_wave_coordinates(x, z, fd)
fdPadding = fd.padding
sourcePhysical = CartesianIndex(argmin(abs.(x .- sourcePosition.x)), argmin(abs.(z .- sourcePosition.z)))
sourceSpecsFD = spatial_source_specs(fdCoordinates.x, fdCoordinates.z)
sourceNodesFD = [(index=CartesianIndex(argmin(abs.(fdCoordinates.x .- s.x)),
    argmin(abs.(fdCoordinates.z .- s.z))), weight=s.weight) for s in sourceSpecsFD]
fdSteps = ceil(Int, duration / fd.dt)
fdOutputStride = max(1, round(Int, outputSampling / fd.dt))
fdFramesX = Matrix{Float32}[copy(fd.ux)]
fdFramesZ = Matrix{Float32}[copy(fd.uz)]
fdTimes = Float64[0.0]
@time for step in 1:fdSteps
    step_elastic_wave_2d!(fd)
    for sourceNode in sourceNodesFD
        elasticWave2D.add_ricker_source!(fd, sourceNode.index;
            f0=sourceFrequency, t0=sourceDelay, amplitude=sourceForce * sourceNode.weight,
            component=:z, source_kind=:force)
    end
    if step % fdOutputStride == 0 || step == fdSteps
        push!(fdFramesX, copy(fd.ux)); push!(fdFramesZ, copy(fd.uz)); push!(fdTimes, fd.time)
    end
end
uxFD = cat(fdFramesX...; dims=3)
uzFD = cat(fdFramesZ...; dims=3)
@show fd.dt size(uzFD) maximum(abs, uxFD) maximum(abs, uzFD)

# ================================================= cells 8/9 (OPT3/OPT5/convFD3 recipes+operators)
dxOPT = dx / optSpatialRefinement
xOPTPhysical = collect(first(x):dxOPT:last(x))
zOPTPhysical = collect(first(z):dxOPT:last(z))
solidOPT = applyFreeSurface ? repeat(reshape(zOPTPhysical .<= 0.0, 1, :),
    length(xOPTPhysical), 1) : trues(length(xOPTPhysical), length(zOPTPhysical))
dtOPT = 0.20 * dxOPT / (sqrt(2) * vpMaximum * temporalRefinement)
vpFieldOPT = fill(vp0, size(solidOPT))
vsFieldOPT = fill(vs0, size(solidOPT))
rhoFieldOPT = fill(rho0, size(solidOPT))
recipeSpacing = (1.0, 1.0, 1.0)
rhoOPT = rhoFieldOPT ./ rho0
muOPT = (rhoFieldOPT .* vsFieldOPT .^ 2 ./ rho0) .* (dtOPT / dxOPT)^2
lambdaOPT = (rhoFieldOPT .* (vpFieldOPT .^ 2 .- 2vsFieldOPT .^ 2) ./ rho0) .* (dtOPT / dxOPT)^2
muOPT[.!solidOPT] .= 0.0
lambdaOPT[.!solidOPT] .= 0.0
optParameters = Dict{String,Any}(
    "famousEquationType" => "2DsismoTimeIsoHeteroSingleForce",
    "Δ" => recipeSpacing, "orderBtime" => 1, "orderBspace" => 1,
    "pointsInSpace" => 3, "pointsInTime" => 3, "supplementaryOrder" => 2,
    "variationalForm" => DEFAULT_FORM[],
    "taylor_inverse_mode" => :weak_operator_optimized,
    "fieldItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "materItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "recipe_backend" => backend,
)
recipeCacheDirectory = "semiSymbolic"
recipeCacheVersion = 9
function cachedOPTRecipe(parameters, prefix; existing=nothing)
    cacheParameters = Dict{String,Any}(k => v for (k, v) in parameters if k != "recipe_backend")
    cacheParameters["recipe_cache_version"] = recipeCacheVersion
    function produceRecipe(config)
        runtimeParameters = Dict{String,Any}(config)
        pop!(runtimeParameters, "hash_id", nothing)
        pop!(runtimeParameters, "recipe_cache_version", nothing)
        runtimeParameters["recipe_backend"] = backend
        isnothing(existing) ? makeOPTsemiSymbolic(runtimeParameters) : existing
    end
    myProduceOrLoad(produceRecipe, cacheParameters, recipeCacheDirectory, prefix)
end
optRecipe = cachedOPTRecipe(optParameters, "elastic2D_OPT3")
clippedSurfaceParameters = copy(optParameters)
clippedSurfaceParameters["pointsInSpace"] = (3, 2)
clippedSurfaceParameters["nuCentre"] = (2, 2)
clippedSurfaceParameters["fieldItpl"] = (ptsSpace=(3, 2), ptsTime=1, offsetSpace=(0.0, 0.0),
    offsetTime=1, YorderBspace=1, YorderBtime=-1)
clippedSurfaceParameters["trialFunctionRefPoints"] = ((1, 2, 3), (1, 2, 3), (1, 2, 3))
clippedSurfaceParameters["testIntegrationBounds"] = ((1.0, 3.0), (1.0, 2.0), (1.0, 3.0))
clippedSurfaceParameters["exactTaylorTotalDegree"] = 2
clippedSurfaceRecipe = cachedOPTRecipe(clippedSurfaceParameters, "elastic2D_OPT3_surface_clipped_available")
opt5Parameters = copy(optParameters)
opt5Parameters["orderBspace"] = 1
opt5Parameters["pointsInSpace"] = 5
opt5Parameters["fieldItpl"] = merge(optParameters["fieldItpl"], (offsetSpace=2.0,))
opt5Parameters["materItpl"] = merge(optParameters["materItpl"], (offsetSpace=2.0,))
opt5Recipe = cachedOPTRecipe(opt5Parameters, "elastic2D_OPT5")
clippedSurfaceParametersOPT5 = copy(opt5Parameters)
clippedSurfaceParametersOPT5["pointsInSpace"] = (5, 3)
clippedSurfaceParametersOPT5["nuCentre"] = (3, 3)
clippedSurfaceParametersOPT5["fieldItpl"] = (ptsSpace=(5, 3), ptsTime=1, offsetSpace=(0.0, 0.0),
    offsetTime=1, YorderBspace=1, YorderBtime=-1)
clippedSurfaceParametersOPT5["trialFunctionRefPoints"] = ((1, 2, 3, 4, 5), (1, 2, 3, 4, 5), (1, 2, 3))
clippedSurfaceParametersOPT5["testIntegrationBounds"] = ((1.0, 5.0), (1.0, 3.0), (1.0, 3.0))
clippedSurfaceParametersOPT5["exactTaylorTotalDegree"] = 4
clippedSurfaceRecipeOPT5 = cachedOPTRecipe(clippedSurfaceParametersOPT5, "elastic2D_OPT5_surface_clipped_available")
conventionalFDParameters = copy(optParameters)
conventionalFDParameters["orderBspace"] = -1
conventionalFDParameters["orderBtime"] = -1
conventionalFDParameters["supplementaryOrder"] = 0
conventionalFDParameters["variationalForm"] = :strong
conventionalFDRecipe = cachedOPTRecipe(conventionalFDParameters, "elastic2D_convFD3")
let coefficientGatePath = joinpath(flexopt_root, "data", "elastic_lhs_coefficient_gate.jld2")
    isfile(coefficientGatePath) || error("Run ElasticLHSOperatorAudit.ipynb before OPT5 propagation")
    coefficientGate = load(coefficientGatePath)
    @assert all(result.passed for result in values(coefficientGate["coefficient_health"]))
end
bcOPT = applyFreeSurface ? boundary_geometry(solidOPT, (dxOPT, dxOPT);
    free_surface_mode=:pinned_void, cerjan=nothing) :
    BoundaryConditionSet(free_surface=nothing, cerjan=nothing,
        material_mask=BitArray(solidOPT), free_surface_mode=:pinned_void)
bcConvFD = applyFreeSurface ? boundary_geometry(solidOPT, (dxOPT, dxOPT);
    free_surface_mode=:zero_traction_flux, cerjan=nothing) : bcOPT
modelsOPT = [rhoOPT, lambdaOPT, muOPT]
pointsOPT = getModelPoints(modelsOPT[1], 3, optRecipe["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
familyOPT = (models=modelsOPT, modelPoints=pointsOPT, Δ=recipeSpacing, modelName="homogeneous_OPT3_scaled_coordinates")
numericalVolume = numericalOperatorConstruction(Dict{String,Any}(
    "optRec" => optRecipe, "modelFam" => familyOPT, "absorbingBoundaries" => nothing,
    "maskedRegionInSpace" => nothing, "boundaryConditions" => bcOPT, "representation" => "matrixfree"))["numOperators"]
preparedVolume = prepareLinearSystem(numericalVolume; free_surface_spacing=recipeSpacing[1:2])
pointsOPT5 = getModelPoints(modelsOPT[1], 5, opt5Recipe["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
familyOPT5 = (models=modelsOPT, modelPoints=pointsOPT5, Δ=recipeSpacing, modelName="homogeneous_OPT5_scaled_coordinates")
numericalOPT5 = numericalOperatorConstruction(Dict{String,Any}(
    "optRec" => opt5Recipe, "modelFam" => familyOPT5, "absorbingBoundaries" => nothing,
    "maskedRegionInSpace" => nothing, "boundaryConditions" => bcOPT, "representation" => "matrixfree"))["numOperators"]
preparedOPT5Volume = prepareLinearSystem(numericalOPT5; free_surface_spacing=recipeSpacing[1:2])
timedConvFDAssembly = numericalOperatorConstruction(Dict{String,Any}(
    "optRec" => conventionalFDRecipe, "modelFam" => familyOPT, "absorbingBoundaries" => nothing,
    "maskedRegionInSpace" => nothing, "boundaryConditions" => bcConvFD, "representation" => "matrixfree"))["numOperators"]
preparedConvFD = prepareLinearSystem(timedConvFDAssembly; free_surface_spacing=recipeSpacing[1:2])
@assert optRecipe["recette"].form === :weak
@assert optRecipe["recette"].variationalForm === :natural_weak
@assert conventionalFDRecipe["recette"].form === :strong
function prepare_surface_geometry(recipe, name; models=modelsOPT, points_in_space=3)
    surfacePoints = getModelPoints(models[1], points_in_space,
        recipe["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
    surfaceFamily = (models=models, modelPoints=surfacePoints, Δ=recipeSpacing, modelName=name)
    numerical = numericalOperatorConstruction(Dict{String,Any}(
        "optRec" => recipe, "modelFam" => surfaceFamily, "absorbingBoundaries" => nothing,
        "maskedRegionInSpace" => bcOPT.free_surface.points, "boundaryConditions" => nothing,
        "representation" => "matrixfree"))["numOperators"]
    prepareLinearSystem(numerical)
end
preparedClippedSurface = prepare_surface_geometry(clippedSurfaceRecipe, "surface_clipped_available")
surfaceWhole = numericalVolume.numericalOperators.left.geometry.freeSurfaceBoundary.points
preparedOPT = overlapBoundaryLinearSystem(preparedVolume, preparedClippedSurface, surfaceWhole; mode=:replace)
@assert preparedOPT.boundary_overlap_mode === :replace
@assert !hasproperty(preparedOPT, :boundary_future_operator)
surfaceWholeOPT5 = numericalOPT5.numericalOperators.left.geometry.freeSurfaceBoundary.points
preparedClippedSurfaceOPT5 = prepare_surface_geometry(clippedSurfaceRecipeOPT5,
    "surface_clipped_available_OPT5"; points_in_space=5)
preparedOPT5 = overlapBoundaryLinearSystem(preparedOPT5Volume, preparedClippedSurfaceOPT5, surfaceWholeOPT5; mode=:replace)
@assert preparedOPT5.boundary_overlap_mode === :replace
@assert !hasproperty(preparedOPT5, :boundary_future_operator)
@assert preparedOPT.NForceField == 2
@assert preparedOPT5.NForceField == 2
@assert preparedConvFD.NForceField == 2
function surface_audit(prepared, name)
    A = flexOPT.materializeConstantMatrix(prepared)
    rows = hasproperty(prepared, :boundary_overlap_rows) ? prepared.boundary_overlap_rows : Int[]
    rowNorms = [norm(A[row, :]) for row in rows]
    factorization = try
        lu(A)
        :regular
    catch err
        err isa SingularException ? :singular : rethrow()
    end
    audit = (name=name, factorization=factorization, boundary_row_norms=extrema(rowNorms),
        zero_boundary_rows=count(iszero, rowNorms))
    display(audit)
    @assert audit.factorization === :regular
end
surface_audit(preparedOPT, :OPT3)
surface_audit(preparedOPT5, :OPT5)
@assert applyFreeSurface
@assert length(bcOPT.free_surface.points) > 0
surfaceNormals = bcOPT.free_surface.normals
@assert all(n -> isapprox(n[1], 0.0; atol=1e-12), surfaceNormals)
@assert all(n -> n[2] > 0.99, surfaceNormals)

# ============================================= cell 10 (OPT3/OPT5/convFD3 propagation)
optPadding = zeros(Int, 2, 2)
xOPT = range(first(xOPTPhysical); step=dxOPT, length=preparedOPT.spaceShape[1])
zOPT = range(first(zOPTPhysical); step=dxOPT, length=preparedOPT.spaceShape[2])
sourceSpecsOPT = spatial_source_specs(xOPTPhysical, zOPTPhysical)
sourceNodesOPT = [(linear=LinearIndices(preparedOPT.spaceShape)[
        CartesianIndex(argmin(abs.(xOPT .- s.x)), argmin(abs.(zOPT .- s.z)))], weight=s.weight)
    for s in sourceSpecsOPT]
@assert isapprox(sum(s.weight for s in sourceNodesOPT), 1.0; atol=1e-14)
optSteps = ceil(Int, duration / dtOPT)
optOutputStride = max(1, round(Int, outputSampling / dtOPT))
nPastForceLevels = preparedOPT.timePointsUsedForOneStep - 1
sourceTimesOPT = ((1-nPastForceLevels):optSteps) .* dtOPT
@assert sourceTimesOPT[preparedOPT.timePointsUsedForOneStep] == dtOPT
waveletOPT = rickerSource.(sourceTimesOPT)
sourceScaleOPT = sourceForce * dtOPT^2 / (rho0 * dxOPT^2)
function assemble_source(prepared)
    s = zeros(Float64, prepared.NforcePoints, prepared.NForceField, length(sourceTimesOPT))
    for sourceNode in sourceNodesOPT
        s[sourceNode.linear, 2, :] .+= sourceNode.weight .* sourceScaleOPT .* waveletOPT
    end
    s
end
sourceOPT = assemble_source(preparedOPT)
@time propagationOPT = propagateLinearSystem(preparedOPT, optSteps, dtOPT;
    sourceFull=sourceOPT, output_stride=optOutputStride, blowup_limit=1e8,
    solver_name="OPT3 clipped-available", scheme=:direct)
@assert !propagationOPT.stopped_early
sourceOPT5 = assemble_source(preparedOPT5)
@time propagationOPT5 = propagateLinearSystem(preparedOPT5, optSteps, dtOPT;
    sourceFull=sourceOPT5, output_stride=optOutputStride, blowup_limit=1e8,
    solver_name="OPT5", scheme=:direct)
@assert !propagationOPT5.stopped_early
uxOPT, uzOPT, optTimes = propagationOPT.history[:, :, 1, :], propagationOPT.history[:, :, 2, :], propagationOPT.times
uxOPT5, uzOPT5, opt5Times = propagationOPT5.history[:, :, 1, :], propagationOPT5.history[:, :, 2, :], propagationOPT5.times
sourceConvFD = assemble_source(preparedConvFD)
convFDBoundaryGuard = 0.5 / sourceFrequency
convFDBoundarySafeTime = max(0.0, earliestPEdgeArrival - convFDBoundaryGuard)
convFDDuration = min(duration, convFDBoundarySafeTime)
convFDSteps = ceil(Int, convFDDuration / dtOPT)
@time propagationConvFD = propagateLinearSystem(preparedConvFD, convFDSteps, dtOPT;
    sourceFull=sourceConvFD, output_stride=optOutputStride, blowup_limit=1e8,
    solver_name="convFD3 (orderB=-1)", scheme=:direct)
uxConvFD = propagationConvFD.history[:, :, 1, :]
uzConvFD = propagationConvFD.history[:, :, 2, :]
convFDTimes = propagationConvFD.times
propagationConvFD.stopped_early && @warn("convFD3 stopped early; its partial history remains available",
    final_time=last(convFDTimes), requested_duration=duration)
@show convFDBoundarySafeTime convFDDuration propagationConvFD.stopped_early
@show dtOPT size(uzOPT) maximum(abs, uxOPT) maximum(abs, uzOPT)
@show size(uzOPT5) maximum(abs, uxOPT5) maximum(abs, uzOPT5)
@show size(uzConvFD) maximum(abs, uxConvFD) maximum(abs, uzConvFD)

# ======================================================= cell 16 (traces)
function displacement_traces(history, times, xaxis, zindex, receivers)
    traces = Matrix{Float64}(undef, length(times), length(receivers))
    for (j, receiver) in enumerate(receivers)
        ix = argmin(abs.(xaxis .- receiver))
        traces[:, j] .= Float64.(history[ix, zindex, :])
    end
    (time=Float64.(times), values=traces)
end
fdReceiverZIndex = argmin(abs.(fdCoordinates.z .- receiverZ))
optReceiverZIndex = argmin(abs.(zOPT .- receiverZ))
fdDisplacementZ = displacement_traces(uzFD, fdTimes, fdCoordinates.x, fdReceiverZIndex, receiverX)
opt3DisplacementZ = displacement_traces(uzOPT, optTimes, xOPT, optReceiverZIndex, receiverX)
opt5DisplacementZ = displacement_traces(uzOPT5, opt5Times, xOPT, optReceiverZIndex, receiverX)
convFDDisplacementZ = displacement_traces(uzConvFD, convFDTimes, xOPT, optReceiverZIndex, receiverX)

# =============================================== cell 18 (SPECFEM2D reference)
# Reuse the already-cached run at dx=500 m with 23 receivers (confirmed to
# include exactly these 5 stations at this receiverZ -- see DATA/STATIONS)
# rather than re-meshing/re-solving SPECFEM2D for this script.
specfemOutput = joinpath(flexopt_root, "data", "specfem2d_benchmarks",
    "homogeneous_flat_free_surface_dx500m_nrec23", "OUTPUT_FILES")
specfemStationNames = ["S0016", "S0022", "S0017", "S0023", "S0018"]  # x=-30,-20,-10,0,10 km @ z=-10km
specfemTimeShift = 1.2 / sourceFrequency   # external-STF convention, see src/specfemBenchmark.jl
specfemFiles = find_specfem2d_traces(specfemOutput; component=:z, seismotype=2, network="FX", stations=specfemStationNames)
@assert length(specfemFiles) == 5
# find_specfem2d_traces (src/specfemBenchmark.jl) always sorts its RETURN by
# station-name field, regardless of the order `stations` was given in -- so
# the 5 paths above come back alphabetically (S0016,S0017,S0018,S0022,S0023,
# i.e. x=-30,-10,10,-20,0 km), NOT in specfemStationNames/receiverX order
# (x=-30,-20,-10,0,10 km). Re-order explicitly rather than assume order
# preservation (this was wrong in an earlier version of this script and
# silently scrambled every station past the first).
specfemFilesByStation = Dict(split(basename(path), '.')[2] => path for path in specfemFiles)
orderedSpecfemFiles = [specfemFilesByStation[name] for name in specfemStationNames]
specfemRawTracesZ = [read_specfem2d_trace(file; time_shift=specfemTimeShift) for file in orderedSpecfemFiles]
function integrate_velocity_trace(trace)
    displacement = zeros(Float64, length(trace.values))
    displacement[2:end] .= cumsum(((trace.values[1:end-1] .+ trace.values[2:end]) ./ 2) .* diff(trace.time))
    (time=Float64.(trace.time), values=displacement)
end
specfemWaveformsZ = integrate_velocity_trace.(specfemRawTracesZ)

# ======================================= SeismicWavePC.jl (isolated env, read-only)
# Produced by dev/freeSurfacePointSource/seismicWavePCEnv/runSeismicWavePC.jl in
# its OWN isolated Julia environment (see that script's header: SeismicWavePC.jl
# pulls in a GeophysicalModelGenerator version that conflicts with flexOPT's own
# pinned Makie/GLMakie, so it cannot live in this shared environment). The saved
# JLD2 holds only plain NamedTuples of Float64 vectors, so it loads fine here
# without SeismicWavePC itself being installed in this environment.
seismicWavePCPath = joinpath(@__DIR__, "outputs", "seismicWavePC", "seismicwavepc_traces.jld2")
haveSeismicWavePC = isfile(seismicWavePCPath)
if haveSeismicWavePC
    seismicWavePCRawTraces = load(seismicWavePCPath, "traces")
    @assert length(seismicWavePCRawTraces) == length(receiverX)
    # Two conventions do not translate from this benchmark's other solvers
    # (see application/forBoris/README.md and runSeismicWavePC.jl's header):
    # force-amplitude normalization, and the Ricker source-time-function's
    # exact time-centering. Rather than guess SeismicWavePC's internal STF
    # centering formula, find the time lag that best cross-correlates its
    # trace at the station directly above the source against the SPECFEM2D
    # reference there, and apply that SAME lag (a pure STF-convention offset,
    # not a propagation effect, so it must be station-independent) everywhere.
    alignmentStation = something(findfirst(==(sourcePosition.x), receiverX), cld(length(receiverX), 2))
    function best_alignment_lag(reference, candidate; max_lag=4.0, step=0.02)
        best_lag, best_corr = 0.0, -Inf
        for lag in -max_lag:step:max_lag
            shifted = (time=candidate.time .+ lag, values=candidate.values)
            corr = try
                waveform_metrics(reference, shifted; samples=501).correlation
            catch
                -Inf
            end
            corr > best_corr && ((best_corr, best_lag) = (corr, lag))
        end
        (lag=best_lag, correlation=best_corr)
    end
    seismicWavePCAlignment = best_alignment_lag(
        specfemWaveformsZ[alignmentStation], seismicWavePCRawTraces[alignmentStation])
    @show alignmentStation seismicWavePCAlignment
    # Per-station amplitude match to the SPECFEM2D peak at that station (same
    # spirit as the forBoris README's suggested comparison strategy: rescale
    # by the ratio of peak amplitudes so relative_error/relative_rmse mean
    # something despite SeismicWavePC's own, non-physical M0 bookkeeping;
    # correlation is scale-invariant and unaffected either way).
    seismicWavePCDisplacementZ = map(eachindex(receiverX)) do station
        raw = seismicWavePCRawTraces[station]
        reference = specfemWaveformsZ[station]
        peakReference = max(maximum(abs, reference.values), eps(Float64))
        peakCandidate = max(maximum(abs, raw.values), eps(Float64))
        (time=raw.time .+ seismicWavePCAlignment.lag, values=raw.values .* (peakReference / peakCandidate))
    end
else
    @warn("SeismicWavePC traces not found at $seismicWavePCPath -- run seismicWavePCEnv/runSeismicWavePC.jl first; continuing without it")
end

# ================================================ cells 20/21 (waveform fit)
function waveform_metrics_compact(reference, candidate)
    full = waveform_metrics(reference, candidate; samples=1001)
    (correlation=full.correlation, relative_error=full.relative_error,
        relative_rmse=full.relative_rmse)
end
function solver_traces_for_station(station)
    base = (FD3=(time=fdDisplacementZ.time, values=fdDisplacementZ.values[:, station]),
        convFD3=(time=convFDDisplacementZ.time, values=convFDDisplacementZ.values[:, station]),
        OPT3=(time=opt3DisplacementZ.time, values=opt3DisplacementZ.values[:, station]),
        OPT5=(time=opt5DisplacementZ.time, values=opt5DisplacementZ.values[:, station]))
    haveSeismicWavePC ? merge(base, (SeismicWavePC=seismicWavePCDisplacementZ[station],)) : base
end
stationErrorMetrics = map(eachindex(receiverX)) do station
    reference = specfemWaveformsZ[station]
    candidates = solver_traces_for_station(station)
    metrics = map(trace -> waveform_metrics_compact(reference, trace), candidates)
    (station=station, x_km=receiverX[station] / 1e3, metrics=metrics)
end
foreach(display, stationErrorMetrics)
methodNames = keys(first(stationErrorMetrics).metrics)
globalErrorMetrics = NamedTuple{methodNames}(map(methodNames) do method
    rows = [getproperty(item.metrics, method) for item in stationErrorMetrics]
    (mean_correlation=mean(getproperty.(rows, :correlation)),
        rms_relative_error=sqrt(mean(abs2, getproperty.(rows, :relative_error))),
        maximum_relative_error=maximum(getproperty.(rows, :relative_error)))
end)
display(globalErrorMetrics)

waveformFigure = Figure(size=(1100, 230 * length(receiverX)))
for station in eachindex(receiverX)
    axis = Axis(waveformFigure[station, 1]; ylabel="x=$(receiverX[station]/1e3) km",
        xlabel=station == lastindex(receiverX) ? "time (s)" : "")
    traces = merge(solver_traces_for_station(station), (SPECFEM2D=specfemWaveformsZ[station],))
    for (name, trace) in pairs(traces)
        scale = max(maximum(abs, trace.values), eps(Float64))
        if name === :SPECFEM2D
            lines!(axis, trace.time, trace.values ./ scale; label=string(name), color=:black, linewidth=2)
        elseif name === :SeismicWavePC
            lines!(axis, trace.time, trace.values ./ scale; label=string(name), linestyle=:dot, linewidth=2)
        else
            lines!(axis, trace.time, trace.values ./ scale; label=string(name))
        end
    end
    station == 1 && axislegend(axis; position=:rb, nbanks=3)
end
Label(waveformFigure[0, 1], "Point force, flat free surface, u_z (normalized): FD3 / convFD3 / OPT3 / OPT5 / SPECFEM2D" *
    (haveSeismicWavePC ? " / SeismicWavePC" : ""); fontsize=14)

# ---------------------------------------------------------------- outputs
outdir = joinpath(DEV_OUTPUT_DIR[], DEV_STEP[])
mkpath(outdir)
if haveSeismicWavePC
    jldsave(joinpath(outdir, "all_solvers_traces.jld2");
        vp0, vs0, rho0, sourcePosition, sourceForce, sourceFrequency, sourceDelay,
        dx, dxOPT, dtOPT, receiverX, receiverZ,
        fdDisplacementZ, convFDDisplacementZ, opt3DisplacementZ, opt5DisplacementZ, specfemWaveformsZ,
        haveSeismicWavePC, seismicWavePCDisplacementZ, seismicWavePCAlignment,
        stationErrorMetrics, globalErrorMetrics)
else
    jldsave(joinpath(outdir, "all_solvers_traces.jld2");
        vp0, vs0, rho0, sourcePosition, sourceForce, sourceFrequency, sourceDelay,
        dx, dxOPT, dtOPT, receiverX, receiverZ,
        fdDisplacementZ, convFDDisplacementZ, opt3DisplacementZ, opt5DisplacementZ, specfemWaveformsZ,
        haveSeismicWavePC,
        stationErrorMetrics, globalErrorMetrics)
end
devSave(waveformFigure, "waveform_fit_all_solvers")
println("DONE")
