# dev/freeSurfacePointSource/run.jl -- OPT3 / OPT5 only, point-force (vertical,
# Gaussian-smoothed), flat free surface, homogeneous half-space.
#
# This is a faithful, TRIMMED port of
# notebooks/HomogeneousElastic2DBenchmark_FreeSurface.ipynb (cells 3, 8, 9, 10,
# 16 -- read in full before writing this), kept for exactly one purpose:
# `application/forBoris/` gives a collaborator (Boris Kaus, SeismicWavePC.jl)
# a clean OPT3/OPT5 reference to benchmark his own code against, without
# needing SPECFEM2D, FD3, or any of the notebook's own diagnostic/audit
# machinery for OPT4, the predictor-corrector marching scheme, the
# shifted_full surface closure, the Takeuchi & Geller A90 mass-hybrid
# experiment, or the Gamma-redistribution ablation (runOPTWithoutGamma).
# Those all still exist, validated, in the notebook; this script only runs
# the two solvers and the one free-surface closure (:clipped_available,
# :direct marching) the notebook itself prefers by default.
#
# The semi-symbolic OPT3/OPT5 recipes and their surface closures are cached
# in data/semiSymbolic/ under the SAME parameter dicts the notebook uses, so
# a first run here reuses whatever the notebook has already built there
# instead of recomputing from scratch (confirmed on disk 2026-11-24:
# elastic2D_OPT3*, elastic2D_OPT3_surface_clipped_available*,
# elastic2D_OPT5*, elastic2D_OPT5_surface_clipped_available* already cached).
#
# Kept from the notebook as correctness guards (cheap, no propagation):
# the coefficient-health gate (data/elastic_lhs_coefficient_gate.jld2),
# the free-surface factorization audit (must be :regular, not :singular),
# and the boundary/interior same-field mass ratio (printed, not asserted --
# the notebook's own comment: a ratio far from OPT3's flags OPT5's surface
# node as under-massed, the signature behind "OPT5 arrives early").
#
# Run:
#   julia --project=. dev/freeSurfacePointSource/run.jl run
#   julia --project=. dev/freeSurfacePointSource/run.jl run --full
# REPL:
#   julia> include("dev/freeSurfacePointSource/run.jl")
#   julia> pointSource("run", "--full")

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

pointSource(args...) = (global APP_ARGS = collect(String, string.(args)); include(joinpath(@__DIR__, "run.jl")))

DEV = parseDevArgs(@isdefined(APP_ARGS) ? APP_ARGS : ARGS; steps=["run"], defaultStep="run")
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
DEV_STEP[] = "run"

# ---------------------------------------------------------------- bootstrap
pointSourceLoadEnv = !@isdefined(POINTSOURCE_ENV)
pointSourceLoadEnv && !@isdefined(commonBatchs) && include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
using .commonBatchs
pointSourceLoadEnv && !@isdefined(backend) && include(joinpath(flexopt_root, "src", "batchFiles", "batchGPU.jl"))
backend isa KernelAbstractions.CPU && @warn("No GPU backend is available; OPT recipes will be built on CPU")
pointSourceLoadEnv && !@isdefined(flexOPT) && include(joinpath(flexopt_root, "src", "flexOPT.jl"))
using .flexOPT
set_default_form!(:weak)
@assert DEFAULT_FORM[] === :weak
using CairoMakie, Statistics, LinearAlgebra, SparseArrays, JLD2, DelimitedFiles
CairoMakie.activate!(type="png")
POINTSOURCE_ENV = true

quickRun = !DEV.full

# =========================================================== cell 3 (config)
# Physical setup copied byte-for-byte from the notebook's cell 3, dropping
# only the FD3/convFD3/SPECFEM2D-specific pieces (fdConfig, sourceSpecsFD,
# waveformQuantity branch that touches FD3) and the full 2-D receiverGrid
# (only ever used by the FD3/SPECFEM2D comparison cells; every assertion
# below that matters is already expressed with the fixed probe point
# `probe10kmMinus5km`, independent of receiverGrid -- confirmed by reading
# cell 3 in full).
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
vpMaximum, vsMinimum = vp0, vs0
receiverZ = -10e3
receiverX = quickRun ? [-30e3, -20e3, -10e3, 0.0, 10e3] : [-20e3, -10e3, 0.0, 10e3, 20e3]
duration = parse(Float64, get(ENV, "FLEXOPT_BENCHMARK_DURATION", string(quickRun ? 9.0 : 14.0)))
outputSampling = 0.10
referencePointsPerSWavelength = 16
rickerSharpness = 1.5
minimumAcceptedPointsPerSWavelength = 10
coarsestSpacing = referenceDx
sourceFrequency = rickerSharpness * vs0 / (referencePointsPerSWavelength * coarsestSpacing)
sourceDelay = 1.5 / sourceFrequency
sourceForce = 1.0e10          # N/m: total vertical line force in a unit-thickness slice
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
pointsPerSWavelength = (OPT3=vsMinimum / (sourceFrequency * (dx / optSpatialRefinement)),
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
@show pointsPerSWavelength earliestPEdgeArrival distanceToNearestArtificialEdge

# ======================================================= cell 8 (OPT recipes)
dxOPT = dx / optSpatialRefinement
xOPTPhysical = collect(first(x):dxOPT:last(x))
zOPTPhysical = collect(first(z):dxOPT:last(z))
solidOPT = applyFreeSurface ? repeat(reshape(zOPTPhysical .<= 0.0, 1, :),
    length(xOPTPhysical), 1) : trues(length(xOPTPhysical), length(zOPTPhysical))
dtOPT = 0.20 * dxOPT / (sqrt(2) * vpMaximum * temporalRefinement)
vpFieldOPT = fill(vp0, size(solidOPT))
vsFieldOPT = fill(vs0, size(solidOPT))
rhoFieldOPT = fill(rho0, size(solidOPT))
recipeSpacing = (1.0, 1.0, 1.0)      # :scaled_coordinates, the notebook's benchmark default
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
let coefficientGatePath = joinpath(flexopt_root, "data", "elastic_lhs_coefficient_gate.jld2")
    isfile(coefficientGatePath) || error("Run ElasticLHSOperatorAudit.ipynb before OPT5 propagation")
    coefficientGate = load(coefficientGatePath)
    @assert all(result.passed for result in values(coefficientGate["coefficient_health"]))
end

# ====================================================== cell 9 (operators)
bcOPT = applyFreeSurface ? boundary_geometry(solidOPT, (dxOPT, dxOPT);
    free_surface_mode=:pinned_void, cerjan=nothing) :
    BoundaryConditionSet(free_surface=nothing, cerjan=nothing,
        material_mask=BitArray(solidOPT), free_surface_mode=:pinned_void)
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
@assert optRecipe["recette"].form === :weak
@assert optRecipe["recette"].variationalForm === :natural_weak
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
function surface_audit(prepared, surfacePoints, name)
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
    return A
end
A_surface = surface_audit(preparedOPT, surfaceWhole, :OPT3)
A_surfaceOPT5 = surface_audit(preparedOPT5, surfaceWholeOPT5, :OPT5)
function same_field_row_sum(A, prepared, point; field=1)
    linear = LinearIndices(prepared.spaceShape)[point]
    row = LinearIndices((prepared.NField, prepared.NpointsSpace))[field, linear]
    unknown = CartesianIndices((prepared.NpointsSpace, prepared.NField))
    total = 0.0
    for (column, value) in zip(findnz(A[row, :])...)
        _, inputField = Tuple(unknown[column])
        inputField == field && (total += value)
    end
    total
end
interiorProbe = CartesianIndex(argmin(abs.(xOPTPhysical .- sourcePosition.x)), argmin(abs.(zOPTPhysical .- sourcePosition.z)))
surfaceProbe = surfaceWhole[cld(length(surfaceWhole), 2)]
boundaryMassAuditOPT3 = (boundary_mass=same_field_row_sum(A_surface, preparedOPT, surfaceProbe),
    interior_mass=same_field_row_sum(A_surface, preparedOPT, interiorProbe))
surfaceProbeOPT5 = surfaceWholeOPT5[cld(length(surfaceWholeOPT5), 2)]
boundaryMassAuditOPT5 = (boundary_mass=same_field_row_sum(A_surfaceOPT5, preparedOPT5, surfaceProbeOPT5),
    interior_mass=same_field_row_sum(A_surfaceOPT5, preparedOPT5, interiorProbe))
@show boundaryMassAuditOPT3 boundaryMassAuditOPT3.boundary_mass/boundaryMassAuditOPT3.interior_mass
@show boundaryMassAuditOPT5 boundaryMassAuditOPT5.boundary_mass/boundaryMassAuditOPT5.interior_mass
@assert applyFreeSurface
@assert length(bcOPT.free_surface.points) > 0
surfaceNormals = bcOPT.free_surface.normals
@assert all(n -> isapprox(n[1], 0.0; atol=1e-12), surfaceNormals)
@assert all(n -> n[2] > 0.99, surfaceNormals)

# ===================================================== cell 10 (propagation)
optPadding = zeros(Int, 2, 2)      # no Cerjan layer (applyCerjan=false)
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
sourceScaleOPT = sourceForce * dtOPT^2 / (rho0 * dxOPT^2)     # :scaled_coordinates convention
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
@show dtOPT size(uzOPT) maximum(abs, uxOPT) maximum(abs, uzOPT)
@show size(uzOPT5) maximum(abs, uxOPT5) maximum(abs, uzOPT5)

# ======================================================= cell 16 (traces)
function displacement_traces(history, times, xaxis, zindex, receivers)
    traces = Matrix{Float64}(undef, length(times), length(receivers))
    for (j, receiver) in enumerate(receivers)
        ix = argmin(abs.(xaxis .- receiver))
        traces[:, j] .= Float64.(history[ix, zindex, :])
    end
    (time=Float64.(times), values=traces)
end
optReceiverZIndex = argmin(abs.(zOPT .- receiverZ))
opt3DisplacementX = displacement_traces(uxOPT, optTimes, xOPT, optReceiverZIndex, receiverX)
opt3DisplacementZ = displacement_traces(uzOPT, optTimes, xOPT, optReceiverZIndex, receiverX)
opt5DisplacementX = displacement_traces(uxOPT5, opt5Times, xOPT, optReceiverZIndex, receiverX)
opt5DisplacementZ = displacement_traces(uzOPT5, opt5Times, xOPT, optReceiverZIndex, receiverX)

# ---------------------------------------------------------------- outputs
outdir = joinpath(DEV_OUTPUT_DIR[], DEV_STEP[])
mkpath(outdir)
jldsave(joinpath(outdir, "forBoris_opt3_opt5_traces.jld2");
    vp0, vs0, rho0, sourcePosition, sourceForce, sourceFrequency, sourceDelay,
    sourceMode, sourceSpatialSigma, dx, dxOPT, dtOPT, receiverX, receiverZ,
    opt3DisplacementX, opt3DisplacementZ, opt5DisplacementX, opt5DisplacementZ)
header = vcat(["time"], ["OPT3_ux_x=$(round(Int,r/1e3))km" for r in receiverX],
    ["OPT3_uz_x=$(round(Int,r/1e3))km" for r in receiverX],
    ["OPT5_ux_x=$(round(Int,r/1e3))km" for r in receiverX],
    ["OPT5_uz_x=$(round(Int,r/1e3))km" for r in receiverX])
writedlm(joinpath(outdir, "forBoris_opt3_opt5_traces.csv"),
    vcat(permutedims(header), hcat(opt3DisplacementZ.time, opt3DisplacementX.values,
        opt3DisplacementZ.values, opt5DisplacementX.values, opt5DisplacementZ.values)), ',')
println("saved -> ", joinpath(outdir, "forBoris_opt3_opt5_traces.jld2"))
println("saved -> ", joinpath(outdir, "forBoris_opt3_opt5_traces.csv"))

# Self-contained figure caption (requested after an external collaborator
# could not reproduce an earlier plot from the image alone): state the
# model, source geometry/"incidence", and wavelet directly on the figure,
# not only in README.md, so the plot stands on its own.
waveformCaptionLine1 = "Point-force, flat free surface, u_z: OPT3 vs OPT5"
waveformCaptionLine2 = "Model: vp=$(round(Int,vp0)) m/s, vs=$(round(Int,vs0)) m/s, " *
    "rho=$(round(Int,rho0)) kg/m3 (homogeneous half-space, flat free surface at z=0)"
waveformCaptionLine3 = "Source: vertical point force, embedded at x=$(round(Int,sourcePosition.x/1e3)) km, " *
    "z=$(round(Int,sourcePosition.z/1e3)) km depth (NOT at the surface, so there is no single " *
    "incidence angle -- receivers span near-vertical, directly above the source, to oblique, at " *
    "the farthest stations); Gaussian-smoothed in space, sigma=$(round(Int,sourceSpatialSigma)) m"
waveformCaptionLine4 = "Wavelet: Ricker, f0=$(round(sourceFrequency;digits=4)) Hz, delay=$(round(sourceDelay;digits=3)) s " *
    "-- (1-2a^2)exp(-a^2), a=pi*f0*(t-delay); total force=$(sourceForce) N/m (2-D line force)"
waveformFigure = Figure(size=(1000, 170 * length(receiverX) + 110))
for (i, xr) in enumerate(receiverX)
    ax = Axis(waveformFigure[i + 1, 1]; ylabel="x=$(round(Int, xr/1e3)) km",
        xlabel=(i == length(receiverX) ? "time (s)" : ""))
    lines!(ax, opt3DisplacementZ.time, opt3DisplacementZ.values[:, i]; label="OPT3", color=:black)
    lines!(ax, opt5DisplacementZ.time, opt5DisplacementZ.values[:, i]; label="OPT5", color=:crimson, linestyle=:dash)
    i == 1 && axislegend(ax; position=:rt)
end
Label(waveformFigure[0, 1], waveformCaptionLine1; fontsize=14, font=:bold)
Label(waveformFigure[1, 1], waveformCaptionLine2 * "\n" * waveformCaptionLine3 * "\n" * waveformCaptionLine4;
    fontsize=10, justification=:left, halign=:left)
devSave(waveformFigure, "opt3_vs_opt5_waveforms")
println("DONE")
