# dev/cutCells/freeEnd1D.jl -- step `freeEnd1D` of dev/cutCells/run.jl
#
# 1-D elastic bar (famous equation "1DsismoTime", OPT3, weak form) whose free
# end x_s falls a fraction η of a cell beyond the last solid node. Point force
# (a few nodes, Gaussian weights) 40 Δx before the end, receivers before and
# behind the source. Reference: the exact solution, direct wave + even image
# (free end), u = F/(2ρc) Σ_j w_j [S(t - |x - x_j|/c) + S(t - |x - x_j*|/c)],
# S the running integral of the Ricker wavelet.
#
# What is checked (all printed, figures under outputs/freeEnd1D/):
#   1. a node with an uncut support gives back the volume row (self-check);
#   2. waveform fit to the exact solution, per receiver;
#   3. timing error of the reflection (cross-correlation against the exact
#      reflected pulse), compared with the 2ηΔx/c a stair-step end would make;
#   4. reflected/direct amplitude ratio (convention-free; exact value 1).
#
# Units: Δx = 1, c = 1, ρ = 1 (scaled coordinates, as the 2-D benchmarks).

if !@isdefined(CUTCELLS_FULL_ENV)
    !@isdefined(commonBatchs) && include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
    @eval using .commonBatchs
    !@isdefined(backend) && include(joinpath(flexopt_root, "src", "batchFiles", "batchGPU.jl"))
    !@isdefined(flexOPT) && include(joinpath(flexopt_root, "src", "flexOPT.jl"))
    @eval using .flexOPT: set_default_form!, DEFAULT_FORM, makeOPTsemiSymbolic, getModelPoints,
        numericalOperatorConstruction, prepareLinearSystem, propagateLinearSystem
    global CUTCELLS_FULL_ENV = true
end
flexOPT.set_default_form!(:weak)

const FE = (dx=1.0, c=1.0, rho=1.0, pointsPerWavelength=16, sourceDistance=40.0,
    leftMargin=120.0, voidNodes=6, receiverDistances=(5.0, 20.0, 60.0), cfl=0.1,
    force=1.0, sourceSigma=1.0, sourceCutoff=3.0)
etaOfCase(case) = Dict("eta0" => 0.0, "eta025" => 0.25, "eta05" => 0.5, "eta075" => 0.75)[case]

f0 = FE.c / (FE.pointsPerWavelength * FE.dx)
t0 = 1.5 / f0
duration = t0 + (FE.sourceDistance + maximum(FE.receiverDistances)) / FE.c + 1.2 / f0
dt = FE.cfl * FE.dx / FE.c
ricker(t) = (a = π * f0 * (t - t0); (1 - 2a^2) * exp(-a^2))
rickerIntegral(t) = (a = π * f0 * (t - t0); (t - t0) * exp(-a^2))      # S' = ricker

# ------------------------------------------------------------------ recipe
optParameters1D = Dict{String,Any}(
    "famousEquationType" => "1DsismoTime",
    "Δ" => (1.0, 1.0), "orderBtime" => 1, "orderBspace" => 1,
    "pointsInSpace" => 3, "pointsInTime" => 3, "supplementaryOrder" => 2,
    "variationalForm" => DEFAULT_FORM[],
    "taylor_inverse_mode" => :weak_operator_optimized,
    "fieldItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "materItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "recipe_backend" => backend,
)
function cachedRecipe1D(parameters, prefix)
    cacheParameters = Dict{String,Any}(k => v for (k, v) in parameters if k != "recipe_backend")
    cacheParameters["recipe_cache_version"] = 9
    myProduceOrLoad(cacheParameters, "semiSymbolic", prefix) do config
        runtime = Dict{String,Any}(config)
        pop!(runtime, "hash_id", nothing)
        pop!(runtime, "recipe_cache_version", nothing)
        runtime["recipe_backend"] = backend
        flexOPT.makeOPTsemiSymbolic(runtime)
    end
end
recipe1D = cachedRecipe1D(optParameters1D, "elastic1D_OPT3")
scheme1D = flexOPT.CutCellScheme(optParameters1D)
@show scheme1D.maxDegree scheme1D.pointsShape

# ------------------------------------------------------------------ cases
resultsFE = Dict{String,Any}()
for case in DEV.cases
    η = etaOfCase(case)
    println("\n== free end between nodes, η = $η")
    xs = 0.0 + η * FE.dx                          # free end; last solid node at x = 0
    x = collect(-(FE.sourceDistance + FE.leftMargin):FE.dx:(FE.voidNodes * FE.dx))
    nx = length(x)
    ls = flexOPT.levelSetFromFunction(r -> r[1] - xs, (x,))
    geom = flexOPT.CutCellGeometry(ls; maxDegree=scheme1D.maxDegree)
    solid = flexOPT.sideMask(ls)
    cut = flexOPT.cutNodes(geom)
    aggregation = flexOPT.sliverAggregation(geom)
    farVoid = [I for I in CartesianIndices(ls.φ) if !solid[I] && !(I in cut)]
    @show length(cut) [x[I[1]] for I in cut] [(x[a.node[1]], x[a.target[1]]) for a in aggregation]

    # material: homogeneous, extended over the void (the cut rows integrate the solid part only)
    rhoS = fill(FE.rho / FE.rho, nx)
    muS = fill(FE.rho * FE.c^2 * (dt / FE.dx)^2 / FE.rho, nx)
    models = [rhoS, muS]
    timeMarching = recipe1D["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching
    family = (models=models, modelPoints=getModelPoints(models[1], 3, timeMarching),
        Δ=(1.0, 1.0), modelName="cutcells_freeEnd1D_homogeneous")
    numerical = numericalOperatorConstruction(Dict{String,Any}(
        "optRec" => recipe1D, "modelFam" => family, "absorbingBoundaries" => nothing,
        "maskedRegionInSpace" => nothing, "boundaryConditions" => nothing,
        "representation" => "matrixfree"))["numOperators"]
    volume = prepareLinearSystem(numerical)
    volume.spaceShape == (nx,) || error("operator grid $(volume.spaceShape) differs from the model grid ($nx,)")

    probe = CartesianIndex(argmin(abs.(x .+ FE.sourceDistance .+ 20)))
    selfCheck = flexOPT.cutCellSelfCheck(scheme1D, volume, models, probe)
    @show selfCheck
    maximum(values(selfCheck)) < 1e-8 || error("cut-cell rows of an uncut node differ from the volume rows: $selfCheck")

    rows = flexOPT.cutCellRows(scheme1D, geom, cut, models, volume)
    system = flexOPT.installRows(volume, rows)
    system = flexOPT.applyAggregation(system, aggregation)
    system = flexOPT.pinRows(system, farVoid)

    # source: Gaussian weights on nodes around x_src
    xsrc = xs - FE.sourceDistance
    weights = [exp(-(xi - xsrc)^2 / (2FE.sourceSigma^2)) for xi in x]
    weights[abs.(x .- xsrc) .> FE.sourceCutoff * FE.sourceSigma] .= 0.0
    weights ./= sum(weights)
    steps = ceil(Int, duration / dt)
    nPast = system.timePointsUsedForOneStep - 1
    sourceTimes = ((1-nPast):steps) .* dt
    source = zeros(system.NforcePoints, system.NForceField, length(sourceTimes))
    scale = FE.force * dt^2 / (FE.rho * FE.dx)        # force density × dt²/ρ, scaled coordinates
    for i in eachindex(weights)
        weights[i] == 0.0 && continue
        source[i, 1, :] .= weights[i] .* scale .* ricker.(sourceTimes)
    end
    stride = max(1, round(Int, 0.25 / dt))
    @time propagation = propagateLinearSystem(system, steps, dt; sourceFull=source,
        output_stride=stride, blowup_limit=1e8, solver_name="OPT3 cut-cell 1-D η=$η", scheme=:direct)
    propagation.stopped_early && error("propagation stopped early (blow-up)")
    times = propagation.times
    history = propagation.history                 # (nx, NField, ntimes)

    exact(xr, t) = FE.force / (2 * FE.rho * FE.c) * sum(weights[j] * (
        rickerIntegral(t - abs(xr - x[j]) / FE.c) + rickerIntegral(t - abs(xr - (2xs - x[j])) / FE.c))
        for j in eachindex(x) if weights[j] != 0)
    exactDirect(xr, t) = FE.force / (2 * FE.rho * FE.c) * sum(weights[j] *
        rickerIntegral(t - abs(xr - x[j]) / FE.c) for j in eachindex(x) if weights[j] != 0)

    receivers = [x[argmin(abs.(x .- (xs - d)))] for d in FE.receiverDistances]
    numericAt(trace, t) = begin
        k = clamp(searchsortedlast(times, t), 1, length(times) - 1)
        w = clamp((t - times[k]) / (times[k+1] - times[k]), 0.0, 1.0)
        (1 - w) * trace[k] + w * trace[k+1]
    end
    # amplitude convention (source normalisation): fitted once, on the direct
    # wave at the farthest receiver, where it is well separated from the reflection
    farthest = argmin(receivers)
    iFar = argmin(abs.(x .- receivers[farthest]))
    tDirectFar = t0 + abs(receivers[farthest] - xsrc) / FE.c
    directWindow = findall(t -> abs(t - tDirectFar) <= 1.0 / f0, times)
    directFar = [exactDirect(receivers[farthest], t) for t in times[directWindow]]
    amplitude = sum(Float64.(history[iFar, 1, directWindow]) .* directFar) / sum(abs2, directFar)
    @printf("  source amplitude factor (numeric / exact, direct wave): %.5f\n", amplitude)

    caseResult = Dict{String,Any}("eta" => η, "selfCheck" => selfCheck, "receivers" => receivers,
        "amplitude" => amplitude)
    traces = []
    for xr in receivers
        i = argmin(abs.(x .- xr))
        numericTrace = Float64.(history[i, 1, :]) ./ amplitude
        exactTrace = [exact(xr, t) for t in times]
        relError = sqrt(sum(abs2, numericTrace .- exactTrace) / sum(abs2, exactTrace))
        # reflection: numeric minus the exact direct wave, against the exact reflected wave
        tReflect = t0 + (xs - xsrc + xs - xr) / FE.c
        tw = range(tReflect - 0.8 / f0, tReflect + 0.8 / f0; length=401)
        exactReflected = [exact(xr, t) - exactDirect(xr, t) for t in tw]
        lags = -2.0:0.01:2.0
        scores = [sum((numericAt(numericTrace, t + lag) - exactDirect(xr, t + lag)) * e
                      for (t, e) in zip(tw, exactReflected)) for lag in lags]
        bestLag = lags[argmax(scores)]
        numericReflected = [numericAt(numericTrace, t) - exactDirect(xr, t) for t in tw]
        peakRatio = maximum(abs, numericReflected) / max(maximum(abs, exactReflected), eps())
        @printf("  receiver x = %6.1f  rel. L2 error %.3e  reflection lag %+.3f (stair-step: %+.3f)  reflected peak ratio %.4f\n",
            xr, relError, bestLag, -2η * FE.dx / FE.c, peakRatio)
        push!(traces, (x=xr, times=times, numeric=numericTrace, exact=exactTrace, relError=relError,
            lag=bestLag, peakRatio=peakRatio))
    end
    caseResult["traces"] = traces
    resultsFE[case] = caseResult

    fig = Figure(size=(1100, 260 * length(traces) + 60))
    Label(fig[0, 1], "1-D free end at x_s = $(xs) (η = $η Δx beyond the last solid node), OPT3 cut-cell vs exact"; fontsize=14)
    for (k, tr) in enumerate(traces)
        ax = Axis(fig[k, 1]; ylabel="x = $(tr.x)", xlabel=k == length(traces) ? "time" : "")
        lines!(ax, tr.times, tr.exact; color=:black, linewidth=2, label="exact")
        lines!(ax, tr.times, tr.numeric; color=:crimson, linestyle=:dash, label="OPT3 cut-cell")
        k == 1 && axislegend(ax; position=:lt)
    end
    devSave(fig, "freeEnd1D_$(case)")
end

println("\n== summary: reflection timing error vs stair-step (time units, Δx = c = 1)")
@printf("  %-8s %-10s %-14s %-14s %-12s\n", "case", "receiver", "cut-cell lag", "stair-step", "peak ratio")
for case in DEV.cases
    r = resultsFE[case]
    for tr in r["traces"]
        @printf("  %-8s %-10.1f %+-14.3f %+-14.3f %-12.4f\n", case, tr.x, tr.lag, -2r["eta"], tr.peakRatio)
    end
end
