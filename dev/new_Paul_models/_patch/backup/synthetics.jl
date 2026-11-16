# dev/new_Paul_models/stages/synthetics.jl
# Models (PREM, PREM×(1+δ) for each SEISMO_MODELS entry), sources SEISMO_POLES,
# FD3 runs cached in data/dev/new_Paul_models/fd3_<model>_<pole>_<hash>.jld2.
include(joinpath(PAUL_STAGES, "seismo_env.jl"))
S = SEISMO
planetRadius = planet1D.my1DDSMmodel.averagedPlanetRadiusInKilometer * 1e3
cmbRadius = planet1D.my1DDSMmodel.averagedPlanetCMBInKilometer * 1e3
println("planet1D model: ", planet1D.DEFAULT_1D_MODEL[], "  R = ", planetRadius / 1e3, " km, CMB = ", cmbRadius / 1e3, " km")
seismoDiskGrid = seismoDisk(S.dx, planetRadius, cmbRadius)
prem = premOnDisk(seismoDiskGrid, planetRadius, S.mohoRadius)
@show size(seismoDiskGrid.X) count(seismoDiskGrid.material)

# ---- models
SEISMO_MODEL_NAMES = ["PREM"]
seismoModels = Dict{String,Any}("PREM" => prem)
seismoAnomalies = Dict{String,Any}()
for entry in SEISMO_MODELS
    files = listStagYYSnapshots(joinpath(PAUL_MODELS_DIR, entry.model, "Output"))[entry.snap]
    an(code) = stagyyLateralAnomaly(files[code], seismoDiskGrid, planetRadius, cmbRadius;
        clip = S.clip, seamBlendDeg = PAUL_SEAM_BLEND_DEG)
    (δvp, exvp), (δvs, exvs), (δρ, exρ) = an("Pvel"), an("Svel"), an("rho")
    name = "$(entry.model)_$(entry.snap)"
    @printf("%s: raw δlnVp %+.3f..%+.3f, δlnVs %+.3f..%+.3f, δlnρ %+.3f..%+.3f (clip ±%.2f)\n",
        name, exvp[1], exvp[2], exvs[1], exvs[2], exρ[1], exρ[2], S.clip)
    seismoAnomalies[name] = (δvp = δvp, δvs = δvs, δρ = δρ)
    seismoModels[name] = (ρ = prem.ρ .* (1 .+ δρ), vp = prem.vp .* (1 .+ δvp), vs = prem.vs .* (1 .+ δvs))
    push!(SEISMO_MODEL_NAMES, name)
end

# ---- one Δt for every run (identical time axes): PREM vp max × (1 + clip)
vpMax = 1e3 * maximum(prem.vp) * (1 + S.clip)
seismoDt = S.cfl / (vpMax * sqrt(2 / S.dx^2))
Mloc = momentTensorLocal(S.M0, S.dip, S.slip)
seismoSources = Dict(pole => poleSource(seismoDiskGrid, pole, S.sourceDepth, planetRadius, Mloc)
    for pole in SEISMO_POLES)
seismoReceivers = ringReceivers(seismoDiskGrid, planetRadius, S.receiverStepDeg, S.receiverDepth, S.dx)
@printf("dt = %.4f s, %d steps per run, %d stations, Mloc = %s\n", seismoDt,
    ceil(Int, S.duration / seismoDt), length(seismoReceivers.φ), string(Mloc))
runCfg = (dx = S.dx, dt = seismoDt, duration = S.duration, period = S.period,
    outputSampling = S.outputSampling, surfaceDamping = S.surfaceDamping, snapshotTimes = S.snapshotTimes)

# ---- runs (cached)
seismoResults = Dict{Tuple{String,Symbol},Any}()
for name in SEISMO_MODEL_NAMES, pole in SEISMO_POLES
    params = Dict{String,Any}(
        "stage" => "fd3", "codeVersion" => CODE_VERSION["fd3"],
        "model" => name, "pole" => string(pole), "planet1D" => string(planet1D.DEFAULT_1D_MODEL[]),
        "dx" => S.dx, "dt" => seismoDt, "duration" => S.duration, "period" => S.period,
        "cfl" => S.cfl, "surfaceDamping" => S.surfaceDamping, "outputSampling" => S.outputSampling,
        "snapshotTimes" => collect(S.snapshotTimes), "mohoRadius" => S.mohoRadius,
        "clip" => S.clip, "seamBlendDeg" => PAUL_SEAM_BLEND_DEG,
        "sourceDepth" => S.sourceDepth, "M0" => S.M0, "dip" => S.dip, "slip" => string(S.slip),
        "receiverStepDeg" => S.receiverStepDeg, "receiverDepth" => S.receiverDepth,
        "dataDir" => name == "PREM" ? "" : basename(PAUL_MODELS_DIR))
    prefix = "fd3_$(name)_$(pole)"
    devRedo!("fd3", params, SEISMO_CACHE_DIR, prefix)
    println("---- ", prefix, "  -> ", devCacheFile(params, SEISMO_CACHE_DIR, prefix))
    seismoResults[(name, pole)] = myProduceOrLoad(params, SEISMO_CACHE_DIR, prefix) do cfg
        runSeismoFD3(seismoModels[name], seismoDiskGrid, seismoSources[pole], seismoReceivers, runCfg;
            label = "$name $pole")
    end
    # runs cached before the 2026-11-15 fix have an uninitialised first sample
    # (uz = similar(ux)); u(t = 0) = 0 exactly, so set it instead of recomputing.
    for key in ("uVert", "uHoriz")
        seismoResults[(name, pole)][key][1, :] .= 0
    end
end

# ---- all synthetics in one file (outputs/synthetics/synthetics.jld2)
let dir = joinpath(DEV_OUTPUT_DIR[], "synthetics")
    mkpath(dir)
    file = joinpath(dir, "synthetics.jld2")
    jldopen(file, "w") do f
        f["README"] = "u (m) at stations every $(S.receiverStepDeg)°; phi: 0 = north pole, 90 = +x, " *
            "180 = south pole; uVert = u·e_r (up), uHoriz = u·e_phi (phi increasing); " *
            "rows = time samples, columns = stations; moment $(S.M0) N m/m, $(S.dip)° $(S.slip), " *
            "depth $(S.sourceDepth / 1e3) km, smooth step, period $(S.period) s"
        f["phi"] = seismoReceivers.φ
        f["models"] = SEISMO_MODEL_NAMES
        for name in SEISMO_MODEL_NAMES, pole in SEISMO_POLES
            r = seismoResults[(name, pole)]
            f["$(name)/$(pole)/time"] = r["time"]
            f["$(name)/$(pole)/uVert"] = r["uVert"]
            f["$(name)/$(pole)/uHoriz"] = r["uHoriz"]
        end
    end
    println("synthetics -> ", file)
end

# ---- figures (PREM travel-time curves from dev/raytheory, same PREM polynomials)
SEISMO_RAY_PREM = readDSMModel(joinpath(flexopt_root, "dataInput", "model", "earth", "PREM.psv.model"))
SEISMO_DEPTH_KM = S.sourceDepth / 1e3
DEV_STEP[] = "synthetics"
paulSave(seismoModelFigure(seismoAnomalies, seismoDiskGrid, seismoSources, planetRadius, cmbRadius), "models")
for pole in SEISMO_POLES, comp in (:vert, :horiz)
    paulSave(anrSectionFigure(seismoResults, SEISMO_MODEL_NAMES, pole, comp), "section_$(pole)_$(comp)")
    paulSave(recordSectionFigure(seismoResults, SEISMO_MODEL_NAMES, pole, comp), "section_overlay_$(pole)_$(comp)")
end
for pole in SEISMO_POLES
    paulSave(snapshotFigure(seismoResults, SEISMO_MODEL_NAMES, pole, seismoDiskGrid), "snapshots_$(pole)")
end
