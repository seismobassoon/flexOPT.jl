# dev/new_Paul_models/stages/features.jl -- which feature separates sb11 from sb21?
@isdefined(seismoResults) || include(joinpath(PAUL_STAGES, "synthetics.jl"))
DEV_STEP[] = "features"
# PREM phase times: same cached τ–p table as the `phases` step
seismoΔgrid = 0.0:S.receiverStepDeg:180.0
taupParams = Dict{String,Any}("stage" => "taup", "codeVersion" => CODE_VERSION["taup"],
    "model" => "PREM.psv.model", "depthKm" => SEISMO_DEPTH_KM, "phases" => SEISMO_PHASES,
    "Δstep" => S.receiverStepDeg)
seismoTaup = myProduceOrLoad(taupParams, SEISMO_CACHE_DIR, "taup") do cfg
    Dict{String,Any}("table" => taupTable(SEISMO_RAY_PREM, SEISMO_PHASES, SEISMO_DEPTH_KM, seismoΔgrid))
end
seismoTable = seismoTaup["table"]

featureCfg = merge(FEATURES, (duration = S.duration, period = S.period))
featureWins = featureWindows(seismoTable, seismoΔgrid, featureCfg)
featureModels = SEISMO_MODEL_NAMES[2:end]
length(featureModels) == 2 || error("the feature step compares exactly two models")
println("feature windows: ", join([w.label for w in featureWins], ", "))
@time "features" featureValues = computeFeatures(seismoResults, featureModels, seismoReceivers.φ,
    featureWins, featureCfg)
featureRows = scoreFeatures(featureValues, featureModels, featureWins, featureCfg)

let dir = joinpath(DEV_OUTPUT_DIR[], "features")
    mkpath(dir)
    open(joinpath(dir, "feature_scores.csv"), "w") do io
        println(io, "rank,source,component,window,kind,feature,n_stations,n_detectable,detectable_fraction,rms_z,median_abs_z,rms_diff_noise_free")
        for (i, r) in enumerate(featureRows)
            @printf(io, "%d,%s,%s,%s,%s,%s,%d,%d,%.4f,%.4f,%.4f,%.6g\n", i, r.source, compLabel(r.comp), r.window,
                r.kind, r.feature, r.nStations, r.nDetectable, r.detectable, r.rmsZ, r.medianAbsZ, r.rmsDiffNoiseFree)
        end
    end
    jldsave(joinpath(dir, "feature_values.jld2");
        values = Dict(join(string.(k), "/") => v for (k, v) in featureValues),
        featureNames = collect(string.(FEATURE_NAMES)), phi = seismoReceivers.φ,
        README = "values[source/comp/model/window]: (station, feature, 1 + nReal); slice 1 noise-free")
    println("feature scores -> ", joinpath(dir, "feature_scores.csv"))
end
println("\nTop 15 (stations with |z| > $(featureCfg.zThreshold)):")
for (i, r) in enumerate(featureRows[1:min(15, end)])
    @printf("%2d  %-6s %-10s %-12s %-15s %3d / %3d stations  rms z %.1f\n", i, r.source, compLabel(r.comp),
        r.window, r.feature, r.nDetectable, r.nStations, r.rmsZ)
end

for pole in SEISMO_POLES
    paulSave(featureScoreFigure(featureRows, featureWins, pole, featureCfg), "features_scores_$(pole)")
end
paulSave(featureRankingFigure(featureRows), "features_00_ranking")
paulSave(topFeaturesFigure(featureRows, featureValues, featureModels, seismoReceivers.φ), "features_01_top")
