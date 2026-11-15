# dev/new_Paul_models/stages/phases.jl -- phase windows (needs the synthetics: cache hits)
@isdefined(seismoResults) || include(joinpath(PAUL_STAGES, "synthetics.jl"))
DEV_STEP[] = "phases"
seismoΔgrid = 0.0:S.receiverStepDeg:180.0
taupParams = Dict{String,Any}("stage" => "taup", "codeVersion" => CODE_VERSION["taup"],
    "model" => "PREM.psv.model", "depthKm" => SEISMO_DEPTH_KM, "phases" => SEISMO_PHASES,
    "Δstep" => S.receiverStepDeg)
devRedo!("taup", taupParams, SEISMO_CACHE_DIR, "taup")
seismoTaup = myProduceOrLoad(taupParams, SEISMO_CACHE_DIR, "taup") do cfg
    println("PREM τ–p first arrivals, source depth ", SEISMO_DEPTH_KM, " km:")
    Dict{String,Any}("table" => taupTable(SEISMO_RAY_PREM, SEISMO_PHASES, SEISMO_DEPTH_KM, seismoΔgrid))
end
seismoTable = seismoTaup["table"]

seismoMeasurements = measurePhases(seismoResults, SEISMO_MODEL_NAMES, seismoTable, seismoΔgrid,
    seismoReceivers.φ, SEISMO_WINDOW)
let dir = joinpath(DEV_OUTPUT_DIR[], "phases")
    mkpath(dir)
    file = joinpath(dir, "phase_measurements.csv")
    open(file, "w") do io
        println(io, "phase,pole,model,phi_deg,delta_deg,t_PREM_s,dt_s,dlnA,cc,clean")
        for m in seismoMeasurements
            @printf(io, "%s,%s,%s,%.1f,%.1f,%.3f,%.3f,%.4f,%.4f,%s\n", m.phase, m.pole, m.model,
                m.φ, m.Δ, m.tPred, m.δt, m.dlnA, m.cc, m.clean)
        end
    end
    println("measurements -> ", file)
end

for ph in SEISMO_PHASES
    try
        paulSave(phaseFigure(ph, seismoResults, SEISMO_MODEL_NAMES, seismoAnomalies, seismoDiskGrid,
            seismoTable, seismoΔgrid, seismoMeasurements, planetRadius, cmbRadius), "phase_$(ph)")
    catch err
        @warn "phase figure failed" ph exception = (err, catch_backtrace())
    end
end
paulSave(phaseSummaryFigure(seismoMeasurements, SEISMO_PHASES, SEISMO_MODEL_NAMES[2:end]), "phase_00_summary")
