# dev/new_Paul_models/stages/neutrinos.jl -- 4 surface detectors × 2 sides × (PREM, sb11, sb21)
@isdefined(seismoModels) || include(joinpath(PAUL_STAGES, "synthetics.jl"))   # models (cache hits)
include(joinpath(PAUL_STAGES, "neutrino_env.jl"))
DEV_STEP[] = "neutrinos"
icbRadius = 1_221_500.0

# ---- nₑ fields (PREM is the spherical average of both models)
YePREM = nuElectronFraction(seismoDiskGrid, cmbRadius; yeCore = NU.yeCore, yeMantle = NU.yeMantle)
neFields = Dict{String,Matrix{Float64}}("PREM" => electronDensity(prem.ρ, YePREM))
for entry in SEISMO_MODELS
    name = "$(entry.model)_$(entry.snap)"
    files = listStagYYSnapshots(joinpath(PAUL_MODELS_DIR, entry.model, "Output"))[entry.snap]
    δlnYe = yeLateralAnomaly(files["bs"], seismoDiskGrid, planetRadius, cmbRadius;
        Ye_hz = PAUL_YE_HZ, Ye_morb = PAUL_YE_MORB, seamBlendDeg = PAUL_SEAM_BLEND_DEG)
    Ye = nuElectronFraction(seismoDiskGrid, cmbRadius; yeCore = NU.yeCore, yeMantle = NU.yeMantle, δlnYe = δlnYe)
    neFields[name] = electronDensity(seismoModels[name].ρ, Ye)
    eY = extrema(δlnYe[seismoDiskGrid.mantle])
    eN = extrema((neFields[name] ./ max.(neFields["PREM"], eps()) .- 1)[seismoDiskGrid.mantle])
    @printf("%s: δlnYₑ %+.4f..%+.4f, δln nₑ %+.3f..%+.3f (mantle)\n", name, eY[1], eY[2], eN[1], eN[2])
end

# ---- flux, cross sections, oscillation parameters (cached Daemonflux table, same binning as the notebook)
set_default_neutrino_flux!(:Daemonflux; flux_mode = :total, return_uncertainties = false)
nuFluxTable = produce_or_load_neutrino_flux(NU.binning)
nuXsecTable = evaluate_neutrino_cross_sections(nuFluxTable.energies)
nuOsc = set_oscillation_parameters()
NU_SIDES = (-1, 1)

nuResults = Dict{Tuple{String,Symbol,Int},Any}()
for name in SEISMO_MODEL_NAMES, det in SEISMO_POLES, side in NU_SIDES
    params = Dict{String,Any}("stage" => "nu", "codeVersion" => CODE_VERSION["nu"], "model" => name,
        "detector" => string(det), "angle" => sourceAngle(det), "side" => side,
        "dx" => S.dx, "mohoRadius" => S.mohoRadius, "clip" => S.clip, "seamBlendDeg" => PAUL_SEAM_BLEND_DEG,
        "yeCore" => NU.yeCore, "yeMantle" => NU.yeMantle, "YeHz" => PAUL_YE_HZ, "YeMorb" => PAUL_YE_MORB,
        "energies" => collect(NU.binning.energies), "cosθ" => collect(NU.binning.cosθgrid),
        "Δbaseline" => NU.Δbaseline, "detectorDepth" => NU.detectorDepth, "massKg" => NU.massKg, "years" => NU.years)
    prefix = "nu_$(name)_$(det)_$(side > 0 ? "p" : "m")"
    devRedo!("nu", params, SEISMO_CACHE_DIR, prefix)
    nuResults[(name, det, side)] = myProduceOrLoad(params, SEISMO_CACHE_DIR, prefix) do cfg
        @info "neutrinos" name det side
        runNeutrinoDetector(neFields[name], seismoDiskGrid, planetRadius, sourceAngle(det), side, NU,
            nuOsc, nuFluxTable, nuXsecTable)
    end
end

nuRows = neutrinoTable(nuResults, SEISMO_MODEL_NAMES, SEISMO_POLES, NU_SIDES, NU, planetRadius, cmbRadius, icbRadius)
let dir = joinpath(DEV_OUTPUT_DIR[], "neutrinos")
    mkpath(dir)
    open(joinpath(dir, "neutrino_table.csv"), "w") do io
        println(io, "detector,side,channel,class,band_GeV,N_PREM,N_sb11,N_sb21,dNN_sb11,dNN_sb21,chi2_sb11_vs_sb21,chi2_sb11_vs_PREM,chi2_sb21_vs_PREM")
        for r in nuRows
            @printf(io, "%s,%s,%s,%s,%g-%g,%.4g,%.4g,%.4g,%.5f,%.5f,%.4g,%.4g,%.4g\n", r.detector, r.side > 0 ? "+" : "-",
                r.channel, r.class, r.band[1], r.band[2], r.N_PREM, r.N_1, r.N_2, r.dNN_1, r.dNN_2, r.chi2_12, r.chi2_1P, r.chi2_2P)
        end
    end
    println("neutrino table -> ", joinpath(dir, "neutrino_table.csv"))
end

paulSave(neutrinoMapFigure(neFields, seismoDiskGrid, nuResults, SEISMO_MODEL_NAMES, SEISMO_POLES, NU_SIDES,
    planetRadius, cmbRadius), "nu_00_maps")
for ch in ("νe", "νμ")
    tag = ch == "νe" ? "e" : "mu"
    paulSave(neutrinoBarFigure(nuRows, ch, SEISMO_MODEL_NAMES, SEISMO_POLES, NU), "nu_01_bars_$(tag)")
    for det in SEISMO_POLES
        paulSave(neutrinoDetectorFigure(nuResults, SEISMO_MODEL_NAMES, det, NU_SIDES, ch), "nu_02_$(det)_$(tag)")
    end
end
