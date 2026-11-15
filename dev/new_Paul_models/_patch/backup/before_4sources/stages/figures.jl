# dev/new_Paul_models/stages/figures.jl -- every field of every snapshot of `paulModel`
paulOutput = joinpath(PAUL_MODELS_DIR, paulModel, "Output")
paulSnaps = listStagYYSnapshots(paulOutput)
println("model ", paulModel, ": snapshots ", join(sort(collect(keys(paulSnaps))), ", "))
paulFailures = String[]
for snap in sort(collect(keys(paulSnaps)))
    files = paulSnaps[snap]
    hRef = stagyyHeader(files[first(sort(collect(keys(files))))])
    rmax = hRef.rcmb + hRef.zc[1, hRef.nz+1]
    paulGrid = diskGrid(rmax; n = PAUL_GRID_N)
    processed = Dict{String,Any}()
    for code in sort(collect(keys(files)))
        code in PAUL_SKIP_CODES && continue
        file = files[code]
        try
            h = stagyyHeader(file)
            if h.nval > 1
                fig = vectorFigure(file, paulModel, snap, paulGrid)
            elseif h.nz <= 2
                fig = surfaceFigure(file, code, paulModel, snap)
            else
                p = processScalar(file, code, paulGrid; seamBlendDeg = PAUL_SEAM_BLEND_DEG)
                processed[code] = p
                fig = scalarFigure(p, code, paulModel, snap, paulGrid)
            end
            paulSave(fig, "$(snap)_$(code)")
        catch err
            push!(paulFailures, "$snap $code: " * sprint(showerror, err))
            @warn "failed" snap code exception = (err, catch_backtrace())
        end
    end
    isempty(processed) || paulSave(overviewFigure(processed, paulModel, snap, paulGrid),
        "$(snap)_00_overview")
    if haskey(files, "rho") && haskey(files, "bs")
        try
            paulSave(electronDensityFigure(files, paulModel, snap, paulGrid;
                    Ye_hz = PAUL_YE_HZ, Ye_morb = PAUL_YE_MORB, bsIsVolume = PAUL_BS_IS_VOLUME,
                    ΔρMorb = PAUL_DRHO_MORB, seamBlendDeg = PAUL_SEAM_BLEND_DEG),
                "$(snap)_01_electron_density")
        catch err
            push!(paulFailures, "$snap electron density: " * sprint(showerror, err))
            @warn "failed" snap exception = (err, catch_backtrace())
        end
    end
end
isempty(paulFailures) ? println("all figures done for ", paulModel) :
    (println("FAILED:"); foreach(println, paulFailures))
