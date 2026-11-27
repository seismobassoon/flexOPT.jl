# dev/moonRSGT/export_h5.jl -- catalogue (JLD2) -> compact HDF5 ("moonRSGT-h5 v2") for the students
#
# Run on a compute node (srun / sbatch), never on stella*:
#   julia --project=. dev/moonRSGT/export_h5.jl --test=00001               # size vs error table, choose keepbits
#   julia --project=. dev/moonRSGT/export_h5.jl --test=00001 --exact --plot   # + exact solution, waveform figure
#   julia --project=. dev/moonRSGT/export_h5.jl --task=3 --per-task=10     # Slurm array (slurm/export_h5.sbatch)
#   julia --project=. dev/moonRSGT/export_h5.jl --models=1:5 --keepbits=12
# --in DIR          catalogue folder (default ENV MOON_RSGT_DIR, else /gpfs/scratch/nobuaki/moonRSGT/catalogue)
# --out DIR         export folder   (default ENV MOON_RSGT_H5_DIR, else /gpfs/scratch/nobuaki/moonRSGT/catalogue_h5)
# --keepbits=12     float32 mantissa bits kept (23 = lossless); relative error ≤ 2^-(keepbits+1)
# --compressor=deflate|zstd|none   deflate (+shuffle) needs nothing in Python; zstd needs `import hdf5plugin`
# --theta-block=16  distances per stored block (a reader decompresses only the blocks it needs)
# --with-model      also copy the .nd table (CONFIDENTIAL, off by default)
# --force           overwrite existing files
# with --test:  --plot   waveform figure (wiggles + residuals) in dev/moonRSGT/outputs/export_test/
#               --exact  also compare with the EXACT source-side solution (solver run at the true source
#                        position, no interpolation, as run.jl check; needs the .nd file in --nd-dir)
# --nd-dir DIR      .nd models (default ENV MOON_ND_DIR, else /gpfs/scratch/nobuaki/moonRSGT/OUTPUT_ND)
# Layout: see lib/moonExportH5.jl.  The same step can run right after each catalogue file is built
# (run.jl build --h5).

using JLD2, Printf, LinearAlgebra
include(joinpath(@__DIR__, "lib", "moonSynthetics.jl"))     # moonModelIds, loadMoonRSGT, moonSynthetics
include(joinpath(@__DIR__, "lib", "moonExportH5.jl"))

function exportArgs(args)
    o = Dict{String,String}()
    for a in args
        startswith(a, "--") || error("unexpected argument $a")
        kv = split(a[3:end], '='; limit=2)
        o[kv[1]] = length(kv) == 2 ? kv[2] : "true"
    end
    return o
end

o = exportArgs(ARGS)
inDir = get(o, "in", moonRSGTDir())
outDir = get(o, "out", get(ENV, "MOON_RSGT_H5_DIR", "/gpfs/scratch/nobuaki/moonRSGT/catalogue_h5"))
keepbits = parse(Int, get(o, "keepbits", "12"))
compressor = Symbol(get(o, "compressor", "deflate"))
θblock = parse(Int, get(o, "theta-block", "16"))
withModel = haskey(o, "with-model")
ids = moonModelIds(; dir=inDir)
ndDir = get(o, "nd-dir", get(ENV, "MOON_ND_DIR", "/gpfs/scratch/nobuaki/moonRSGT/OUTPUT_ND"))
EXPORT_EXACT = haskey(o, "test") && haskey(o, "exact")
EXPORT_PLOT = haskey(o, "test") && haskey(o, "plot")
EXPORT_EXACT && include(joinpath(@__DIR__, "..", "common", "devBootstrap.jl"))
EXPORT_EXACT && include(joinpath(@__DIR__, "env.jl"))            # solver (CPU), lib/*, moonRecipes, shConfigs
EXPORT_PLOT && @eval using CairoMakie

if haskey(o, "test")
    # ------------------------------------------------------------------ size vs error (+ waveforms)
    include(joinpath(@__DIR__, "..", "..", "educationPurposePythonScripts", "RSGTreader", "moon_rsgt_h5.jl"))
    tid = o["test"]
    catl = loadMoonRSGT(tid; dir=inDir)
    raw = 8 * (length(catl.disk.strain) + length(catl.disk.disp))          # plain Float32 re + im
    tdir = mktempdir(get(ENV, "TMPDIR", tempdir()))
    cases = [(stla=-3.0, stlo=-23.4, evla=20.0, evlo=15.0, evdp=300.0),
             (stla=-3.0, stlo=-23.4, evla=-25.0, evlo=70.0, evdp=913.0),
             (stla=-3.0, stlo=-23.4, evla=10.0, evlo=-60.0, evdp=9.0),
             (stla=-3.0, stlo=-23.4, evla=-1.0, evlo=-21.0, evdp=40.0)]
    M = (Mrr=3e14, Mtt=-1e15, Mpp=7e14, Mrt=2e14, Mrp=-4e14, Mtp=5e14)
    Mcat = fromGCMT(M.Mrr, M.Mtt, M.Mpp, M.Mrt, M.Mrp, M.Mtp)
    comps = ("Z", "R", "T")
    tref = Float64[]
    ref = map(cases) do c                                                    # JLD2 catalogue (uncompressed)
        s_ = moonSynthetics(catl; station=(c.stla, c.stlo), source=(c.evla, c.evlo, c.evdp), M=Mcat,
                            comps=(:Z, :R, :T), quantity=:velocity)
        global tref = s_.t
        Dict(k => s_.data[Symbol(k)] for k in comps)
    end
    exact = nothing
    if EXPORT_EXACT                                                          # source-side solver, exact position
        nd = joinpath(ndDir, string(tid) * ".nd")
        mm = moonModelFromND(nd; fitTol=catl.params["fitTol"])
        psv, sh = moonSystems(mm.model, catl.params; recipes=moonRecipes, configs=shConfigs)
        R = mm.info.R; Mcode = collect(Float64, Mcat) ./ 1e18
        exact = map(cases) do c
            Rcv = pointFromLatLon(R, 0.0, c.stla, c.stlo); Q = pointFromLatLon(R, c.evdp, c.evla, c.evlo)
            Rs, _ = sourceToPole(:moment, Mcode, Q); θp, _, _ = frameOf(Rs, Rcv)
            dS = computeSGTDisk(psv, sh, mm.model; family=:moment, rs_km=R - c.evdp,
                                radii=(r_km=[R], side=[:upper]), θdeg=[θp], tlen=catl.params["tlen"],
                                fmax=catl.params["fmax"], omegaiTlen=catl.params["omegaiTlen"],
                                p=sourceOrderOf(catl.params["scheme"]), lmaxOf=moonLmaxOf(catl.params, R),
                                lcap=5000, verbose=false)
            Dict(k => moonTimeSeries(dS.ω, dS.ωI, dS.tlen,
                                     [1e3 * dot(collect(receiverDirection(Symbol(k), Q, Rcv)),
                                                sgtAt(dS, i, Q, Mcode, Rcv; which=:disp)) for i in eachindex(dS.ω)];
                                     halfDuration=10.0, quantity=:velocity)[2] for k in comps)
        end
    end
    relerr(a, b) = norm(a .- b) / max(norm(b), 1e-300)
    @printf("model %s: plain Float32 = %.0f MB\n", tid, raw / 1e6)
    @printf("%-9s %-8s %10s %7s %16s %16s\n", "keepbits", "compr.", "size (MB)", "ratio",
            "err vs JLD2", EXPORT_EXACT ? "err vs exact" : "")
    if EXPORT_EXACT
        e0 = maximum(relerr(ref[k][c], exact[k][c]) for k in eachindex(cases), c in comps)
        @printf("%-9s %-8s %10s %7s %16s %16.2e\n", "JLD2", "-", "-", "-", "-", e0)
    end
    kept = Dict{Int,Any}()
    for kb in (23, 16, 12, 10, 8)
        f = joinpath(tdir, "rsgt_$(tid).h5")
        sz = exportMoonH5(catl, f; keepbits=kb, compressor=compressor, θblock=θblock)
        _MOON_H5_LAST[] = nothing
        tr = [moonH5Synthetics(f; c..., M..., comps="ZRT", quantity="velocity")[2] for c in cases]
        kept[kb] = tr
        ec = maximum(relerr(tr[k][c], ref[k][c]) for k in eachindex(cases), c in comps)
        ex = EXPORT_EXACT ? maximum(relerr(tr[k][c], exact[k][c]) for k in eachindex(cases), c in comps) : NaN
        @printf("%-9d %-8s %10.1f %7.2f %16.2e %16s\n", kb, string(compressor), sz / 1e6, raw / sz, ec,
                EXPORT_EXACT ? @sprintf("%.2e", ex) : "")
        rm(f; force=true)
    end
    println("(max relative L2 misfit of Z, R, T velocity over 4 sources: 300, 913, 9 and 40 km; hd 10 s)")
    if EXPORT_PLOT
        base = EXPORT_EXACT ? exact : ref
        bname = EXPORT_EXACT ? "exact (source side)" : "JLD2 catalogue"
        showkb = [(12, :red), (8, :orange)]
        fig = Main.CairoMakie.Figure(size=(1500, 300 * length(cases)))
        for (k, c) in enumerate(cases), (j, comp) in enumerate(comps)
            ax = Main.CairoMakie.Axis(fig[k, j]; xlabel=k == length(cases) ? "time (s)" : "", ylabel=j == 1 ? "velocity (m/s)" : "",
                title=@sprintf("%s, evdp %.0f km, Δ = %.1f°", comp, c.evdp, moonH5Distance(c.evla, c.evlo, c.stla, c.stlo)))
            b = base[k][comp]; amp = maximum(abs, b)
            Main.CairoMakie.lines!(ax, tref, b; color=:black, linewidth=2, label=bname)
            EXPORT_EXACT && Main.CairoMakie.lines!(ax, tref, ref[k][comp]; color=:dodgerblue, linewidth=1, label="JLD2 catalogue")
            for (kb, col) in showkb
                Main.CairoMakie.lines!(ax, tref, kept[kb][k][comp]; color=col, linewidth=0.8, label="compressed, $kb bits")
                # residual, offset below and magnified x10
                Main.CairoMakie.lines!(ax, tref, -1.6amp .+ 10 .* (kept[kb][k][comp] .- b); color=col, linewidth=0.6)
            end
            EXPORT_EXACT && Main.CairoMakie.lines!(ax, tref, -1.6amp .+ 10 .* (ref[k][comp] .- b); color=:dodgerblue, linewidth=0.6)
            Main.CairoMakie.text!(ax, tref[1], -1.25amp; text="residual × 10", fontsize=10)
            (k == 1 && j == 1) && Main.CairoMakie.axislegend(ax; position=:rt, labelsize=9)
        end
        od = joinpath(@__DIR__, "outputs", "export_test"); mkpath(od)
        fp = joinpath(od, "export_test_$(tid)" * (EXPORT_EXACT ? "_exact" : "") * ".png")
        Main.CairoMakie.save(fp, fig); println("figure -> ", fp)
    end
else
    # ------------------------------------------------------------------ export
    sel = if haskey(o, "task")
        k = parse(Int, o["task"]); n = parse(Int, get(o, "per-task", "10"))
        ((k-1)*n+1):min(k * n, length(ids))
    elseif haskey(o, "models")
        m1, m2 = parse.(Int, split(o["models"], ':')); m1:min(m2, length(ids))
    else
        1:length(ids)
    end
    mkpath(outDir)
    for k in sel
        mid = ids[k]
        out = joinpath(outDir, "rsgt_$(mid).h5")
        if isfile(out) && !haskey(o, "force")
            println("skip (exists): ", out); continue
        end
        t0 = time()
        sz = exportMoonH5(loadMoonRSGT(mid; dir=inDir), out; keepbits=keepbits, compressor=compressor,
                          θblock=θblock, withModel=withModel)
        @printf("%s -> %s  %.0f MB (keepbits %d, %s)  %.0f s\n", mid, out, sz / 1e6, keepbits, compressor, time() - t0)
    end
end
