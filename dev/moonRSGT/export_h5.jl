# dev/moonRSGT/export_h5.jl -- catalogue (JLD2) -> compact HDF5 ("moonRSGT-h5 v2") for the students
#
# Run on a compute node (srun / sbatch), never on stella*:
#   julia --project=. dev/moonRSGT/export_h5.jl --test=00001               # size vs error table, choose keepbits
#   julia --project=. dev/moonRSGT/export_h5.jl --task=3 --per-task=10     # Slurm array (slurm/export_h5.sbatch)
#   julia --project=. dev/moonRSGT/export_h5.jl --models=1:5 --keepbits=12
# --in DIR          catalogue folder (default ENV MOON_RSGT_DIR, else /gpfs/scratch/nobuaki/moonRSGT/catalogue)
# --out DIR         export folder   (default ENV MOON_RSGT_H5_DIR, else /gpfs/scratch/nobuaki/moonRSGT/catalogue_h5)
# --keepbits=12     float32 mantissa bits kept (23 = lossless); relative error ≤ 2^-(keepbits+1)
# --compressor=deflate|zstd|none   deflate (+shuffle) needs nothing in Python; zstd needs `import hdf5plugin`
# --theta-block=16  distances per stored block (a reader decompresses only the blocks it needs)
# --with-model      also copy the .nd table (CONFIDENTIAL, off by default)
# --force           overwrite existing files
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

if haskey(o, "test")
    # ------------------------------------------------------------------ size vs error
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
    ref = [moonSynthetics(catl; station=(c.stla, c.stlo), source=(c.evla, c.evlo, c.evdp),
                          M=fromGCMT(M.Mrr, M.Mtt, M.Mpp, M.Mrt, M.Mrp, M.Mtp), comps=(:Z, :R, :T),
                          quantity=:velocity).data for c in cases]
    @printf("model %s: plain Float32 = %.0f MB\n", tid, raw / 1e6)
    @printf("%-10s %-9s %10s %8s %14s %14s\n", "keepbits", "compr.", "size (MB)", "ratio", "max trace err", "write (s)")
    for kb in (23, 16, 12, 10, 8), comp in (compressor,)
        f = joinpath(tdir, "rsgt_$(tid).h5")
        t0 = time()
        sz = exportMoonH5(catl, f; keepbits=kb, compressor=comp, θblock=θblock)
        tw = time() - t0
        _MOON_H5_LAST[] = nothing
        err = 0.0
        for (k, c) in enumerate(cases)
            _, tr = moonH5Synthetics(f; c..., M..., comps="ZRT", quantity="velocity")
            for comp_ in ("Z", "R", "T")
                a = ref[k][Symbol(comp_)]
                err = max(err, norm(tr[comp_] .- a) / max(norm(a), 1e-300))
            end
        end
        @printf("%-10d %-9s %10.1f %8.2f %14.2e %14.1f\n", kb, string(comp), sz / 1e6, raw / sz, err, tw)
        rm(f; force=true)
    end
    println("(trace error: relative L2 of Z, R, T velocity vs the JLD2 catalogue, 4 sources incl. 9 and 40 km)")
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
