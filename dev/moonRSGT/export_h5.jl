# dev/moonRSGT/export_h5.jl -- catalogue (JLD2, Julia NamedTuples) -> plain HDF5 for Python (h5py) users
#
#   julia --project=. dev/moonRSGT/export_h5.jl --task=3 --per-task=10     # Slurm array (slurm/export_h5.sbatch)
#   julia --project=. dev/moonRSGT/export_h5.jl --models=1:5
# --in  catalogue folder (default ENV MOON_RSGT_DIR, else /gpfs/scratch/nobuaki/moonRSGT/catalogue)
# --out export folder   (default ENV MOON_RSGT_H5_DIR, else /gpfs/scratch/nobuaki/moonRSGT/catalogue_h5)
# --with-model          also copy the .nd table (depth, vp, vs, rho): CONFIDENTIAL, off by default
#
# Layout of rsgt_<id>.h5 (shapes as seen from h5py / numpy, C order; indices 1-based in pairs_*):
#   strain_re, strain_im  float32 (nω, 10, nr, nθ)   RSGT strain pairs, pole frame, φ = 0
#   disp_re,   disp_im    float32 (nω,  5, nr, nθ)   displacement pairs
#   pairs_strain_src, pairs_strain_comp   int (10)   force component s (r,θ,φ) and strain component c (rr,θθ,φφ,rθ,rφ,θφ)
#   pairs_disp_src,   pairs_disp_comp     int (5)    force component s and displacement component c (r,θ,φ)
#   omega (nω) rad/s, omega_i, tlen (s), r_km (nr, increasing), layer (nr), theta_deg (nθ)
#   R_km, max_depth_km, model_id, format = "moonRSGT-h5 v1"
# Units: strain per unit force (code units, see educationPurposePythonScripts/RSGTreader/moon_rsgt_py.py), displacement in km.

using JLD2, Printf
include(joinpath(@__DIR__, "lib", "moonSynthetics.jl"))     # moonModelIds, loadMoonRSGT

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
withModel = haskey(o, "with-model")
ids = moonModelIds(; dir=inDir)
sel = if haskey(o, "task")
    k = parse(Int, o["task"]); n = parse(Int, get(o, "per-task", "10"))
    ((k-1)*n+1):min(k * n, length(ids))
elseif haskey(o, "models")
    a, b = parse.(Int, split(o["models"], ':')); a:min(b, length(ids))
else
    1:length(ids)
end
mkpath(outDir)
"(nω, np, nr, nθ) Julia array -> written so that h5py sees (nω, np, nr, nθ) in C order."
toC(A) = permutedims(A, (4, 3, 2, 1))

for k in sel
    id = ids[k]
    out = joinpath(outDir, "rsgt_$(id).h5")
    if isfile(out)
        println("skip (exists): ", out); continue
    end
    t0 = time()
    c = loadMoonRSGT(id; dir=inDir)
    d = c.disk
    tmp = out * ".part"
    jldopen(tmp, "w") do f
        f["format"] = "moonRSGT-h5 v1"
        f["model_id"] = string(id)
        f["strain_re"] = toC(Float32.(real.(d.strain))); f["strain_im"] = toC(Float32.(imag.(d.strain)))
        f["disp_re"] = toC(Float32.(real.(d.disp)));     f["disp_im"] = toC(Float32.(imag.(d.disp)))
        f["pairs_strain_src"] = [p[1] for p in d.pairs.strain]; f["pairs_strain_comp"] = [p[2] for p in d.pairs.strain]
        f["pairs_disp_src"] = [p[1] for p in d.pairs.disp];     f["pairs_disp_comp"] = [p[2] for p in d.pairs.disp]
        f["omega"] = collect(Float64, d.ω); f["omega_i"] = Float64(d.ωI); f["tlen"] = Float64(d.tlen)
        f["r_km"] = collect(Float64, d.r_km); f["layer"] = collect(Int, d.seg); f["theta_deg"] = collect(Float64, d.θdeg)
        f["R_km"] = Float64(c.info.R); f["max_depth_km"] = Float64(c.maxDepth)
        if withModel
            f["model_depth"] = c.model.depth; f["model_vp"] = c.model.vp; f["model_vs"] = c.model.vs; f["model_rho"] = c.model.ρ
        end
    end
    mv(tmp, out; force=true)
    @printf("%s -> %s (%.0f s)\n", id, out, time() - t0)
end
