# dev/moonRSGT/synth.jl -- command-line synthetics from the Moon RSGT catalogue
# (used by educationPurposePythonScripts/RSGTreader/moon_rsgt.py; also handy on its own).  Run it on a compute
# node (srun / sbatch / Jupyterhub), never on stella*.
#
# One synthetic:
#   julia --project=. dev/moonRSGT/synth.jl --models=00817 \
#       --stla=-3.0 --stlo=-23.4 --evla=20 --evlo=15 --evdp=900 \
#       --mrr=1e15 --mtt=-5e14 --mpp=-5e14 --mrt=0 --mrp=0 --mtp=0 --out=synth.jld2
# Many (one per line of a CSV with header stla,stlo,evla,evlo,evdp,mrr,mtt,mpp,mrt,mrp,mtp):
#   julia --project=. dev/moonRSGT/synth.jl --models=all --requests=req.csv --out=synth.jld2
#
# --models   all | 00817 | 00001,00817 | 1:50 (positions in the sorted list of catalogue ids)
# --comps    ZRT (default) or any of Z N E R T
# --dt 1.0   --hd 10.0 (Gaussian half duration, s)   --quantity displacement|velocity|acceleration
# --dir      catalogue folder (default ENV MOON_RSGT_DIR, else /gpfs/scratch/nobuaki/moonRSGT/catalogue)
# Mij: GCMT order and convention (r = up, t = south, p = east), N·m.
#
# Output (JLD2 = HDF5, readable with h5py):
#   t        Vector (nt), s
#   data     Array (nt, ncomp, nmodel, nreq) in Julia = (nreq, nmodel, ncomp, nt) in h5py, m (m/s, m/s²)
#   models   "00001,00002,..."   comps "Z,R,T"   requests (nreq × 11) as given   delta_deg (nreq)

using JLD2, DelimitedFiles, Printf
include(joinpath(@__DIR__, "lib", "moonSynthetics.jl"))

function synthArgs(args)
    o = Dict{String,String}()
    for a in args
        startswith(a, "--") || error("unexpected argument $a")
        kv = split(a[3:end], '='; limit=2)
        o[kv[1]] = length(kv) == 2 ? kv[2] : "true"
    end
    return o
end

o = synthArgs(ARGS)
dir = get(o, "dir", moonRSGTDir())
allIds = moonModelIds(; dir=dir)
isempty(allIds) && error("no catalogue files in $dir")
msel = get(o, "models", "all")
ids = if msel == "all"
    allIds
elseif occursin(':', msel)
    a, b = parse.(Int, split(msel, ':')); allIds[a:min(b, end)]
else
    [length(s) < 5 && all(isdigit, s) ? lpad(s, 5, '0') : String(s) for s in split(msel, ',')]
end
const COLS = ["stla", "stlo", "evla", "evlo", "evdp", "mrr", "mtt", "mpp", "mrt", "mrp", "mtp"]
req = if haskey(o, "requests")
    A = readdlm(o["requests"], ',', Float64; skipstart=1)
    size(A, 2) == 11 || error("requests file: 11 columns expected ($(join(COLS, ',')))")
    A
else
    reshape([parse(Float64, o[c]) for c in COLS], 1, 11)
end
comps = Tuple(Symbol(string(c)) for c in get(o, "comps", "ZRT"))
dt = parse(Float64, get(o, "dt", "1.0")); hd = parse(Float64, get(o, "hd", "10.0"))
quantity = Symbol(get(o, "quantity", "displacement"))
out = get(o, "out", "synth.jld2")

nreq = size(req, 1)
data = nothing; t = nothing; Δ = zeros(nreq)
t0 = time()
for (im_, id) in enumerate(ids)
    cat = loadMoonRSGT(id; dir=dir)
    for k in 1:nreq
        r = req[k, :]
        tr = moonSynthetics(cat; station=(r[1], r[2]), source=(r[3], r[4], r[5]),
                            M=fromGCMT(r[6], r[7], r[8], r[9], r[10], r[11]),
                            comps=comps, dt=dt, halfDuration=hd, quantity=quantity)
        if isnothing(data)
            global t = tr.t
            global data = zeros(length(t), length(comps), length(ids), nreq)
        end
        for (ic, c) in enumerate(comps)
            data[:, ic, im_, k] .= tr.data[c]
        end
        Δ[k] = tr.Δdeg
    end
    (im_ % 50 == 0 || im_ == length(ids)) && @printf("%d / %d models (%.0f s)\n", im_, length(ids), time() - t0)
end
jldopen(out, "w") do f
    f["t"] = t
    f["data"] = data
    f["models"] = join(ids, ",")
    f["comps"] = join(string.(comps), ",")
    f["requests"] = req
    f["delta_deg"] = Δ
end
println("wrote ", out, ": ", length(ids), " models × ", nreq, " source-station pairs × ", length(comps), " components")
