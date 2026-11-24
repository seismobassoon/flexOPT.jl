# dev/moonRSGT/lib/moonCatalogue.jl
#
# Receiver-side strain Green's tensor (RSGT) catalogue of one Moon model:
# by reciprocity, the strain at any source point Q due to a unit force at a
# surface station R gives the station displacement for any moment tensor at
# Q (u_n(R) = ε^{R,n}(Q) : M), and the displacement at Q gives it for any
# force at Q (impacts, u_n(R) = u^{R,n}(Q) · F).  In a 1-D model every
# station is equivalent, so ONE pole disk per model is enough:
#
#   computeSGTDisk(...; family = :force, rs_km = R)   (dev/spherical1DelasticKernel/lib/sgtDisk1D.jl)
#   stored on   depths 0, dz, 2dz, ... down to the top of the fluid core,
#               + both sides of every discontinuity above it,
#           ×   epicentral distances 0:dθ:180°,
#           ×   frequencies 0 ... fmax (Δf = 1/tlen),
#   10 strain pairs + 5 displacement pairs, ComplexF32.
#
# Reconstruction for any station / source / mechanism: lib/moonSynthetics.jl
# (sgtAtAllω: local Lagrange, order 3 in depth inside each physical layer,
# order 5 in distance, explicit 3×3 rotations).

using JLD2, Printf

const MOON_RSGT_VERSION = "v1-2026-11-24"

"Default catalogue parameters (fmax 0.05 Hz, tlen 3276.8 s, elastic, chosen 2026-11-24)."
moonDefaultParams() = Dict{String,Any}(
    "version" => MOON_RSGT_VERSION,
    "scheme" => "OPT3",          # radial operator (DSMconv also possible)
    "h" => 5.0,                  # km, radial grid step of the solver (≥ 12 intervals per zone)
    "dz" => 20.0,                # km, catalogue depth step
    "dθ" => 0.5,                 # degrees, catalogue distance step
    "tlen" => 3276.8,            # s
    "fmax" => 0.05,              # Hz
    "omegaiTlen" => 1e-2,        # wrap-around damping, ω_I = -ln(omegaiTlen)/tlen (classic DSM value)
    "fitTol" => 0.005,           # km/s, g/cm³: cubic-zone fit of the .nd lines (lib/moonModel.jl)
    "lSlowness" => 2.0,          # lmax(ω) = ω R / lSlowness + lExtra  (min phase velocity, km/s)
    "lExtra" => 100,
)

"""
    catalogueRadii(info; dz) -> (r_km, side, layer)

Disk radii (increasing r): depths 0:dz:… inside every solid layer above the
fluid core, plus the top and bottom of each layer; a discontinuity appears
twice (side :lower = layer below, :upper = layer above).  `layer` is the
physical-layer index (1 = top) used as the interpolation block, so the
depth interpolation never crosses a discontinuity but does cross the
artificial cubic-fit zone boundaries inside a layer.
"""
function catalogueRadii(info; dz::Real)
    R = info.R
    rows = NamedTuple[]
    for (k, L) in enumerate(info.layers)
        L.fluid && break                         # stop at the top of the fluid core
        d = [L.dtop]
        for x in dz*ceil(L.dtop / dz):dz:L.dbot
            (x - L.dtop > 1e-6 && L.dbot - x > 1e-6) && push!(d, x)
        end
        push!(d, L.dbot)
        for x in d
            # top of a layer: take the zone below the boundary (:lower);
            # elsewhere (inside, or its bottom = discontinuity seen from above): :upper
            push!(rows, (r=R - x, side=(x == L.dtop ? :lower : :upper), layer=k))
        end
    end
    sort!(rows; by=q -> (q.r, -q.layer))         # increasing r; at a discontinuity the deeper layer first
    return (r_km=[q.r for q in rows], side=[q.side for q in rows], layer=[q.layer for q in rows])
end

"Deepest catalogue depth (km) of a model: top of its fluid core (or of the deepest solid layer stored)."
catalogueMaxDepth(info) = isnan(info.cmbDepth) ? info.layers[end].dbot : info.cmbDepth

"""
    moonSystems(model, params; recipes, configs) -> (psv, sh)

P-SV (+ fluid core in q = P/ω) and SH global systems of `params["scheme"]`
on a uniform-per-zone grid of step `params["h"]` (dev/spherical1Delastic/lib).
"""
function moonSystems(model, params; recipes, configs)
    gP = elasticRadialGrid(model; h=params["h"], region=:all)
    gS = elasticRadialGrid(model; h=params["h"], region=:outerSolid)
    kw = (scheme=params["scheme"], configs=configs, recipes=recipes, label=params["scheme"])
    return elasticGlobalSystem(:PSV, model, gP; kw...), elasticGlobalSystem(:SH, model, gS; kw...)
end

"lmax(ω) used for every disk of the catalogue."
moonLmaxOf(params, R) = ω -> ceil(Int, ω * R / params["lSlowness"]) + params["lExtra"]

"Same parameters (ignoring bookkeeping keys)?"
_sameParams(a, b) = all(k -> get(a, k, nothing) == get(b, k, nothing), keys(moonDefaultParams()))

"""
    buildMoonRSGT(ndFile, outFile; params=moonDefaultParams(), recipes, configs, force=false)

Compute and save the RSGT catalogue of one model (JLD2: `catalogue`, `info`,
`model` (= the .nd table), `params`, `timing`).  Skips the model when
`outFile` already exists with the same parameters (resume after a crash or
time limit), unless `force`.  The file is written to `outFile * ".part"`
and renamed at the end, so a killed job never leaves a truncated catalogue.
"""
function buildMoonRSGT(ndFile::AbstractString, outFile::AbstractString;
                       params=moonDefaultParams(), recipes, configs, force::Bool=false, verbose::Bool=true)
    if !force && isfile(outFile)
        old = try
            JLD2.load(outFile, "params")
        catch
            nothing
        end
        if !isnothing(old) && _sameParams(old, params)
            verbose && println("skip (already done): ", outFile)
            return :skipped
        end
    end
    t0 = time()
    mm = moonModelFromND(ndFile; fitTol=params["fitTol"])
    info = mm.info
    verbose && println(moonModelSummary(info))
    psv, sh = moonSystems(mm.model, params; recipes=recipes, configs=configs)
    t1 = time()
    rad = catalogueRadii(info; dz=params["dz"])
    disk = computeSGTDisk(psv, sh, mm.model; family=:force, rs_km=info.R,
                          radii=(r_km=rad.r_km, side=rad.side),
                          θdeg=0:params["dθ"]:180, tlen=params["tlen"], fmax=params["fmax"],
                          omegaiTlen=params["omegaiTlen"], p=sourceOrderOf(params["scheme"]),
                          lmaxOf=moonLmaxOf(params, info.R), lcap=5000, verbose=verbose)
    t2 = time()
    # store: Float32, interpolation blocks = physical layers
    cat = merge(disk, (seg=rad.layer, strain=ComplexF32.(disk.strain), disp=ComplexF32.(disk.disp)))
    tl = disk.tail[end, :]
    timing = Dict("systems_s" => t1 - t0, "disk_s" => t2 - t1, "tailMedianFmax" => sort(tl)[(length(tl) + 1) ÷ 2],
                  "tailMaxFmax" => maximum(tl))
    mkpath(dirname(abspath(outFile)))
    tmp = outFile * ".part"
    jldsave(tmp; catalogue=cat, info=info, model=mm.table, params=params, timing=timing,
            maxDepth=catalogueMaxDepth(info))
    mv(tmp, outFile; force=true)
    nb = 8 * (length(cat.strain) + length(cat.disp))
    verbose && @printf("saved %s: %d ω × %d r × %d θ, %.0f MB; systems %.0f s, disk %.0f s; tail at fmax median %.1e max %.1e\n",
                       outFile, length(cat.ω), length(cat.r_km), length(cat.θdeg), nb / 1e6,
                       t1 - t0, t2 - t1, timing["tailMedianFmax"], timing["tailMaxFmax"])
    return :done
end
