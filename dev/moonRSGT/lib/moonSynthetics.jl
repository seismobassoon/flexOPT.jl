# dev/moonRSGT/lib/moonSynthetics.jl -- synthetics from the Moon RSGT catalogue
#
# Light on purpose (for students): needs only JLD2 + the stdlib.  It reuses
# the geometry / interpolation code of dev/spherical1DelasticKernel/lib
# (single source of truth); none of the solver code is loaded.
#
#   include("dev/moonRSGT/lib/moonSynthetics.jl")
#   cat = loadMoonRSGT(817)                                  # rsgt_00817.jld2 in moonRSGTDir()
#   M   = doubleCouple(strike=30, dip=60, rake=90, M0=1e15)  # N·m, (Mrr, Mrθ, Mrφ, Mθθ, Mθφ, Mφφ)
#   tr  = moonSynthetics(cat; station=(lat, lon), source=(lat, lon, depth_km), M=M)
#   tr.t, tr.data[:Z], tr.data[:R], tr.data[:T]               # s, m
#
# Conventions: local basis (r, θ, φ) = (up, south, east), i.e. Harvard/GCMT
# (r, t, p).  Moment tensor in N·m, force in N, output in m (or m/s, m/s²).
# Components :Z (up), :N, :E, :R (horizontal, away from the source), :T
# (normal to the source-station great circle, see receiverDirection).

using JLD2, LinearAlgebra, Printf

const _MOON_KER_LIB = normpath(joinpath(@__DIR__, "..", "..", "spherical1DelasticKernel", "lib"))

if !isdefined(@__MODULE__, :fornbergWeights)
    # same algorithm as dev/spherical1Delastic/lib/sphericalTI_SH_step3.jl (not loaded here: it needs the solver)
    "Fornberg (1988) weights for the m-th derivative at x0 on nodes x."
    function fornbergWeights(x0::Real, x::AbstractVector, m::Integer=1)
        n = length(x)
        c = zeros(n, m + 1)
        c1 = 1.0; c4 = x[1] - x0; c[1, 1] = 1.0
        for i in 2:n
            mn = min(i - 1, m); c2 = 1.0; c5 = c4; c4 = x[i] - x0
            for j in 1:i-1
                c3 = x[i] - x[j]; c2 *= c3
                if j == i - 1
                    for k in mn:-1:1
                        c[i, k+1] = c1 * (k * c[i-1, k] - c5 * c[i-1, k+1]) / c2
                    end
                    c[i, 1] = -c1 * c5 * c[i-1, 1] / c2
                end
                for k in mn:-1:1
                    c[j, k+1] = (c4 * c[j, k+1] - k * c[j, k]) / c3
                end
                c[j, 1] = c4 * c[j, 1] / c3
            end
            c1 = c2
        end
        return c[:, m+1]
    end
end
isdefined(@__MODULE__, :sgtPairs) || include(joinpath(_MOON_KER_LIB, "sgtHarmonics1D.jl"))
isdefined(@__MODULE__, :diskθWeights) || include(joinpath(_MOON_KER_LIB, "sgtDisk1D.jl"))
isdefined(@__MODULE__, :PointGeo) || include(joinpath(_MOON_KER_LIB, "sgtGeometry1D.jl"))
isdefined(@__MODULE__, :sgtAtAllω) || include(joinpath(_MOON_KER_LIB, "sgtInterp1D.jl"))

"Catalogue folder: ENV[\"MOON_RSGT_DIR\"], default /gpfs/scratch/nobuaki/moonRSGT/catalogue."
moonRSGTDir() = get(ENV, "MOON_RSGT_DIR", "/gpfs/scratch/nobuaki/moonRSGT/catalogue")

"Catalogue file name of a model id (the .nd file stem, e.g. 817 or \"00817\")."
moonRSGTFile(id::Integer; dir=moonRSGTDir()) = joinpath(dir, @sprintf("rsgt_%05d.jld2", id))
moonRSGTFile(id::AbstractString; dir=moonRSGTDir()) = joinpath(dir, "rsgt_$(id).jld2")

"Sorted list of the model ids (strings) present in the catalogue folder."
moonModelIds(; dir=moonRSGTDir()) =
    sort([m.captures[1] for m in (match(r"^rsgt_(.+)\.jld2$", f) for f in readdir(dir)) if !isnothing(m)])

"""
    loadMoonRSGT(id_or_file; dir=moonRSGTDir())

NamedTuple (disk, info, model, params, maxDepth, file).  `model` is the .nd
table (depth, vp, vs, ρ); `maxDepth` the deepest source depth available
(top of the fluid core).  About 0.5 GB in memory per model.
"""
function loadMoonRSGT(id::Union{Integer,AbstractString}; dir=moonRSGTDir())
    file = (id isa AbstractString && isfile(id)) ? String(id) : moonRSGTFile(id; dir=dir)
    isfile(file) || error("no catalogue file $file")
    d = JLD2.load(file)
    return (disk=d["catalogue"], info=d["info"], model=d["model"], params=d["params"],
            maxDepth=d["maxDepth"], file=file)
end

# ---------------------------------------------------------------------------
# sources
# ---------------------------------------------------------------------------

"Moment tensor in this code's order (Mrr, Mrθ, Mrφ, Mθθ, Mθφ, Mφφ) from GCMT order (Mrr, Mtt, Mpp, Mrt, Mrp, Mtp)."
fromGCMT(Mrr, Mtt, Mpp, Mrt, Mrp, Mtp) = (Mrr, Mrt, Mrp, Mtt, Mtp, Mpp)

"""
    doubleCouple(; strike, dip, rake, M0=1.0)

Moment tensor (N·m if M0 is) of a double couple, in this code's order
(Mrr, Mrθ, Mrφ, Mθθ, Mθφ, Mφφ) (Aki & Richards, in GCMT (r, t, p) axes).
Angles in degrees.
"""
function doubleCouple(; strike, dip, rake, M0=1.0)
    φ, δ, λ = deg2rad(strike), deg2rad(dip), deg2rad(rake)
    Mrr = M0 * sin(2δ) * sin(λ)
    Mtt = -M0 * (sin(δ) * cos(λ) * sin(2φ) + sin(2δ) * sin(λ) * sin(φ)^2)
    Mpp = M0 * (sin(δ) * cos(λ) * sin(2φ) - sin(2δ) * sin(λ) * cos(φ)^2)
    Mrt = -M0 * (cos(δ) * cos(λ) * cos(φ) + cos(2δ) * sin(λ) * sin(φ))
    Mrp = M0 * (cos(δ) * cos(λ) * sin(φ) - cos(2δ) * sin(λ) * cos(φ))
    Mtp = -M0 * (sin(δ) * cos(λ) * cos(2φ) + 0.5 * sin(2δ) * sin(λ) * sin(2φ))
    return fromGCMT(Mrr, Mtt, Mpp, Mrt, Mrp, Mtp)
end

"Epicentral distance (degrees) between two (lat, lon) points (degrees)."
function epicentralDistance(lat1, lon1, lat2, lon2)
    a = sphUnit(deg2rad(90 - lat1), deg2rad(lon1)); b = sphUnit(deg2rad(90 - lat2), deg2rad(lon2))
    return rad2deg(acos(clamp(dot(a, b), -1.0, 1.0)))
end

# ---------------------------------------------------------------------------
# spectra and time series
# ---------------------------------------------------------------------------

"""
    moonSpectra(cat; station, source, M=nothing, force=nothing, comps=(:Z, :R, :T), side=:upper)

Complex spectra (m · s, Dirac source time function) at the catalogue
frequencies `cat.disk.ω` (rad/s, evaluated at ω - iω_I).  `station = (lat,
lon)`, `source = (lat, lon, depth_km)` (degrees, km), `M` (N·m, this
code's order, see `doubleCouple` / `fromGCMT`) and/or `force` (N, (f_r, f_θ,
f_φ) = (up, south, east)).  `side` only matters for a source exactly on a
discontinuity (:upper = the layer above).
"""
function moonSpectra(cat; station, source, M=nothing, force=nothing, comps=(:Z, :R, :T), side::Symbol=:upper)
    (isnothing(M) && isnothing(force)) && error("give a moment tensor M and/or a force")
    R = cat.info.R
    lat_r, lon_r = station
    lat_s, lon_s, dep = source
    (0.0 <= dep <= cat.maxDepth + 1e-9) ||
        error(@sprintf("source depth %.2f km outside the catalogue (0 - %.2f km, top of the fluid core)", dep, cat.maxDepth))
    Δ = epicentralDistance(lat_s, lon_s, lat_r, lon_r)
    if Δ < 1e-3 && any(c -> c in (:R, :T), comps)
        error("source below the station: R/T undefined, use comps = (:Z, :N, :E)")
    end
    Rcv = pointFromLatLon(R, 0.0, lat_r, lon_r)
    Q = pointFromLatLon(R, dep, lat_s, lon_s)
    # code units: M in (g/cm³)(km/s)²km³ = 1e18 N·m, f in (g/cm³)(km/s)²km² = 1e15 N, u in km
    Mm = isnothing(M) ? zeros(3, 3) : mtMatrix(collect(Float64, M) ./ 1e18)
    F = isnothing(force) ? zeros(3) : collect(Float64, force) ./ 1e15
    out = Dict{Symbol,Vector{ComplexF64}}()
    for comp in comps
        n = receiverDirection(comp, Q, Rcv)
        E, U = sgtAtAllω(cat.disk, Rcv, n, Q; side=side)
        out[comp] = [1e3 * (sum(view(E, k, :, :) .* Mm) + sum(view(U, k, :) .* F)) for k in axes(E, 1)]
    end
    return (ω=cat.disk.ω, ωI=cat.disk.ωI, tlen=cat.disk.tlen, spec=out, Δdeg=Δ)
end

"""
    moonTimeSeries(ω, ωI, tlen, spec; dt=1.0, stf=:step, halfDuration=0.0,
                   quantity=:displacement, taper=0.2) -> (t, u)

Real time series of one spectrum (one-sided, exp(+iωt), undamped by
exp(ω_I t), cosine taper on the last `taper` fraction of the band, as
`diskSnapshot`).  `stf = :step` (Heaviside moment / force; default) or
`:impulse` (Dirac, the raw catalogue); `halfDuration` > 0 adds a Gaussian
source time function (σ = halfDuration/1.628, SPECFEM convention).
`quantity = :displacement | :velocity | :acceleration`.
"""
function moonTimeSeries(ω, ωI, tlen, spec::AbstractVector; dt::Real=1.0, stf::Symbol=:step,
                        halfDuration::Real=0.0, quantity::Symbol=:displacement, taper::Real=0.2)
    abs(ω[1]) < 1e-12 || error("the catalogue must start at f = 0")
    stf in (:step, :impulse) || error("stf must be :step or :impulse")
    p = quantity === :displacement ? 0 : quantity === :velocity ? 1 : quantity === :acceleration ? 2 :
        error("quantity must be :displacement, :velocity or :acceleration")
    ωc = ω .- im * ωI
    A = ComplexF64.(spec)
    stf === :step && (A ./= (im .* ωc))
    p > 0 && (A .*= (im .* ωc) .^ p)
    if halfDuration > 0
        σ = halfDuration / 1.628
        A .*= exp.(-(ω .* σ) .^ 2 ./ 2)
    end
    w = _taperWeights(length(ω), taper)
    t = collect(0.0:dt:(tlen - dt / 2))
    u = zeros(length(t))
    Δω = 2π / tlen
    for (i, ti) in enumerate(t)
        s = 0.5 * w[1] * real(A[1])
        for k in 2:length(ω)
            s += w[k] * real(A[k] * cis(ω[k] * ti))
        end
        u[i] = s * (Δω / π) * exp(ωI * ti)
    end
    return t, u
end

"""
    moonSynthetics(cat; station, source, M=nothing, force=nothing, comps=(:Z, :R, :T),
                   dt=1.0, stf=:step, halfDuration=10.0, quantity=:displacement, taper=0.2)

Time-domain synthetics (s; m, m/s or m/s²) of one catalogue model, see
`moonSpectra` (geometry, units) and `moonTimeSeries` (source time
function).  The catalogue is band-limited to fmax (0.05 Hz): keep
halfDuration ≳ 1/(2 fmax) = 10 s, or the taper alone shapes the pulse.
Returns (t, data::Dict{Symbol,Vector}, Δdeg, model, source, station).
"""
function moonSynthetics(cat; station, source, M=nothing, force=nothing, comps=(:Z, :R, :T),
                        dt::Real=1.0, stf::Symbol=:step, halfDuration::Real=10.0,
                        quantity::Symbol=:displacement, taper::Real=0.2, side::Symbol=:upper)
    sp = moonSpectra(cat; station=station, source=source, M=M, force=force, comps=comps, side=side)
    data = Dict{Symbol,Vector{Float64}}()
    t = Float64[]
    for c in comps
        t, data[c] = moonTimeSeries(sp.ω, sp.ωI, sp.tlen, sp.spec[c]; dt=dt, stf=stf,
                                    halfDuration=halfDuration, quantity=quantity, taper=taper)
    end
    return (t=t, data=data, Δdeg=sp.Δdeg, model=cat.info.file, source=source, station=station)
end

"""
    moonSyntheticsAllModels(ids=moonModelIds(); dir=moonRSGTDir(), kwargs...)

`moonSynthetics` for every model id (loaded one at a time, ~0.5 GB each);
Dict id => result.  Same keywords as `moonSynthetics`.
"""
function moonSyntheticsAllModels(ids=moonModelIds(); dir=moonRSGTDir(), verbose::Bool=true, kwargs...)
    out = Dict{String,Any}()
    for (k, id) in enumerate(ids)
        cat = loadMoonRSGT(id; dir=dir)
        out[string(id)] = moonSynthetics(cat; kwargs...)
        verbose && (k % 50 == 0 || k == length(ids)) && println("  $k / $(length(ids)) models")
    end
    return out
end
