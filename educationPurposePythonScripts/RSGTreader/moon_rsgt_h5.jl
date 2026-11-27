# moon_rsgt_h5.jl -- Moon synthetics from the HDF5 RSGT catalogue, standalone Julia (JLD2 + stdlib)
#
# Julia twin of moon_rsgt_py.py: reads rsgt_<id>.h5 (v1 or compressed v2, made by dev/moonRSGT/export_h5.jl) without
# flexOPT.  Same conventions as the Python reader:
#   station (stla, stlo), event (evla, evlo) in degrees, evdp in km;
#   moment tensor GCMT order/convention (r up, t south, p east) in N·m; force (fr, ft, fp) in N;
#   components Z N E R T; quantity "displacement" (m), "velocity" (m/s, default), "acceleration".
#
#   include("moon_rsgt_h5.jl")
#   MOON_H5_DIR[] = expanduser("~/Dropbox/stage_tutore_2026/catalogue_h5")
#   t, tr = moonH5Synthetics("00817"; stla=-3.0, stlo=-23.4, evla=20.0, evlo=15.0, evdp=900.0,
#                            Mrr=1e15, Mtt=-5e14, Mpp=-5e14)
#   tr["Z"], tr["R"], tr["T"]

using JLD2, LinearAlgebra, Printf

const MOON_H5_FORMATS = ("moonRSGT-h5 v1", "moonRSGT-h5 v2")
const MOON_H5_DIR = Ref(get(ENV, "MOON_RSGT_H5_DIR", expanduser("~/Dropbox/stage_tutore_2026/catalogue_h5")))

"Sorted model ids present in the folder (files still being copied are listed too: use finished ones)."
moonH5Ids(dir=MOON_H5_DIR[]) = sort([f[6:end-3] for f in readdir(dir) if startswith(f, "rsgt_") && endswith(f, ".h5")])

moonH5File(id; dir=MOON_H5_DIR[]) =
    (id isa AbstractString && isfile(id)) ? String(id) : joinpath(dir, "rsgt_" * lpad(string(id), 5, '0') * ".h5")

const _MOON_H5_LAST = Ref{Any}(nothing)

"""
Open one model: metadata only (v2: the compressed θ blocks are read on demand, a few MB per
synthetic; v1: the whole arrays, ≈ 0.6 GB).  The last model opened is kept, so loops over
events in one model are cheap.
"""
function moonH5Load(id; dir=MOON_H5_DIR[])
    f = moonH5File(id; dir=dir)
    c = _MOON_H5_LAST[]
    (c !== nothing && c.file == f) && return c
    isfile(f) || error("no file $f")
    c = jldopen(f, "r") do h
        fmt = haskey(h, "format") ? h["format"] : ""
        fmt in MOON_H5_FORMATS || error("$f: not a moonRSGT-h5 file")
        v2 = fmt == "moonRSGT-h5 v2"
        (file=f, v2=v2, id=h["model_id"],
         S=v2 ? nothing : complex.(h["strain_re"], h["strain_im"]),      # v1: (nθ, nr, 10, nω)
         D=v2 ? nothing : complex.(h["disp_re"], h["disp_im"]),
         tstart=v2 ? h["theta_block_start"] : Int[],
         blocks=Dict{Int,Any}(),
         ps=h["pairs_strain_src"], pc=h["pairs_strain_comp"], ds=h["pairs_disp_src"], dc=h["pairs_disp_comp"],
         ω=h["omega"], ωI=h["omega_i"], tlen=h["tlen"],
         r=h["r_km"], layer=h["layer"], θ=h["theta_deg"], R=h["R_km"], maxDepth=h["max_depth_km"])
    end
    _MOON_H5_LAST[] = c
    return c
end

"Strain and displacement samples at θ indices `tidx` (contiguous): arrays (length(tidx), nr, np, nω)."
function _h5samples(c, tidx::UnitRange{Int})
    c.v2 || return (c.S[tidx, :, :, :], c.D[tidx, :, :, :])
    ks = unique([searchsortedlast(c.tstart, i) for i in tidx])
    length(c.blocks) + length(ks) > 8 && empty!(c.blocks)        # keep a few decompressed blocks
    for k in ks
        haskey(c.blocks, k) && continue
        tag = lpad(k, 3, '0')
        c.blocks[k] = jldopen(c.file, "r") do h
            (complex.(h["strain_re_" * tag], h["strain_im_" * tag]), complex.(h["disp_re_" * tag], h["disp_im_" * tag]))
        end
    end
    pS = Any[]; pD = Any[]
    for i in tidx
        k = searchsortedlast(c.tstart, i); j = i - c.tstart[k] + 1
        push!(pS, c.blocks[k][1][j:j, :, :, :]); push!(pD, c.blocks[k][2][j:j, :, :, :])
    end
    S = Base.cat(pS...; dims=1); D = Base.cat(pD...; dims=1)
    return S, D
end

# ---------------------------------------------------------------- geometry (as sgtGeometry1D.jl)
_h5sph(θ, φ) = [sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)]
_h5basis(θ, φ) = hcat(_h5sph(θ, φ), [cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ)], [-sin(φ), cos(φ), 0.0])
_h5Rz(a) = [cos(a) -sin(a) 0.0; sin(a) cos(a) 0.0; 0.0 0.0 1.0]
_h5Ry(a) = [cos(a) 0.0 sin(a); 0.0 1.0 0.0; -sin(a) 0.0 cos(a)]
_h5pole(θ, φ) = _h5Ry(-θ) * _h5Rz(-φ)

function _h5dir(comp::Char, src, rcv)
    comp == 'Z' && return [1.0, 0.0, 0.0]
    comp == 'N' && return [0.0, -1.0, 0.0]
    comp == 'E' && return [0.0, 0.0, 1.0]
    E = _h5basis(rcv...); s = _h5sph(src...); rr = _h5sph(rcv...)
    t = cross(s, rr); norm(t) < 1e-12 && error("event below the station: use comps=\"ZNE\"")
    t ./= norm(t); L = cross(t, rr)
    comp == 'R' && return E' * L
    comp == 'T' && return E' * t
    error("component must be one of Z N E R T")
end

function _h5lagrange(x0, x)
    w = ones(length(x))
    for i in eachindex(x), j in eachindex(x)
        j != i && (w[i] *= (x0 - x[j]) / (x[i] - x[j]))
    end
    return w
end

"Epicentral distance (degrees) between (lat1, lon1) and (lat2, lon2)."
moonH5Distance(lat1, lon1, lat2, lon2) =
    rad2deg(acos(clamp(dot(_h5sph(deg2rad(90 - lat1), deg2rad(lon1)), _h5sph(deg2rad(90 - lat2), deg2rad(lon2))), -1.0, 1.0)))

# ---------------------------------------------------------------- spectra and time series

"Complex spectra (m·s, Dirac source) at the catalogue frequencies: (ω, ωI, tlen, Dict comp => spectrum, Δdeg)."
function moonH5Spectra(id; stla, stlo, evla, evlo, evdp, Mrr=0.0, Mtt=0.0, Mpp=0.0, Mrt=0.0, Mrp=0.0, Mtp=0.0,
                       force=nothing, comps="ZRT", dir=MOON_H5_DIR[])
    c = moonH5Load(id; dir=dir)
    0.0 <= evdp <= c.maxDepth + 1e-9 || error(@sprintf("evdp %.2f km outside 0 - %.2f km", evdp, c.maxDepth))
    rcv = (deg2rad(90 - stla), deg2rad(stlo)); src = (deg2rad(90 - evla), deg2rad(evlo)); rq = c.R - evdp
    Rs = _h5pole(rcv...)
    Bm = _h5basis(0.0, 0.0)' * Rs * _h5basis(rcv...)
    q = Rs * _h5sph(src...)
    θp = acos(clamp(q[3] / norm(q), -1.0, 1.0)); φp = atan(q[2], q[1])
    T = _h5basis(src...)' * Rs' * _h5basis(θp, φp)
    Qφ = [1.0 0.0 0.0; 0.0 cos(φp) sin(φp); 0.0 -sin(φp) cos(φp)]
    # depth: block (physical layer) holding rq, the upper one on a discontinuity
    L = c.layer
    starts = vcat(1, [i for i in 2:length(L) if L[i] != L[i-1]])
    stops = vcat(starts[2:end] .- 1, length(L))
    a, b = last([(a, b) for (a, b) in zip(starts, stops) if c.r[a] - 1e-9 <= rq <= c.r[b] + 1e-9])
    rv = c.r[a:b]; n = length(rv); qo = min(3, n - 1)
    j = clamp(searchsortedlast(rv, rq), 1, n); lo = clamp(j - (qo - 1) ÷ 2, 1, n - qo)
    rows = (a + lo - 1):(a + lo - 1 + qo)
    wr = _h5lagrange(rq, c.r[rows])
    # distance
    θd = rad2deg(θp); g = c.θ; nt = length(g); j = searchsortedlast(g, θd)
    tidx, wt = if 1 <= j <= nt && abs(g[j] - θd) <= 1e-12 * max(1.0, abs(θd))
        (j:j, [1.0])
    else
        lo = clamp(j - 2, 1, nt - 5); (lo:lo+5, _h5lagrange(θd, g[lo:lo+5]))
    end
    nω = length(c.ω)
    Ssel, Dsel = _h5samples(c, first(tidx):last(tidx))            # (nθsel, nr, np, nω)
    Sw = zeros(ComplexF64, nω, size(Ssel, 3)); Dw = zeros(ComplexF64, nω, size(Dsel, 3))
    for (ia, ir) in enumerate(rows), ib in eachindex(tidx)
        w = wr[ia] * wt[ib]
        Sw .+= w .* transpose(view(Ssel, ib, ir, :, :))
        Dw .+= w .* transpose(view(Dsel, ib, ir, :, :))
    end
    Mm = [Mrr Mrt Mrp; Mrt Mtt Mtp; Mrp Mtp Mpp] ./ 1e18           # code units 1e18 N·m
    F = isnothing(force) ? zeros(3) : collect(Float64, force) ./ 1e15
    out = Dict{String,Vector{ComplexF64}}()
    for comp in comps
        bv = Qφ * (Bm * _h5dir(comp, src, rcv))
        e6 = zeros(ComplexF64, nω, 6); u3 = zeros(ComplexF64, nω, 3)
        for ip in eachindex(c.ps); e6[:, c.pc[ip]] .+= bv[c.ps[ip]] .* Sw[:, ip]; end
        for ip in eachindex(c.ds); u3[:, c.dc[ip]] .+= bv[c.ds[ip]] .* Dw[:, ip]; end
        sp = zeros(ComplexF64, nω)
        for k in 1:nω
            e = view(e6, k, :)
            Ep = [e[1] e[4] e[5]; e[4] e[2] e[6]; e[5] e[6] e[3]]
            sp[k] = 1e3 * (sum((T * Ep * T') .* Mm) + sum((T * u3[k, :]) .* F))
        end
        out[string(comp)] = sp
    end
    return (ω=c.ω, ωI=c.ωI, tlen=c.tlen, spec=out, Δdeg=moonH5Distance(evla, evlo, stla, stlo))
end

"Real time series of one spectrum (same transform as moon_rsgt_py.time_series)."
function moonH5TimeSeries(ω, ωI, tlen, spec; dt=1.0, stf="step", half_duration=10.0, quantity="velocity", taper=0.2)
    p = Dict("displacement" => 0, "velocity" => 1, "acceleration" => 2)[quantity]
    ωc = ω .- im * ωI
    A = ComplexF64.(spec)
    stf == "step" ? (A ./= (im .* ωc)) : stf == "impulse" || error("stf must be \"step\" or \"impulse\"")
    p > 0 && (A .*= (im .* ωc) .^ p)
    half_duration > 0 && (A .*= exp.(-(ω .* (half_duration / 1.628)) .^ 2 ./ 2))
    nω = length(ω); w = ones(nω)
    if taper > 0
        k0 = floor(Int, (1 - taper) * (nω - 1))
        for k in k0:nω-1
            w[k+1] = 0.5 * (1 + cos(π * (k - k0) / max(1, nω - 1 - k0)))
        end
    end
    t = collect(0.0:dt:(tlen - dt / 2))
    u = [0.5 * w[1] * real(A[1]) + sum(w[k] * real(A[k] * cis(ω[k] * ti)) for k in 2:nω) for ti in t]
    return t, u .* (2 / tlen) .* exp.(ωI .* t)
end

"Time-domain synthetics of one model: (t, Dict comp => trace)."
function moonH5Synthetics(id; comps="ZRT", quantity="velocity", dt=1.0, half_duration=10.0, stf="step",
                          dir=MOON_H5_DIR[], kwargs...)
    s = moonH5Spectra(id; comps=comps, dir=dir, kwargs...)
    tr = Dict{String,Vector{Float64}}(); t = Float64[]
    for c in comps
        t, tr[string(c)] = moonH5TimeSeries(s.ω, s.ωI, s.tlen, s.spec[string(c)]; dt=dt, stf=stf,
                                            half_duration=half_duration, quantity=quantity)
    end
    return t, tr
end

"GCMT moment tensor (NamedTuple Mrr … Mtp, N·m if M0 is) of a double couple; angles in degrees."
function moonH5DoubleCouple(strike, dip, rake; M0=1.0)
    f, d, l = deg2rad(strike), deg2rad(dip), deg2rad(rake)
    return (Mrr=M0 * sin(2d) * sin(l),
            Mtt=-M0 * (sin(d) * cos(l) * sin(2f) + sin(2d) * sin(l) * sin(f)^2),
            Mpp=M0 * (sin(d) * cos(l) * sin(2f) - sin(2d) * sin(l) * cos(f)^2),
            Mrt=-M0 * (cos(d) * cos(l) * cos(f) + cos(2d) * sin(l) * sin(f)),
            Mrp=M0 * (cos(d) * cos(l) * sin(f) - cos(2d) * sin(l) * cos(f)),
            Mtp=-M0 * (sin(d) * cos(l) * cos(2f) + 0.5 * sin(2d) * sin(l) * sin(2f)))
end

"(lat, lon) at distance Δ (deg) and azimuth az (deg, clockwise from north) from (lat, lon)."
function moonH5PointAt(lat, lon, Δ, az)
    la, lo, d, a = deg2rad.((lat, lon, Δ, az))
    la2 = asin(sin(la) * cos(d) + cos(la) * sin(d) * cos(a))
    lo2 = lo + atan(sin(a) * sin(d) * cos(la), cos(d) - sin(la) * sin(la2))
    return rad2deg(la2), mod(rad2deg(lo2) + 180, 360) - 180
end
