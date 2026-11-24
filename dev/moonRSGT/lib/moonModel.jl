# dev/moonRSGT/lib/moonModel.jl
#
# TauP-style .nd Moon models (depth  vp  vs  ρ per line, km, km/s, g/cm³;
# a repeated depth = a discontinuity; vs = 0 = fluid) -> planet1D model
# (`planet1D.planet1DPSVmodel`, cubic polynomials in ξ = r/R per zone, the
# format every dev/spherical1Delastic routine reads).
#
# Why not the TauP path of src/planet1D (read1DModel)?  It turns every
# 10-km line interval into its own zone (fusion only when the slope changes
# by < 1 %, and the 3-decimal rounding of the .nd files defeats it), which
# gives ~100 zones of ≥ 12 grid intervals each and puts a zone end (one-sided
# strain derivatives) at almost every catalogue depth.  Here each physical
# layer (between two discontinuities) is fitted by ONE cubic in ξ, split
# recursively in halves only where the fit misses a .nd value by more than
# `fitTol` (km/s, g/cm³).  Over the 1000 models: 3-13 zones for fitTol =
# 0.005, 3-10 for 0.01 (crust, mantle pieces, fluid core).
#
# Nothing in src/ is modified; only the public struct constructor is used.

using LinearAlgebra, Printf

"""
    readMoonND(file) -> (depth, vp, vs, ρ)

Numeric lines with at least 4 numbers; anything else (comments, TauP
keywords such as `mantle`, `outer-core`) is skipped.  Depth increasing.
"""
function readMoonND(file::AbstractString)
    d = Float64[]; vp = Float64[]; vs = Float64[]; ρ = Float64[]
    for line in eachline(file)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        w = split(s)
        length(w) < 4 && continue
        x = tryparse.(Float64, w[1:4])
        any(isnothing, x) && continue
        push!(d, x[1]); push!(vp, x[2]); push!(vs, x[3]); push!(ρ, x[4])
    end
    length(d) >= 2 || error("$file: no model lines found")
    issorted(d) || error("$file: depths must be increasing")
    return (depth=d, vp=vp, vs=vs, ρ=ρ)
end

"Row ranges of the physical layers (split at repeated depths), top first."
function ndLayers(depth)
    starts = [1]
    for i in 2:length(depth)
        depth[i] == depth[i-1] && push!(starts, i)
    end
    stops = vcat(starts[2:end] .- 1, length(depth))
    return [a:b for (a, b) in zip(starts, stops) if b > a]      # a 1-line "layer" has no thickness
end

"Least-squares polynomial (degree ≤ 3, padded to 4 coefficients) of y(x); returns (coefs, max |misfit|)."
function _cubicFit(x, y)
    deg = min(3, length(x) - 1)
    A = [xi^k for xi in x, k in 0:deg]
    c = A \ y
    return vcat(c, zeros(3 - deg)), maximum(abs.(A * c .- y))
end

"""
    fitLayer(ξ, Y; tol) -> [(rows, coefs::Matrix 4×3, misfit)]

ξ increasing (bottom → top), Y = [vp vs ρ].  One cubic per piece; a piece
whose misfit exceeds `tol` is split in two halves sharing their middle node
(so the pieces stay continuous up to `tol`).
"""
function fitLayer(ξ, Y; tol)
    n = length(ξ)
    C = zeros(4, 3); mis = 0.0
    for j in 1:3
        c, m = _cubicFit(ξ, Y[:, j]); C[:, j] = c; mis = max(mis, m)
    end
    (mis <= tol || n <= 5) && return [(rows=1:n, coefs=C, misfit=mis)]
    mid = (n + 1) ÷ 2
    lo = fitLayer(ξ[1:mid], Y[1:mid, :]; tol=tol)
    hi = fitLayer(ξ[mid:end], Y[mid:end, :]; tol=tol)
    hi = [(rows=p.rows .+ (mid - 1), coefs=p.coefs, misfit=p.misfit) for p in hi]
    return vcat(lo, hi)
end

"""
    moonModelFromND(file; fitTol=0.005) -> (model, info, table)

`model`: planet1D.planet1DPSVmodel (zones bottom → top, isotropic, no Q).
`info`:  R, the physical layers (top first: dtop, dbot, fluid, zones) and
         the depth of the top of the uppermost fluid layer (`cmbDepth`,
         NaN if none), the worst fit misfit, and the file name.
`table`: the .nd values (depth, vp, vs, ρ).
"""
function moonModelFromND(file::AbstractString; fitTol::Real=0.005)
    t = readMoonND(file)
    R = t.depth[end]                         # TauP convention: the last depth is the centre
    t.depth[1] == 0.0 || error("$file: the first depth must be 0 (surface)")
    layers = ndLayers(t.depth)
    zones = NamedTuple[]                     # collected top → bottom, reversed at the end
    linfo = NamedTuple[]
    worst = 0.0
    for (k, rows) in enumerate(layers)
        fluid = all(iszero, t.vs[rows])
        any(iszero, t.vs[rows]) && !fluid && error("$file: layer $k mixes vs = 0 and vs > 0")
        rr = reverse(collect(rows))          # bottom → top
        ξ = (R .- t.depth[rr]) ./ R
        Y = hcat(t.vp[rr], t.vs[rr], t.ρ[rr])
        pieces = fitLayer(ξ, Y; tol=fitTol)
        nz0 = length(zones)
        for p in reverse(pieces)             # top piece first
            ξb, ξt = ξ[first(p.rows)], ξ[last(p.rows)]
            C = copy(p.coefs)
            fluid && (C[:, 2] .= 0.0)
            push!(zones, (ξb=ξb, ξt=ξt, C=C, fluid=fluid, layer=k))
            worst = max(worst, p.misfit)
        end
        push!(linfo, (dtop=t.depth[first(rows)], dbot=t.depth[last(rows)], fluid=fluid,
                      nzone=length(zones) - nz0))
    end
    reverse!(zones)                          # bottom → top (planet1D order)
    nz = length(zones)
    m = planet1D.planet1DPSVmodel(nz, Float64(R))
    for (iz, z) in enumerate(zones)
        m.solid_or_fluid[iz] = z.fluid ? "F" : "S"
        m.normalisedBottomRadius[iz] = z.ξb; m.normalisedTopRadius[iz] = z.ξt
        m.bottomRadius[iz] = z.ξb * R;       m.topRadius[iz] = z.ξt * R
        m.C_Vpv[iz, :] = z.C[:, 1]; m.C_Vph[iz, :] = z.C[:, 1]
        m.C_Vsv[iz, :] = z.C[:, 2]; m.C_Vsh[iz, :] = z.C[:, 2]
        m.C_ρ[iz, :] = z.C[:, 3]
        m.C_η[iz, :] = [1.0, 0.0, 0.0, 0.0]
        for f in (:C_Qκ, :C_Qμ, :C_QκPower, :C_QμPower)
            getfield(m, f)[iz, :] .= 0.0
        end
    end
    kf = findfirst(l -> l.fluid, linfo)
    cmbDepth = isnothing(kf) ? NaN : linfo[kf].dtop
    m.averagedPlanetCMBInKilometer = isnothing(kf) ? 0.0 : R - cmbDepth
    m.normalisedCMB = m.averagedPlanetCMBInKilometer / R
    m.solid_or_fluid[end] == "S" || error("$file: fluid top layer (ocean) not supported")
    info = (file=basename(file), R=Float64(R), layers=identity.(linfo), cmbDepth=cmbDepth,
            nzone=nz, fitTol=Float64(fitTol), fitMisfit=worst)
    return (model=m, info=info, table=t)
end

"One-line summary of a model (for logs)."
moonModelSummary(info) = @sprintf("%s: R = %.1f km, %d layers (%s), %d zones, CMB %.1f km, fit misfit %.4f",
    info.file, info.R, length(info.layers),
    join([@sprintf("%s %.1f-%.1f", l.fluid ? "F" : "S", l.dtop, l.dbot) for l in info.layers], ", "),
    info.nzone, info.cmbDepth, info.fitMisfit)
