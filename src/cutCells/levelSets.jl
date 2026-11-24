# src/cutCells/levelSets.jl
#
# Signed description of a free surface or an internal discontinuity on a
# Cartesian grid. Convention used everywhere in src/cutCells/:
#
#   φ <= 0   side "−"  (the solid for a free surface; material 0 for a discontinuity)
#   φ >  0   side "+"  (void above a free surface; material 1 for a discontinuity)
#   n = ∇φ/|∇φ|        unit normal pointing from "−" to "+" (outward for a free surface)
#
# φ is sampled at the grid nodes. Inside every simplex of the Kuhn
# decomposition of a cell (cellMoments.jl) φ is linear, so the interface is
# replaced by a plane facet per simplex: the curvature inside one cell is a
# higher-order effect (Mizutani 2002; Ovcharenko 2015, §3.3).

"""
    CutLevelSet(φ, axes)

Node values `φ` on the Cartesian grid `axes = (x, z)` / `(x, y, z)` / `(x,)`,
each axis uniformly spaced.
"""
struct CutLevelSet{N}
    φ::Array{Float64,N}
    axes::NTuple{N,Vector{Float64}}
    spacing::NTuple{N,Float64}
    function CutLevelSet(φ::AbstractArray{<:Real,N}, axes::NTuple{N,AbstractVector}) where {N}
        size(φ) == map(length, axes) ||
            throw(DimensionMismatch("φ has size $(size(φ)), axes give $(map(length, axes))"))
        spacing = map(axes) do a
            length(a) >= 2 || throw(ArgumentError("each axis needs at least two nodes"))
            h = Float64(a[2] - a[1])
            all(isapprox.(diff(a), h; rtol=1e-9)) ||
                throw(ArgumentError("axes must be uniformly spaced"))
            h
        end
        new{N}(Float64.(φ), map(a -> Float64.(collect(a)), axes), spacing)
    end
end

nodeCoordinates(ls::CutLevelSet{N}, I::CartesianIndex{N}) where {N} =
    SVector{N,Float64}(ntuple(d -> ls.axes[d][I[d]], N))

"φ(node) = f(SVector of node coordinates)."
function levelSetFromFunction(f, axes::NTuple{N,AbstractVector}) where {N}
    φ = Array{Float64,N}(undef, map(length, axes))
    for I in CartesianIndices(φ)
        φ[I] = f(SVector{N,Float64}(ntuple(d -> axes[d][I[d]], N)))
    end
    CutLevelSet(φ, axes)
end

"""
    levelSetFromHeight(h, axes)

Free surface / interface given as a height: 2-D `z = h(x)`, 3-D `z = h(x, y)`
(last axis vertical, up). φ = z - h, i.e. solid below. Exact zero crossing on
vertical grid lines; along other grid lines the crossing is the linear
interpolation of φ.
"""
levelSetFromHeight(h, axes::NTuple{2,AbstractVector}) =
    levelSetFromFunction(r -> r[2] - h(r[1]), axes)
levelSetFromHeight(h, axes::NTuple{3,AbstractVector}) =
    levelSetFromFunction(r -> r[3] - h(r[1], r[2]), axes)

"""
    levelSetFromPointsNormals(points, normals, axes)

Interface sampled by points with unit normals (pointing to the "+" side), e.g.
surface points from GeoPoints with their normal vectors. Each node takes
φ = n_p · (x - p) of its nearest interface point p: exact for a plane, and the
signed distance to the tangent plane of the nearest sample otherwise. Brute
force (Nnodes × Npoints); pass only the points near the grid.
"""
function levelSetFromPointsNormals(points::AbstractVector{SVector{N,Float64}},
                                   normals::AbstractVector{SVector{N,Float64}},
                                   axes::NTuple{N,AbstractVector}) where {N}
    length(points) == length(normals) ||
        throw(DimensionMismatch("points and normals differ in length"))
    isempty(points) && throw(ArgumentError("no interface points"))
    unit = [normalize(n) for n in normals]
    levelSetFromFunction(axes) do x
        best = 0
        bestDistance = Inf
        for k in eachindex(points)
            d = sum(abs2, x - points[k])
            d < bestDistance && ((best, bestDistance) = (k, d))
        end
        dot(unit[best], x - points[best])
    end
end

"""
    interfaceFromLocalCoords(coords)

`(points, normals)` from GeoPoints `localCoord2D`/`localCoord3D` values (or any
objects with fields `xz`/`xyz` and `normalVector`), for
`levelSetFromPointsNormals`. The normal must be the interface normal at that
point, oriented towards the "+" side.
"""
function interfaceFromLocalCoords(coords)
    first_ = first(coords)
    field = hasproperty(first_, :xyz) ? :xyz : hasproperty(first_, :xz) ? :xz :
        throw(ArgumentError("coordinates need a field xz or xyz"))
    points = [SVector{length(getproperty(c, field)),Float64}(getproperty(c, field)) for c in coords]
    normals = [SVector{length(c.normalVector),Float64}(c.normalVector) for c in coords]
    return points, normals
end

"Multilinear interpolation of φ at a point `x` inside the grid."
function levelSetValue(ls::CutLevelSet{N}, x::SVector{N,Float64}) where {N}
    base, frac = _cellLocate(ls, x)
    v = 0.0
    for c in CartesianIndices(ntuple(_ -> 0:1, N))
        w = 1.0
        for d in 1:N
            w *= c[d] == 1 ? frac[d] : 1 - frac[d]
        end
        w == 0.0 && continue
        v += w * ls.φ[base + c]
    end
    v
end

"""
    levelSetNormal(ls, x)

Unit normal ∇φ/|∇φ| at `x`: gradient of the multilinear interpolant of φ.
"""
function levelSetNormal(ls::CutLevelSet{N}, x::SVector{N,Float64}) where {N}
    base, frac = _cellLocate(ls, x)
    g = zeros(Float64, N)
    for c in CartesianIndices(ntuple(_ -> 0:1, N))
        value = ls.φ[base + c]
        for d in 1:N
            w = (c[d] == 1 ? 1.0 : -1.0) / ls.spacing[d]
            for e in 1:N
                e == d && continue
                w *= c[e] == 1 ? frac[e] : 1 - frac[e]
            end
            g[d] += w * value
        end
    end
    n = norm(g)
    n > 0 || throw(ArgumentError("∇φ vanishes at $x"))
    SVector{N,Float64}(g ./ n)
end

function _cellLocate(ls::CutLevelSet{N}, x::SVector{N,Float64}) where {N}
    idx = ntuple(N) do d
        t = (x[d] - ls.axes[d][1]) / ls.spacing[d]
        n = length(ls.axes[d])
        (-1e-9 <= t <= n - 1 + 1e-9) || throw(DomainError(x, "outside the grid"))
        clamp(floor(Int, t), 0, n - 2)
    end
    frac = SVector{N,Float64}(ntuple(d -> clamp(
        (x[d] - ls.axes[d][1]) / ls.spacing[d] - idx[d], 0.0, 1.0), N))
    return CartesianIndex(idx .+ 1), frac   # lower corner node of the containing cell
end
