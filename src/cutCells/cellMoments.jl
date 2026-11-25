# src/cutCells/cellMoments.jl
#
# Moments of the part of one grid cell that lies on a given side of the
# interface:
#
#   m[a_1+1, ..., a_N+1] = ∫_{cell ∩ {±φ <= 0}} y_1^a_1 ... y_N^a_N dy,
#   y ∈ [0, 1]^N the cell's local coordinates (multiply by Π spacing and use
#   y = (x - x_cellcorner) / spacing to go back to physical coordinates).
#
# Every OPT integrand (test function × interpolants × Taylor kernels) is a
# polynomial inside one cell, so the cut-cell volume integral is exactly a
# contraction of its polynomial coefficients with these moments. That is the
# only geometric input the cut-cell operators need (README.md, stage 2).
#
# Method: the cell is split into the N! simplices of the Kuhn (Freudenthal)
# decomposition; φ is linear on each, so each simplex is clipped exactly by
# a plane and the clipped piece is again a union/difference of simplices,
# integrated exactly by Grundmann–Möller (simplexIntegration.jl).

function _permutations(v::Vector{Int})
    length(v) <= 1 && return [copy(v)]
    out = Vector{Vector{Int}}()
    for i in eachindex(v)
        rest = vcat(v[1:i-1], v[i+1:end])
        for p in _permutations(rest)
            push!(out, vcat(v[i], p))
        end
    end
    out
end

const _KUHN_CACHE = Dict{Int,Vector{Vector{NTuple}}}()

"Kuhn simplices of [0,1]^N as lists of N+1 corner offsets (0/1 tuples)."
function kuhnSimplices(N::Int)
    get!(_KUHN_CACHE, N) do
        out = Vector{Vector{NTuple}}()
        for perm in _permutations(collect(1:N))
            corner = zeros(Int, N)
            verts = NTuple[Tuple(corner)]
            for d in perm
                corner[d] += 1
                push!(verts, Tuple(corner))
            end
            push!(out, verts)
        end
        out
    end
end

"""
    clippedSimplexPieces(V, φ) -> Vector{(vertices, sign)}

Region {φ <= 0} of the simplex `V` (φ linear, values `φ` at the vertices) as
a signed sum of simplices.
"""
function clippedSimplexPieces(V::Vector{SVector{N,Float64}}, φ::AbstractVector{Float64}) where {N}
    inside = [i for i in 1:N+1 if φ[i] <= 0]
    outside = [i for i in 1:N+1 if φ[i] > 0]
    pieces = Tuple{Vector{SVector{N,Float64}},Float64}[]
    isempty(inside) && return pieces
    isempty(outside) && return push!(pieces, (V, 1.0))
    crossing(i, j) = V[i] + (φ[i] / (φ[i] - φ[j])) * (V[j] - V[i])
    if length(inside) == 1
        a = inside[1]
        push!(pieces, (vcat([V[a]], [crossing(a, o) for o in outside]), 1.0))
    elseif length(outside) == 1
        o = outside[1]
        push!(pieces, (V, 1.0))
        push!(pieces, (vcat([V[o]], [crossing(o, a) for a in inside]), -1.0))
    else
        # N == 3, two vertices on each side: a prism, three tetrahedra.
        N == 3 || error("unexpected clipping configuration")
        a, b = inside
        c, d = outside
        A1, A2, A3 = V[a], crossing(a, c), crossing(a, d)
        B1, B2, B3 = V[b], crossing(b, c), crossing(b, d)
        push!(pieces, ([A1, A2, A3, B3], 1.0))
        push!(pieces, ([A1, A2, B2, B3], 1.0))
        push!(pieces, ([A1, B1, B2, B3], 1.0))
    end
    return pieces
end

"""
    cellMoments(φcorners, p; side=:minus)

Moments (see the file header) for tensor degree `p` per axis of the part of
the unit cell with φ <= 0 (`side=:minus`) or φ > 0 (`side=:plus`).
`φcorners` is the 2×…×2 array of φ at the cell corners (index 1 = lower).
"""
function cellMoments(φcorners::AbstractArray{<:Real,N}, p::Int; side::Symbol=:minus) where {N}
    size(φcorners) == ntuple(_ -> 2, N) ||
        throw(DimensionMismatch("φcorners must be 2×…×2"))
    side in (:minus, :plus) || throw(ArgumentError("side must be :minus or :plus"))
    s = side === :minus ? 1.0 : -1.0       # φ > 0  ⇔  -φ < 0
    m = zeros(Float64, ntuple(_ -> p + 1, N))
    for verts in kuhnSimplices(N)
        V = [SVector{N,Float64}(Float64.(v)) for v in verts]
        # :plus uses -φ; the set {φ == 0} has zero measure either way.
        φV = [s * Float64(φcorners[CartesianIndex(v .+ 1)]) for v in verts]
        for (piece, sign) in clippedSimplexPieces(V, φV)
            addSimplexMonomialMoments!(m, piece, p; sign=sign)
        end
    end
    return m
end

"Moments of the whole unit cell: Π 1/(a_d+1)."
fullCellMoments(N::Int, p::Int) =
    [prod(1.0 / I[d] for d in 1:N) for I in CartesianIndices(ntuple(_ -> p + 1, N))]

"""
    CutCellGeometry(ls; maxDegree=4, tolerance=1e-12, zigzag=false)

Cell-by-cell description of the interface `ls` (a `CutLevelSet`):

* `fraction[c]`  volume fraction of cell `c` on the "−" side (cell `c` spans
  nodes `c` … `c + 1`)
* `state[c]`     0 = entirely "+", 1 = entirely "−", 2 = cut
* `moments[c]`   "−"-side moments (tensor degree `maxDegree`) of every cut
  cell; the "+" side is `fullCellMoments - moments[c]`.

With `zigzag=true` the interface is replaced by its stair-step version
(`zigzagLevelSet`): no cut cells, and `levelSet` is the stair-step level set,
so `sideMask(geom.levelSet)` gives the material nodes of the staircase. The
same operators then give the classical zigzag free surface, for comparison.
"""
struct CutCellGeometry{N}
    levelSet::CutLevelSet{N}
    maxDegree::Int
    fraction::Array{Float64,N}
    state::Array{Int8,N}
    moments::Dict{CartesianIndex{N},Array{Float64,N}}
end

"""
    zigzagLevelSet(ls)

Stair-step ("zigzag") version of the interface: a cell is material iff its
2^N corners all have φ <= 0, a node is material iff it is a corner of such a
cell. Returns a level set with the same |φ| (so normals stay close to the
true ones) and the sign of that node classification.
"""
function zigzagLevelSet(ls::CutLevelSet{N}) where {N}
    raw = ls.φ .<= 0.0
    cellsShape = size(ls.φ) .- 1
    corners = CartesianIndices(ntuple(_ -> 0:1, N))
    solid = falses(size(ls.φ))
    for c in CartesianIndices(cellsShape)
        all(raw[c + k] for k in corners) || continue
        for k in corners
            solid[c + k] = true
        end
    end
    tiny = 1e-12 * maximum(ls.spacing)
    φ = [solid[I] ? -max(abs(ls.φ[I]), tiny) : max(abs(ls.φ[I]), tiny) for I in CartesianIndices(ls.φ)]
    CutLevelSet(φ, ls.axes)
end

function CutCellGeometry(ls0::CutLevelSet{N}; maxDegree::Int=4, tolerance::Float64=1e-12,
                         zigzag::Bool=false) where {N}
    # zigzag = true: stair-step surface, every cell entirely material or entirely void
    ls = zigzag ? zigzagLevelSet(ls0) : ls0
    cellsShape = size(ls.φ) .- 1
    fraction = zeros(Float64, cellsShape)
    state = zeros(Int8, cellsShape)
    moments = Dict{CartesianIndex{N},Array{Float64,N}}()
    corners = CartesianIndices(ntuple(_ -> 0:1, N))
    for c in CartesianIndices(cellsShape)
        φc = [ls.φ[c + k] for k in corners]
        if all(<=(0.0), φc)
            fraction[c] = 1.0
            state[c] = 1
        elseif zigzag || all(>(0.0), φc)
            state[c] = 0
        else
            m = cellMoments(φc, maxDegree)
            f = m[1]
            fraction[c] = f
            if f <= tolerance
                state[c] = 0
            elseif f >= 1 - tolerance
                fraction[c] = 1.0
                state[c] = 1
            else
                state[c] = 2
                moments[c] = m
            end
        end
    end
    CutCellGeometry{N}(ls, maxDegree, fraction, state, moments)
end
