# src/cutCells/cutNodes.jl
#
# Which nodes need a cut-cell operator, the cells of their test-function
# support, and the crossings of the interface with the segments from a node
# to its stencil neighbours (Ovcharenko 2015, Algorithm 1: η0, η1 and the
# normal at the crossing -- the input of the interface transfer, stage 3).

"Nodes on the \"−\" side (φ <= 0): the solid for a free surface."
sideMask(ls::CutLevelSet; side::Symbol=:minus) =
    side === :minus ? ls.φ .<= 0.0 : ls.φ .> 0.0

"""
    cutNodes(geom; side=:minus)

Nodes whose support (the 2^N cells touching the node) holds material of side
`side` AND a cut cell or a cell entirely on the other side. Those are the rows
whose volume integral differs from the uncut operator. Cells outside the grid
are ignored (artificial edges stay as they are).

The list can contain nodes on the other side (φ > 0 for `side=:minus`) whose
support still reaches a sliver of material: the usual small-cut-cell
question (keep them as extended unknowns, or merge them into a neighbour).
That choice belongs to the operators (stage 2), so they are not removed here;
`sideMask` tells which is which.
"""
function cutNodes(geom::CutCellGeometry{N}; side::Symbol=:minus) where {N}
    own = side === :minus ? Int8(1) : Int8(0)
    nodes = CartesianIndex{N}[]
    for I in CartesianIndices(geom.levelSet.φ)
        touchesOwn = false
        touchesOther = false
        for k in CartesianIndices(ntuple(_ -> -1:0, N))
            c = I + k
            checkbounds(Bool, geom.state, c) || continue
            st = geom.state[c]
            if st == 2
                touchesOwn = true
                touchesOther = true
            elseif st == own
                touchesOwn = true
            else
                touchesOther = true
            end
        end
        touchesOwn && touchesOther && push!(nodes, I)
    end
    return nodes
end

"""
    nodeSupportCells(geom, I; side=:minus)

The cells around node `I`, as NamedTuples
`(cell, offset, state, moments)`: `offset ∈ {-1,0}^N` is the cell's lower
corner relative to `I`; `moments` are the side-`side` moments in the cell's
local coordinates (`nothing` for a cell with no material on that side, the
full-cell moments for an uncut cell of that side).
"""
function nodeSupportCells(geom::CutCellGeometry{N}, I::CartesianIndex{N};
                          side::Symbol=:minus) where {N}
    own = side === :minus ? Int8(1) : Int8(0)
    full = fullCellMoments(N, geom.maxDegree)
    out = NamedTuple[]
    for k in CartesianIndices(ntuple(_ -> -1:0, N))
        c = I + k
        checkbounds(Bool, geom.state, c) || continue
        st = geom.state[c]
        m = if st == 2
            side === :minus ? geom.moments[c] : full .- geom.moments[c]
        elseif st == own
            full
        else
            nothing
        end
        push!(out, (cell=c, offset=Tuple(k), state=st, moments=m))
    end
    return out
end

"""
    stencilCrossings(ls, I; reach=1)

For every stencil neighbour `J = I + o`, `o ∈ (-reach:reach)^N \\ 0`, on the
other side of the interface: the crossing of segment I→J, as
`(offset, η0, η1, ηAxes, point, normal)`. `η0` is the fraction of the segment
from `I` (φ linear along it), `η1 = 1 - η0`, `ηAxes = η0 .* o` the
per-axis fractions of the grid spacing (η0x, η0z of the thesis), `normal` the
unit normal (− → +) at the crossing.
"""
function stencilCrossings(ls::CutLevelSet{N}, I::CartesianIndex{N}; reach::Int=1) where {N}
    out = NamedTuple[]
    φ0 = ls.φ[I]
    xI = nodeCoordinates(ls, I)
    for o in CartesianIndices(ntuple(_ -> -reach:reach, N))
        all(iszero, Tuple(o)) && continue
        J = I + o
        checkbounds(Bool, ls.φ, J) || continue
        φ1 = ls.φ[J]
        ((φ0 <= 0) == (φ1 <= 0)) && continue
        η0 = φ0 / (φ0 - φ1)
        point = xI + η0 * (nodeCoordinates(ls, J) - xI)
        push!(out, (offset=Tuple(o), η0=η0, η1=1 - η0,
            ηAxes=η0 .* Tuple(o), point=point, normal=levelSetNormal(ls, point)))
    end
    return out
end
