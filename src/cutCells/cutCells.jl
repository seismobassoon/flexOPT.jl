# src/cutCells/cutCells.jl -- cut-cell volume integrals for flexOPT
#
# Free surfaces and internal discontinuities that cut the Cartesian grid
# (1-D, 2-D, 3-D). See README.md in this folder for the plan:
#
#   stage 1 (this code): interface description (level set, from a function, a
#            height or GeoPoints points + normals), cell classification,
#            exact moments of the cut part of each cell, stencil crossings
#            and normals;
#   stage 2: OPT rows whose test-function integral covers only the material
#            part of the support (non-separable contraction with the moments);
#   stage 3: internal discontinuities -- each side with its own material,
#            trial interpolation across the interface through the jump
#            conditions (Mizutani 2002, Ovcharenko 2015 eq. 37-38).
#
# Included by src/flexOPT.jl inside module flexOPT.

using StaticArrays
using LinearAlgebra
using SparseArrays

include("simplexIntegration.jl")
include("levelSets.jl")
include("cellMoments.jl")
include("cutNodes.jl")
include("aggregation.jl")
# the operator part needs the rest of flexOPT (recipes, @unpack, @strdict);
# dev/cutCells/run.jl geometry loads the geometry files alone.
(@isdefined(CUTCELLS_GEOMETRY_ONLY) && CUTCELLS_GEOMETRY_ONLY) || include("cutCellOperators.jl")

export CutLevelSet, levelSetFromFunction, levelSetFromHeight
export levelSetFromPointsNormals, interfaceFromLocalCoords
export levelSetValue, levelSetNormal, nodeCoordinates, sideMask
export grundmannMollerRule, simplexVolume, addSimplexMonomialMoments!
export kuhnSimplices, clippedSimplexPieces, cellMoments, fullCellMoments
export CutCellGeometry, zigzagLevelSet, cutNodes, nodeSupportCells, stencilCrossings
export sliverAggregation, applyAggregation
export CutCellScheme, cutCellRows, cutCellSelfCheck, installRows, pinRows
