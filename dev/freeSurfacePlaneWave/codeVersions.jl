# dev/freeSurfacePlaneWave/codeVersions.jl
CODE_VERSION = Dict(
    "opt" => "v2-2026-11-15-planewave-ic",  # recipes + operators + OPT3/OPT5
        # :direct propagation from an exact plane-wave INITIAL CONDITION,
        # zero forcing throughout (replaces the old sheet-source/Green's-
        # function-forcing v1)
    "fd"  => "v2-2026-11-15-planewave-ic",  # explicit elasticWave2D FD3,
        # same plane-wave-IC, zero-forcing replacement of v1
)
PLANEWAVE_CACHE_DIR = joinpath("dev", "freeSurfacePlaneWave")
