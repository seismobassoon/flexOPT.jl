# dev/freeSurfacePlaneWave/codeVersions.jl
CODE_VERSION = Dict(
    "opt" => "v3-2026-11-15-pin-lateral-edges",  # v2 + pinning the two x-
        # edge columns to the exact analytical solution every step (v2
        # left them frozen at their t=0 IC snapshot, since nothing ever
        # updates them otherwise -- harmless when the field starts/stays
        # at rest, wrong once the IC is nonzero everywhere)
    "fd"  => "v3-2026-11-15-pin-lateral-edges",  # same fix, FD3 side
)
PLANEWAVE_CACHE_DIR = joinpath("dev", "freeSurfacePlaneWave")
