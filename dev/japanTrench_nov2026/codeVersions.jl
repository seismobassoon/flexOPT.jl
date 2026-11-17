# dev/japanTrench_nov2026/codeVersions.jl
# One version string per cached stage (enters its myProduceOrLoad key).
# Changing a parameter already gives a new cache file; changing the CODE of a
# stage (or src/ behind it) is not seen -> bump the string, or --redo=<stage>.
CODE_VERSION = Dict(
    "model" => "v1-2026-11-16",   # getParamsAndTopo: NIED + planet1D fallback, earth_relief, water, air
    "opt"   => "v1-2026-11-17",   # OPT3 3-D, Cerjan sides+bottom, free surface on rock top, ocean load, plane-wave sheet
)
