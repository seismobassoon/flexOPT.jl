# dev/spherical1DelasticKernel/codeVersions.jl -- one version string per cached stage
# (bump when the code behind it changes: lib/*.jl here or dev/spherical1Delastic/lib).
CODE_VERSION_K = Dict(
    "k_systems" => "v1-2026-10-28",   # elasticGlobalSystem PREM / layered (dev/spherical1Delastic/lib)
    "k1_disks"  => "v1-2026-10-28",   # computeSGTDisk for the checks
    "k2_disks"  => "v1-2026-10-28",   # computeSGTDisk PREM, grid radii
    "k3_born"   => "v1-2026-10-28",   # bornBruteForceCheck
)
KER_CACHE_DIR = joinpath("dev", "spherical1DelasticKernel")
