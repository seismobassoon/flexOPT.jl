# dev/new_Paul_models/codeVersions.jl -- one string per cached stage (enters
# its myProduceOrLoad key). Changing the CODE of a stage is not seen by the
# hash: bump the string, or run once with --redo=<stage>.
CODE_VERSION = Dict(
    "fd3" => "v1-2026-11-15",     # lib/paulSeismo.jl: model, moment tensor, FD3 stepping, receivers
    "taup" => "v1-2026-11-15",
    "nu" => "v1-2026-11-16",      # lib/paulNeutrino.jl: nₑ, rotated detector frames, Pνν, events    # PREM τ–p first-arrival table per phase (dev/raytheory)
)
