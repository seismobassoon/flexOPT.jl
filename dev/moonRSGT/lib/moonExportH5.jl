# dev/moonRSGT/lib/moonExportH5.jl -- catalogue (JLD2 NamedTuples) -> compact HDF5 for Julia AND Python readers
#
# Compression, all in Julia (JLD2 filters, standard HDF5 -> h5py reads it without any plugin):
#   1. bit rounding (lossy, error-bounded): keep `keepbits` of the 23 float32 mantissa bits,
#      relative error ≤ 2^-(keepbits+1) per value (12 bits: 1.2e-4, 10 bits: 4.9e-4);
#      the dropped bits become zeros, which is what the next step compresses;
#   2. Shuffle + Deflate (lossless, built-in HDF5 filters; or Zstd = HDF5 filter 32015, smaller and
#      faster but Python then needs `import hdf5plugin`).
# The (r, θ) disk is split into θ blocks (one dataset each), so a reader decompresses only the
# block(s) holding the 6 distances it needs, not the whole 0.6 GB.
#
# Layout "moonRSGT-h5 v2" (shapes as seen from h5py / numpy, C order):
#   strain_re_<k>, strain_im_<k>   float32 (nω, 10, nr, nθ_k)   θ block k = 1 … nblocks (3 digits)
#   disp_re_<k>,   disp_im_<k>     float32 (nω,  5, nr, nθ_k)
#   theta_block_start (nblocks)    first θ index (1-based) of each block
#   pairs_*, omega, omega_i, tlen, r_km, layer, theta_deg, R_km, max_depth_km, model_id, keepbits, format
# (v1 = the same with single arrays strain_re, strain_im, disp_re, disp_im; still read by the readers.)

using JLD2, Printf

const MOON_H5_V2 = "moonRSGT-h5 v2"

"Round the mantissa of every Float32 to `keep` bits (round to nearest), in place."
function bitround!(A::AbstractArray{Float32}, keep::Integer)
    keep >= 23 && return A
    drop = 23 - keep
    half = UInt32(1) << (drop - 1)
    mask = ~((UInt32(1) << drop) - UInt32(1))
    @inbounds for i in eachindex(A)
        u = reinterpret(UInt32, A[i])
        A[i] = reinterpret(Float32, (u + half) & mask)
    end
    return A
end

"JLD2 filter object by name, wherever this JLD2 version defines it (JLD2 or JLD2.Filters); nothing if absent."
function _jld2filter(name::Symbol, args...)
    for m in (JLD2, isdefined(JLD2, :Filters) ? JLD2.Filters : nothing)
        m === nothing && continue
        isdefined(m, name) && return getfield(m, name)(args...)
    end
    return nothing
end

"`compress=` value for jldopen: :deflate (Shuffle + Deflate), :zstd, :none."
function moonH5Compressor(kind::Symbol)
    kind === :none && return false
    if kind === :deflate
        s = _jld2filter(:Shuffle); d = _jld2filter(:Deflate)
        (s !== nothing && d !== nothing) && return [s, d]
        @warn "this JLD2 has no Shuffle/Deflate filter objects: plain deflate (compress=true) instead"
        return true
    elseif kind === :zstd
        z = _jld2filter(:ZstdFilter)
        z === nothing && error("this JLD2 has no ZstdFilter (needs a recent JLD2 and CodecZstd)")
        s = _jld2filter(:Shuffle)
        return s === nothing ? z : [s, z]
    end
    error("compressor must be :deflate, :zstd or :none")
end

"(nω, np, nr, nθ) Julia array -> what h5py sees as (nω, np, nr, nθ) in C order."
_toC(A) = permutedims(A, (4, 3, 2, 1))

"""
    exportMoonH5(cat, out; keepbits=12, compressor=:deflate, θblock=16, withModel=false)

Write the compact HDF5 file of one loaded catalogue (`loadMoonRSGT(id)`), to `out*".part"`
then renamed.  Returns the file size in bytes.
"""
function exportMoonH5(cat, out::AbstractString; keepbits::Integer=12, compressor::Symbol=:deflate,
                      θblock::Integer=16, withModel::Bool=false)
    d = cat.disk
    nθ = length(d.θdeg)
    starts = collect(1:θblock:nθ)
    tmp = out * ".part"
    prep(A) = bitround!(Float32.(A), keepbits)
    jldopen(tmp, "w"; compress=moonH5Compressor(compressor)) do f
        f["format"] = MOON_H5_V2
        f["model_id"] = string(replace(basename(cat.file), r"^rsgt_|\.jld2$" => ""))
        f["keepbits"] = Int(keepbits)
        f["theta_block_start"] = starts
        for (k, s) in enumerate(starts)
            idx = s:min(s + θblock - 1, nθ)
            tag = lpad(k, 3, '0')
            f["strain_re_" * tag] = _toC(prep(real.(view(d.strain, :, :, :, idx))))
            f["strain_im_" * tag] = _toC(prep(imag.(view(d.strain, :, :, :, idx))))
            f["disp_re_" * tag] = _toC(prep(real.(view(d.disp, :, :, :, idx))))
            f["disp_im_" * tag] = _toC(prep(imag.(view(d.disp, :, :, :, idx))))
        end
        f["pairs_strain_src"] = [p[1] for p in d.pairs.strain]; f["pairs_strain_comp"] = [p[2] for p in d.pairs.strain]
        f["pairs_disp_src"] = [p[1] for p in d.pairs.disp];     f["pairs_disp_comp"] = [p[2] for p in d.pairs.disp]
        f["omega"] = collect(Float64, d.ω); f["omega_i"] = Float64(d.ωI); f["tlen"] = Float64(d.tlen)
        f["r_km"] = collect(Float64, d.r_km); f["layer"] = collect(Int, d.seg); f["theta_deg"] = collect(Float64, d.θdeg)
        f["R_km"] = Float64(cat.info.R); f["max_depth_km"] = Float64(cat.maxDepth)
        if withModel
            f["model_depth"] = cat.model.depth; f["model_vp"] = cat.model.vp
            f["model_vs"] = cat.model.vs; f["model_rho"] = cat.model.ρ
        end
    end
    mv(tmp, out; force=true)
    return filesize(out)
end
