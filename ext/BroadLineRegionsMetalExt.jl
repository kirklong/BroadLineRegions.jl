module BroadLineRegionsMetalExt

# Apple-silicon (Metal) hooks for the backend-generic GPU entry points in core (`gpu`,
# `gpuDiskWindModel`, `gpuCloudModel`, `defaultGPUBackend` — see src/gpu_arrays.jl). All kernels are
# the shared KernelAbstractions ones; this file only registers the backend and supplies the pieces
# that differ on Metal. Every method dispatches on a Metal-specific type, so it never clashes with
# BroadLineRegionsCUDAExt.
#
# Metal differences handled here:
#   * no Float64 on Apple GPUs -> `_check_gpu_eltype` rejects T=Float64 with an actionable error;
#   * Metal.jl cannot `sortperm` a vector of tuples -> the raytrace (key, depth, index) sort runs the
#     CPU reference sort on a host copy and uploads the permutation (bit-identical order to the CPU
#     and CUDA paths; on unified memory the copy is cheap, and only raytrace! uses it).

using Adapt
using BroadLineRegions
using Metal

const BLR = BroadLineRegions

BLR._gpu_backend(::Val{:Metal}) = Metal.MetalBackend()
BLR._gpu_functional(::Metal.MetalBackend) = Metal.functional()
BLR._gpu_array_type(::Metal.MetalBackend) = Metal.MtlArray

function BLR._check_gpu_eltype(::Metal.MetalBackend, ::Type{T}) where {T}
    T === Float64 && throw(ArgumentError(
        "Metal (Apple GPU) does not support Float64 — use T=Float32 (the default for gpu, " *
        "gpuDiskWindModel and gpuCloudModel; pass T=Float32 to raytrace!/residentDiskWindModel/" *
        "residentCloudModel when backend=MetalBackend())"))
    return nothing
end

BLR._rt_backend_adapt(x, ma::BLR.ModelArrays{T,<:Metal.MtlVector{T}}) where {T<:Real} = Adapt.adapt(Metal.MtlArray, x)

# Host sort + upload: reuse the CPU reference (`_rt_sortperm_by_key_depth(::AbstractVector{Int}, ...)`,
# stable MergeSort by depth desc then key asc == total order (key, -depth, index)), so the permutation
# bit-matches the CPU backend and the CUDA extension's tuple sort.
function BLR._rt_sortperm_by_key_depth(keys::Metal.MtlVector{Int}, x::Metal.MtlVector)
    perm = BLR._rt_sortperm_by_key_depth(Array(keys), Array(x))
    return Metal.MtlArray(perm)
end

BLR._rt_sortperm_by_key_depth(ma::BLR.ModelArrays{T,<:Metal.MtlVector{T}}, keys::Metal.MtlVector{Int}) where {T<:Real} =
    BLR._rt_sortperm_by_key_depth(keys, ma.x)

end
