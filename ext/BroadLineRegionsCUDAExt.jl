module BroadLineRegionsCUDAExt

# CUDA hooks for the backend-generic GPU entry points in core (`gpu`, `gpuDiskWindModel`,
# `gpuCloudModel`, `defaultGPUBackend` — see src/gpu_arrays.jl). Every method here dispatches on a
# CUDA-specific type (`Val{:CUDA}`, `CUDABackend`, `CuVector`), so this extension never shares a
# signature with another GPU extension (e.g. BroadLineRegionsMetalExt) and both can load together.

using Adapt
using BroadLineRegions
using CUDA

const BLR = BroadLineRegions

BLR._gpu_backend(::Val{:CUDA}) = CUDA.CUDABackend()
BLR._gpu_functional(::CUDA.CUDABackend) = CUDA.functional()
BLR._gpu_array_type(::CUDA.CUDABackend) = CUDA.CuArray

BLR._rt_backend_adapt(x, ma::BLR.ModelArrays{T,<:CUDA.CuVector{T}}) where {T<:Real} = Adapt.adapt(CUDA.CuArray, x)

# Order points the same way the CPU reference (`BLR._rt_sortperm_by_key_depth`) does:
# pixel key ascending, then depth `x` descending (front-to-back), then original index ascending.
# A single sort over the `(key, -depth, index)` tuple makes the order *total* -- the trailing index
# breaks every tie -- so the result bit-matches the CPU MergeSort reference and does NOT depend on
# whether CUDA's `sortperm` is stable. Non-finite depth sorts last within its key (matches the CPU
# `-Inf` sentinel). This also collapses the previous two passes into one sort.
BLR._rt_sortperm_by_key_depth(keys::CUDA.CuVector{Int}, x::CUDA.CuVector) = begin
    depth = ifelse.(isfinite.(x), x, convert(eltype(x), -Inf))
    idx = CUDA.CuArray(collect(1:length(keys)))
    sortperm(tuple.(keys, .-depth, idx))
end

BLR._rt_sortperm_by_key_depth(ma::BLR.ModelArrays{T,<:CUDA.CuVector{T}}, keys::CUDA.CuVector{Int}) where {T<:Real} =
    BLR._rt_sortperm_by_key_depth(keys, ma.x)

end
