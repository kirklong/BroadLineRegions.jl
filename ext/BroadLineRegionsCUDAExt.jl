module BroadLineRegionsCUDAExt

using Adapt
using BroadLineRegions
using CUDA

const BLR = BroadLineRegions

BLR.gpu(ma::BLR.ModelArrays) = Adapt.adapt(CUDA.CuArray, ma)
# User entry point: a device-resident handle (ModelArrays on the GPU + the CUDA backend) that the
# observable methods (getProfile/getΨ/getΨt/phase/secondMoment) reuse across calls.
BLR.gpu(m::BLR.model; T=Float32) =
    BLR.ResidentModel(BLR.gpu(BLR.flatten(m; T=T)), CUDA.CUDABackend(), length(m.subModelStartInds))
# Bare device ModelArrays for the raytrace backend (no wrapper -- the scan needs the columns directly).
BLR._rt_backend_model_arrays(m::BLR.model, ::CUDA.CUDABackend; T=Float64) = BLR.gpu(BLR.flatten(m; T=T))
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
