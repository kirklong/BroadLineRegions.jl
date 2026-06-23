module BroadLineRegionsCUDAExt

using Adapt
using BroadLineRegions
using CUDA

const BLR = BroadLineRegions

BLR.gpu(ma::BLR.ModelArrays) = Adapt.adapt(CUDA.CuArray, ma)
# User entry point: a device-resident handle (ModelArrays on the GPU + the CUDA backend) that the
# observable methods (getProfile/getΨ/getΨt/phase/secondMoment) reuse across calls. For combined models
# (>1 submodel) it also carries device-resident raytrace metadata (output grid + per-point submodel /
# discrete info) built once here on the host, so `raytrace!(rm)` runs entirely on the GPU.
function BLR.gpu(m::BLR.model; T=Float32)
    ma = BLR.gpu(BLR.flatten(m; T=T))
    meta = length(m.subModelStartInds) > 1 ? Adapt.adapt(CUDA.CuArray, BLR._rt_build_meta(m; T=T)) : nothing
    return BLR.ResidentModel(ma, CUDA.CUDABackend(), length(m.subModelStartInds), meta)
end
# Bare device ModelArrays for the raytrace backend (no wrapper -- the scan needs the columns directly).
BLR._rt_backend_model_arrays(m::BLR.model, ::CUDA.CUDABackend; T=Float64) = BLR.gpu(BLR.flatten(m; T=T))

# On-device DiskWind construction: build the ModelArrays columns directly on the GPU (one fused
# kernel, no host rings), returning a device-resident handle. Float32 by default (GeForce FP64 is
# ~1/64 rate). Mirrors `residentDiskWindModel` with `backend=CUDABackend()`.
BLR.gpuDiskWindModel(rMin::Real, rMax::Real, i::Real; T=Float32, kwargs...) =
    BLR.residentDiskWindModel(rMin, rMax, i; backend=CUDA.CUDABackend(), T=T, kwargs...)
BLR.gpuDiskWindModel(r̄::Real, rFac::Real, Sα::Real, i::Real; T=Float32, kwargs...) =
    BLR.residentDiskWindModel(r̄, rFac, Sα, i; backend=CUDA.CUDABackend(), T=T, kwargs...)

# On-device cloud construction (counter-based Philox substreams, no host rings). Float32 default.
BLR.gpuCloudModel(nClouds::Int, seed::Integer; T=Float32, kwargs...) =
    BLR.residentCloudModel(nClouds, seed; backend=CUDA.CUDABackend(), T=T, kwargs...)
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
