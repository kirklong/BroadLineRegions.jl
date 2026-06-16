module BroadLineRegionsCUDAExt

using Adapt
using BroadLineRegions
using CUDA

const BLR = BroadLineRegions

BLR.gpu(ma::BLR.ModelArrays) = Adapt.adapt(CUDA.CuArray, ma)
BLR.gpu(m::BLR.model; T=Float32) = BLR.gpu(BLR.flatten(m; T=T))

BLR._rt_sortperm_by_key_depth(keys::CUDA.CuVector{Int}, x::CUDA.CuVector) = begin
    depth = ifelse.(isfinite.(x), x, convert(eltype(x), -Inf))
    depthOrder = sortperm(depth; rev=true)
    keyOrder = sortperm(keys[depthOrder])
    depthOrder[keyOrder]
end

BLR._rt_sortperm_by_key_depth(ma::BLR.ModelArrays{T,<:CUDA.CuVector{T}}, keys::CUDA.CuVector{Int}) where {T<:Real} =
    BLR._rt_sortperm_by_key_depth(keys, ma.x)

end
