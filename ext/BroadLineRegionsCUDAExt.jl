module BroadLineRegionsCUDAExt

using Adapt
using BroadLineRegions
using CUDA

const BLR = BroadLineRegions

BLR.gpu(ma::BLR.ModelArrays) = Adapt.adapt(CUDA.CuArray, ma)
BLR.gpu(m::BLR.model; T=Float32) = BLR.gpu(BLR.flatten(m; T=T))

end
