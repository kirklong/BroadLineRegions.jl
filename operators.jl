#!/usr/bin/env julia

Base.:+(m1::model,m2::model) = begin
    r1 = m1.rings; r2 = m2.rings
    mCombined = deepcopy(m1)
    mCombined.rings = [r1...; r2...]
    mCombined.subModelStartInds = push!(mCombined.subModelStartInds, length(r1)+1)
    mCombined.profiles = Dict{Symbol,profile}()
    α1 = m1.camera.α; α2 = m2.camera.α
    β1 = m1.camera.β; β2 = m2.camera.β
    mCombined.camera = camera(vcat(α1,α2),vcat(β1,β2),nothing)
    return mCombined
end

# Base.getindex(m::model,i::Int)