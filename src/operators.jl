#!/usr/bin/env julia
"""
    Base.:+(m1::model, m2::model) combines two models by concatenating their rings and camera parameters.

Create a new model with the combined rings and camera parameters, and updates the subModelStartInds accordingly.
"""
Base.:+(m1::model,m2::model) = begin
    r1 = m1.rings; r2 = m2.rings
    mCombined = deepcopy(m1)
    mCombined.rings = [r1...; r2...]
    mCombined.subModelStartInds = push!(mCombined.subModelStartInds, length(r1)+1)
    mCombined.profiles = Dict{Symbol,profile}()
    α1 = m1.camera.α; α2 = m2.camera.α
    β1 = m1.camera.β; β2 = m2.camera.β
    raytraced = m1.camera.raytraced || m2.camera.raytraced
    mCombined.camera = camera(vcat(α1,α2),vcat(β1,β2),raytraced)
    return mCombined
end

"""
    Base.getindex(m::model, i::Int)

Retrieves the i-th submodel from the model `m`.
"""
Base.getindex(m::model, i::Int) = begin
    if i > length(m.subModelStartInds)
        error("Index out of bounds: $i for model with $(length(m.subModelStartInds)) submodels.")
    end
    startInd = m.subModelStartInds[i]
    endInd = i < length(m.subModelStartInds) ? m.subModelStartInds[i+1]-1 : length(m.rings)
    subModelRings = m.rings[startInd:endInd]
    camStartInds = getFlattenedCameraIndices(m)
    camStartInd = camStartInds[i]
    camEndInd = i < length(camStartInds) ? camStartInds[i+1]-1 : length(m.camera.α)
    subModelCamera = camera(m.camera.α[camStartInd:camEndInd],
                            m.camera.β[camStartInd:camEndInd],
                            m.camera.raytraced)
    return model(subModelRings, Dict{Symbol,profile}(), subModelCamera, [1]) 
end