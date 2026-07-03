#!/usr/bin/env julia
"""
    Base.:+(m1::model, m2::model) combines two models by concatenating their rings and camera parameters.

Create a new model with the combined rings and camera parameters, and updates the subModelStartInds accordingly.
"""
Base.:+(m1::model,m2::model) = begin
    r1 = m1.rings; r2 = m2.rings
    mCombined = deepcopy(m1)
    mCombined.rings = [r1...; r2...]
    mCombined.subModelStartInds = vcat(m1.subModelStartInds, m2.subModelStartInds .+ length(r1))
    mCombined.profiles = Dict{Symbol,profile}()
    mCombined.cache = Dict{Any,Array}() #fresh cache -- combined rings differ from m1's
    α1 = m1.camera.α; α2 = m2.camera.α
    β1 = m1.camera.β; β2 = m2.camera.β
    raytraced = m1.camera.raytraced || m2.camera.raytraced
    mCombined.camera = camera(vcat(α1,α2),vcat(β1,β2),raytraced)
    #record construction provenance as a NESTED tree (not flattened) mirroring the actual + calls.
    #Set AFTER the deepcopy (which copied m1.params). Either side may be nothing (operand built
    #without a public constructor). See the model `params` field docstring and `rebuild`.
    mCombined.params = (; constructor=:+, left=m1.params, right=m2.params)
    return mCombined
end

#number of submodel slots a params (sub)tree spans, or nothing when it cannot be determined.
#All four public leaf constructors build single-submodel models ([1]); raytrace! compacts/regroups
#slots (and may append a free-clouds slot, see _rt_finalize_model), so a :raytrace! node's slot
#count is NOT derivable from its parent record.
_paramsSubmodelCount(::Nothing) = nothing
_paramsSubmodelCount(p::NamedTuple) = begin
    if p.constructor === :+
        l = _paramsSubmodelCount(p.left); r = _paramsSubmodelCount(p.right)
        (l === nothing || r === nothing) ? nothing : l + r
    elseif p.constructor === :raytrace!
        nothing
    else
        1
    end
end

#the params subtree describing submodel i of a params tree spanning n slots, or nothing when the
#extracted submodel has no faithful standalone construction record. When one side of a :+ has an
#unknown slot count the other side pins it (the node spans exactly n slots by construction); on any
#count mismatch with the actual slots, return nothing rather than guess.
_submodelParams(::Nothing, i::Int, n::Int) = nothing
_submodelParams(p::NamedTuple, i::Int, n::Int) = begin
    if p.constructor === :+
        l = _paramsSubmodelCount(p.left); r = _paramsSubmodelCount(p.right)
        (l === nothing && r === nothing) && return nothing
        l === nothing && (l = n - r)
        r === nothing && (r = n - l)
        (l >= 1 && r >= 1 && l + r == n) || return nothing
        return i <= l ? _submodelParams(p.left, i, l) : _submodelParams(p.right, i - l, r)
    end
    #non-:+ node (leaf constructor or :raytrace!) spanning this whole range: extracting the only
    #slot of a single-slot model IS the recorded model, so the record transfers; a submodel sliced
    #out of a multi-slot :raytrace! model has no standalone record (re-raytracing a slice alone is
    #NOT equivalent to slicing the raytraced whole -- occlusion couples the submodels)
    return n == 1 ? p : nothing
end

"""
    Base.getindex(m::model, i::Int)

Retrieves the i-th submodel from the model `m`.

The extracted submodel carries the matching side of a `:+` `params` provenance tree when that side
can be identified unambiguously (so e.g. `(a+b)[2].params == b.params` and the result is
[`rebuild`](@ref)-able); otherwise its `params` is `nothing` -- in particular for operands built
without a public constructor and for submodels sliced out of a multi-slot `raytrace!`d model
(re-raytracing a slice alone is not equivalent to slicing the raytraced whole, so no faithful
construction record exists).
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
    return model(subModelRings, Dict{Symbol,profile}(), subModelCamera, [1],
                 _submodelParams(m.params, i, length(m.subModelStartInds)))
end