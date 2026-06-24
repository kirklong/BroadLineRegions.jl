#!/usr/bin/env julia
"""
    response(r::Float64; ηₒ::Float64=0.5, η₁::Float64=0.5, αRM::Float64=0.0, rNorm::Float64=1.0, _...)

Calculate response function for use in reverberation mapping calculations.
"""
function response(r::Float64; ηₒ::Float64=0.5, η₁::Float64=0.5, αRM::Float64=0.0, rNorm::Float64=1.0, _...)
    η = ηₒ + η₁*(r/rNorm)^αRM
    return η
end

"""
    response(r::Float64, ηₒ::Float64, η₁::Float64, αRM::Float64, rNorm::Float64)

Positional-argument version of `response` for hot loops — avoids the per-call keyword splat.
"""
response(r::Float64, ηₒ::Float64, η₁::Float64, αRM::Float64, rNorm::Float64) = ηₒ + η₁*(r/rNorm)^αRM

"""
    getΨ(m::model,vEdges::Array{Float64},tEdges::Array{Float64})

Calculate the 2D transfer function Ψ for a model `m` over specified velocity and time bins, whose edges are given by `vEdges` and `tEdges`.

A device-resident method (`getΨ(::ResidentModel, …)`) runs this on the GPU without re-flattening the
model each call; see [`gpu`](@ref BLR.gpu) and [`ResidentModel`](@ref BLR.ResidentModel).
"""
function getΨ(m::model,vEdges::Array{Float64},tEdges::Array{Float64})
    I = getVariable(m,:I)
    ΔA = getVariable(m,:ΔA)
    v = getVariable(m,:v)
    delays = getVariable(m,t) #pass t itself (not a closure) so the result is memoized in m.cache
    Ψ = zeros(length(vEdges)-1,length(tEdges)-1)
    ΨAccumulate!(Ψ, v, delays, I, ΔA, vec(vEdges), vec(tEdges))
    @. Ψ = ifelse(Ψ > 0, Ψ, 1e-30) #empty (or NaN-poisoned) bins floored, exactly as the old per-bin `s > 0 ? s : 1e-30`
    return Ψ
end

"""
    ΨAccumulate!(Ψ, v, delays, I, ΔA, vEdges, tEdges)

Single-pass accumulation for `getΨ` (function barrier): one sweep over all points, locating each
point's (velocity, delay) bin instead of building a boolean mask per bin pair (the old approach was
O(nbins² × npoints) in time and memory traffic).

Bin assignment uses the **same convention as `binnedSum`** (`searchsortedfirst`): interior edges are
left-EXCLUSIVE / right-inclusive — a point exactly on an interior edge goes to the bin below it — and
points at or below the first edge or at or above the last edge are dropped (Ψ has no overflow bins).
NaN velocities/delays never land in a bin; NaN intensities/areas still accumulate into their bin so
the final floor maps the poisoned bin to 1e-30. This matches the line profile and every other
`binnedSum`-based observable, so a shared edge set bins identically everywhere.
"""
function ΨAccumulate!(Ψ::Matrix{Float64}, v, delays, I, ΔA, vEdges::AbstractVector{Float64}, tEdges::AbstractVector{Float64})
    nV = length(vEdges)-1; nT = length(tEdges)-1
    @inbounds for (vi, di, Ii, Ai) in zip(v, delays, I, ΔA)
        (isfinite(vi) && isfinite(di)) || continue
        (vi <= vEdges[1] || vi >= vEdges[end]) && continue   # binnedSum convention, no overflow bins
        (di <= tEdges[1] || di >= tEdges[end]) && continue
        bv = searchsortedfirst(vEdges, vi) - 1
        bt = searchsortedfirst(tEdges, di) - 1
        Ψ[bv, bt] += Ii*Ai
    end
    return Ψ
end

"""
    getΨ(m::model,vBins::Int64,tBins::Int64)

Calculate the 2D transfer function Ψ for a model `m` over specified number of velocity bins `vBins` and time bins `tBins`.
The velocity and time edges are automatically calculated based on the minimum and maximum values for velocity and delays in the model.
"""
function getΨ(m::model,vBins::Int64,tBins::Int64)
    v = getVariable(m,:v)
    delays = getVariable(m,t)
    maxV =  maximum(i for i in v if !isnan(i))
    minV =  minimum(i for i in v if !isnan(i))
    maxT =  maximum(i for i in delays if !isnan(i))
    minT =  minimum(i for i in delays if !isnan(i))
    vEdges = collect(range(minV,stop=maxV,length=vBins+1))
    tEdges = collect(range(minT,stop=maxT,length=tBins+1))
    vCenters = @. (vEdges[1:end-1] + vEdges[2:end])/2
    tCenters = @. (tEdges[1:end-1] + tEdges[2:end])/2
    return vCenters,tCenters,getΨ(m,vEdges,tEdges)
end

"""
    getΨt(m::model,tEdges::Array{Float64},overflow::Bool=false;)

Calculate the 1D transfer function Ψ(t) for a model `m` over specified time edges `tEdges`.
The `overflow` parameter determines whether to include contributions from delays outside the specified edges in the edge bins.

A device-resident method (`getΨt(::ResidentModel, …)`) runs this on the GPU without re-flattening the
model each call; see [`gpu`](@ref BLR.gpu) and [`ResidentModel`](@ref BLR.ResidentModel).
"""
function getΨt(m::model,tEdges::Array{Float64},overflow::Bool=false;)
    I = getVariable(m,:I)
    ΔA = getVariable(m,:ΔA)
    delays = getVariable(m,t)
    Ψt = zeros(length(tEdges)-1)
    sUnder, sOver = ΨtAccumulate!(Ψt, delays, I, ΔA, vec(tEdges))
    @. Ψt = ifelse(Ψt > 0, Ψt, 1e-30) #empty (or NaN-poisoned) bins floored BEFORE overflow is added, as before
    if overflow
        Ψt[1] += sUnder > 0 ? sUnder : 1e-30 #note: adds 1e-30 even when nothing under/overflowed -- pre-existing behavior, kept
        Ψt[end] += sOver > 0 ? sOver : 1e-30
    end
    return Ψt
end

"""
    ΨtAccumulate!(Ψt, delays, I, ΔA, tEdges) -> (sUnder, sOver)

Single-pass accumulation for `getΨt` (function barrier -- the gathered arrays are not inferrable in
the caller). Same `binnedSum` (`searchsortedfirst`, left-exclusive interior edge) convention as
`ΨAccumulate!`; additionally collects the underflow/overflow sums (`delays <= tEdges[1]` /
`delays >= tEdges[end]`, NaN in neither) in the same sweep, matching `binnedSum`'s `overflow=true`
edge-bin folding.
"""
function ΨtAccumulate!(Ψt::Vector{Float64}, delays, I, ΔA, tEdges::AbstractVector{Float64})
    nT = length(tEdges)-1
    sUnder = 0.0; sOver = 0.0
    @inbounds for (di, Ii, Ai) in zip(delays, I, ΔA)
        isnan(di) && continue #NaN delays belong to neither a bin nor the overflow buckets (as before)
        if di <= tEdges[1]
            sUnder += Ii*Ai
        elseif di >= tEdges[end]
            sOver += Ii*Ai
        else
            Ψt[searchsortedfirst(tEdges, di) - 1] += Ii*Ai
        end
    end
    return sUnder, sOver
end

"""
    getΨt(m::model,tBins::Int64,maxT::Float64=Inf,overflow::Bool=false)

Calculate the 1D transfer function Ψ(t) for a model `m` over specified number of time bins `tBins`.
The `maxT` parameter specifies the maximum time delay to consider, and `overflow` determines whether to include contributions from delays outside the specified edges in the edge bins.
"""
function getΨt(m::model,tBins::Int64,maxT::Float64=Inf,overflow::Bool=false)
    delays = getVariable(m,t)
    if isinf(maxT)
        maxT =  maximum(i for i in delays if !isnan(i))
    end
    minT =  minimum(i for i in delays if !isnan(i))
    tEdges = collect(range(minT,stop=maxT,length=tBins+1))
    tCenters = @. (tEdges[1:end-1] + tEdges[2:end])/2
    return tCenters,getΨt(m,tEdges,overflow)
end