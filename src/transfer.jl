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
point's (velocity, delay) bin with `searchsortedlast`, instead of building a boolean mask per bin
pair (the old approach was O(nbins² × npoints) in time and memory traffic).

Bin assignment is identical to the original `>= left, < right` comparisons: interior edges are
left-inclusive, values at the right-most edge are excluded, NaN velocities/delays never land in a
bin (`searchsortedlast` sorts them past the end), and `+0.0` normalization keeps `-0.0` on the same
side as the old IEEE `>=` comparisons. NaN intensities/areas accumulate into their bin so the final
floor maps the poisoned bin to 1e-30, exactly as before.
"""
function ΨAccumulate!(Ψ::Matrix{Float64}, v, delays, I, ΔA, vEdges::AbstractVector{Float64}, tEdges::AbstractVector{Float64})
    nV = length(vEdges)-1; nT = length(tEdges)-1
    @inbounds for (vi, di, Ii, Ai) in zip(v, delays, I, ΔA)
        bv = searchsortedlast(vEdges, vi + 0.0)
        (1 <= bv <= nV) || continue
        bt = searchsortedlast(tEdges, di + 0.0)
        (1 <= bt <= nT) || continue
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
the caller). Same conventions as `ΨAccumulate!`; additionally collects the underflow/overflow sums
(`delays < tEdges[1]` / `delays >= tEdges[end]`, NaN in neither) in the same sweep.
"""
function ΨtAccumulate!(Ψt::Vector{Float64}, delays, I, ΔA, tEdges::AbstractVector{Float64})
    nT = length(tEdges)-1
    sUnder = 0.0; sOver = 0.0
    @inbounds for (di, Ii, Ai) in zip(delays, I, ΔA)
        isnan(di) && continue #NaN delays belong to neither a bin nor the overflow buckets (as before)
        bt = searchsortedlast(tEdges, di + 0.0) #+0.0: -0.0 -> +0.0 so isless-based search matches the old >=/< comparisons
        if bt < 1
            sUnder += Ii*Ai
        elseif bt > nT
            sOver += Ii*Ai
        else
            Ψt[bt] += Ii*Ai
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