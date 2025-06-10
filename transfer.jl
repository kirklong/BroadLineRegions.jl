#!/usr/bin/env julia
function response(r::Float64; ηₒ::Float64=0.5, η₁::Float64=0.5, αRM::Float64=0.0, rNorm::Float64=1.0, _...)
    η = ηₒ + η₁*(r/rNorm)^αRM
    return η
end

function getΨ(m::model,vEdges::Array{Float64},tEdges::Array{Float64};)
    I = getVariable(m,:I)
    ΔA = getVariable(m,:ΔA)
    v = getVariable(m,:v)
    d(ring::ring) = t(ring)
    delays = getVariable(m,d)
    Ψ = Array{Float64}(undef,length(vEdges)-1,length(tEdges)-1)
    for i in 1:length(vEdges)-1
        for j in 1:length(tEdges)-1
            mask = (v .>= vEdges[i]) .& (v .< vEdges[i+1]) .& (delays .>= tEdges[j]) .& (delays .< tEdges[j+1])
            s = sum(I[mask].*ΔA[mask])
            Ψ[i,j] = s > 0 ? s : 1e-30
        end
    end
    return Ψ
end

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

function getΨt(m::model,tEdges::Array{Float64},overflow::Bool=false;)
    I = getVariable(m,:I)
    ΔA = getVariable(m,:ΔA)
    d(ring::ring) = (typeof(ring.r) == Float64 && typeof(ring.ϕ) == Float64) ? tCloud(ring) : tDisk(ring)
    delays = getVariable(m,d)
    Ψt = Array{Float64}(undef,length(tEdges)-1)
    for j in 1:length(tEdges)-1
        mask = (delays .>= tEdges[j]) .& (delays .< tEdges[j+1])
        s = sum(I[mask].*ΔA[mask])
        Ψt[j] = s > 0 ? s : 1e-30
    end
    if overflow
        mask = delays .< tEdges[1]
        s = sum(I[mask].*ΔA[mask])
        Ψt[1] += s > 0 ? s : 1e-30
        mask = delays .>= tEdges[end]
        s = sum(I[mask].*ΔA[mask])
        Ψt[end] += s > 0 ? s : 1e-30
    end
    return Ψt
end

function getΨt(m::model,tBins::Int64,maxT::Float64=Inf,overflow::Bool=false)
    t(ring::ring) = (typeof(ring.r) == Float64 && typeof(ring.ϕ) == Float64) ? tCloud(ring) : tDisk(ring)
    delays = getVariable(m,t)
    if isinf(maxT)
        maxT =  maximum(i for i in delays if !isnan(i))
    end
    minT =  minimum(i for i in delays if !isnan(i))
    tEdges = collect(range(minT,stop=maxT,length=tBins+1))
    tCenters = @. (tEdges[1:end-1] + tEdges[2:end])/2
    return tCenters,getΨt(m,tEdges,overflow)
end