#!/usr/bin/env julia

using Distributions

function getGamma(;μ::Float64,β::Float64,F::Float64) #parameters as in Pancoast+ 2014 equation 7-10
    α = β^(-2)
    θ = μ*β^2*(1-F)
    return Gamma(α,θ) #note -- doesn't work for real things? appears too small by factor of μF i.e. need to add μF to result to "shift"?
end

function getR(rₛ::Float64, γ::Gamma{Float64}) #draw random number from gamma distribution
    return rₛ + rand(γ)
end

function getG(β::Float64) #equation 12
    return Gamma(1/β^2,1.0)
end

function getR(rₛ::Float64,μ::Float64,β::Float64,F::Float64,g::Float64) #equation 12, assumes g already set
    return rₛ+μ*F+μ*β^2*(1-F)*g
end

function getR(rₛ::Float64,μ::Float64,β::Float64,F::Float64,g::Gamma{Float64}) #equation 12, but draw from g
    return rₛ+μ*F+μ*β^2*(1-F)*rand(g)
end

function drawCloud(;μ::Float64=500.,β::Float64=1.0,F::Float64=0.5,ϕₒ::Float64=0.0,i::Float64=π/4,rot::Float64=0.0,rₛ::Float64=1.0,θₒ::Float64=0.0,θₒSystem::Float64=0.0,
    I::Union{Function,Float64}=IsotropicIntensity,v::Union{Function,Float64}=vCircularCloud,ξ::Float64=1.0, kwargs...)
    #assuming e = 0.0, vCircular, IsotropicIntensity for now
    # γ = getGamma(μ=μ,β=β,F=F)
    # r = getR(rₛ,γ)
    g = getG(β)
    r = getR(rₛ,μ,β,F,g)
    xyzSys = rotate3D(r,ϕₒ,i,rot,θₒ) #system coordinates xyz
    ϕ = atan(xyzSys[2],xyzSys[1]) # so that ϕ = 0 is along x
    # ϕ = acos(xyzSys[1]/r)
    reflect = (xyzSys[3] > midPlaneXZ(xyzSys[1],i)) && (rand() > ξ) #reflect particles from back of disk across disk midplane to front, leaving ξ fraction of particles in back
    if reflect
        xyzSys = reflect!(xyzSys,i)
        ϕ = atan(xyzSys[2],xyzSys[1])
    end
    return ring(r=r,i=i,v=v,I=I,ϕ=ϕ,ΔA=1.0,rot=rot,θₒ=θₒ,ϕₒ=ϕₒ,reflect=reflect;kwargs...)
end