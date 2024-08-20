#!/usr/bin/env julia

using Distributions, Random

function getGamma(;μ::Float64,β::Float64,F::Float64) #parameters as in Pancoast+ 2014 equation 7-10
    α = β^(-2)
    θ = μ*β^2*(1-F)
    return Gamma(α,θ) #note -- doesn't work for real things? appears too small by factor of μF i.e. need to add μF to result to "shift"?
end

function getR(rₛ::Float64, γ::Gamma{Float64},rng::AbstractRNG=Random.GLOBAL_RNG) #draw random number from gamma distribution
    return rₛ + rand(rng,γ)
end

function getG(β::Float64) #equation 12
    return Gamma(1/β^2,1.0)
end

function getR(rₛ::Float64,μ::Float64,β::Float64,F::Float64,g::Float64) #equation 12, assumes g already set
    return rₛ+μ*F+μ*β^2*(1-F)*g
end

function getR(rₛ::Float64,μ::Float64,β::Float64,F::Float64,g::Gamma{Float64},rng::AbstractRNG=Random.GLOBAL_RNG) #equation 12, but draw from g
    return rₛ+μ*F+μ*β^2*(1-F)*rand(rng,g)
end

function drawCloud(;μ::Float64=500.,β::Float64=1.0,F::Float64=0.5,ϕₒ::Float64=0.0,i::Float64=π/4,rot::Float64=0.0,rₛ::Float64=1.0,θₒ::Float64=0.0,θₒSystem::Float64=0.0,
    I::Union{Function,Float64}=IsotropicIntensity,v::Union{Function,Float64}=vCircularCloud,ξ::Float64=1.0,rng::AbstractRNG=Random.GLOBAL_RNG,kwargs...)
    #assuming e = 0.0, vCircular, IsotropicIntensity for now
    # γ = getGamma(μ=μ,β=β,F=F)
    # r = getR(rₛ,γ)
    g = getG(β)
    r = getR(rₛ,μ,β,F,g,rng)
    xyzSys = rotate3D(r,ϕₒ,i,rot,θₒ) #system coordinates xyz, #flip i to match convention of +z being up, relic
    reflect = (xyzSys[3] > midPlaneXZ(xyzSys[1],i)) && (rand(rng) > ξ) #reflect particles from back of disk across disk midplane to front, leaving ξ fraction of particles in back
    if reflect
        xyzSys = reflect!(xyzSys,i)
    end
    undo_tilt = [sin(i) 0.0 cos(i); 0.0 1.0 0.0; -cos(i) 0.0 sin(i)]
    xyzPlane = undo_tilt*xyzSys
    ϕ = atan(xyzPlane[2],xyzPlane[1]) #ϕ after rotation, measured from +x in disk plane
    return ring(r=r,i=i,v=v,I=I,ϕ=ϕ,ΔA=1.0,rot=rot,θₒ=θₒ,ϕₒ=ϕₒ,reflect=reflect,rng=rng;kwargs...)
end