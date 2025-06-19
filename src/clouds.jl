#!/usr/bin/env julia

using Distributions, Random

"""
    getGamma(;μ::Float64,β::Float64,F::Float64)

Returns a Gamma distribution with parameters derived from the input parameters as in Pancoast+ 2014 equations 7-10.
"""
function getGamma(;μ::Float64,β::Float64,F::Float64) #parameters as in Pancoast+ 2014 equation 7-10
    α = β^(-2)
    θ = μ*β^2*(1-F)
    return Gamma(α,θ) #note -- doesn't work for real things? appears too small by factor of μF i.e. need to add μF to result to "shift"?
end
"""
    getR(rₛ::Float64,γ::Gamma{Float64},rng::AbstractRNG=Random.GLOBAL_RNG)

Returns a random number drawn from the Gamma distribution `γ` and shifted by `rₛ`.
"""
function getR(rₛ::Float64, γ::Gamma{Float64},rng::AbstractRNG=Random.GLOBAL_RNG) #draw random number from gamma distribution
    return rₛ + rand(rng,γ)
end

"""
    getG(β::Float64)

Returns a Gamma distribution with parameters derived from the input parameter `β` as in Pancoast+ 2014 equation 12.
"""
function getG(β::Float64) #equation 12
    return Gamma(1/β^2,1.0)
end

"""
    getR(rₛ::Float64,μ::Float64,β::Float64,F::Float64,g::Float64)

Returns the radius `r` calculated using Pancoast+ 2014 equation 12, where `g` is already set.
"""
function getR(rₛ::Float64,μ::Float64,β::Float64,F::Float64,g::Float64) #equation 12, assumes g already set
    return rₛ+μ*F+μ*β^2*(1-F)*g
end

"""
    getR(rₛ::Float64,μ::Float64,β::Float64,F::Float64,
        g::Gamma{Float64},rng::AbstractRNG=Random.GLOBAL_RNG)

Returns the radius `r` calculated using Pancoast+ 2014 equation 12, where `g` is a Gamma distribution and a random number is drawn from it.
"""
function getR(rₛ::Float64,μ::Float64,β::Float64,F::Float64,g::Gamma{Float64},rng::AbstractRNG=Random.GLOBAL_RNG) #equation 12, but draw from g
    return rₛ+μ*F+μ*β^2*(1-F)*rand(rng,g)
end

"""
    drawCloud(;μ::Float64=500.,β::Float64=1.0,F::Float64=0.5,
            ϕ₀::Float64=0.0,i::Float64=π/4,rot::Float64=0.0,rₛ::Float64=1.0,
            θₒ::Float64=0.0,θₒSystem::Float64=0.0,I::Union{Function,Float64}=IsotropicIntensity,
            v::Union{Function,Float64}=vCircularCloud,ξ::Float64=1.0,
            rng::AbstractRNG=Random.GLOBAL_RNG,kwargs...)

Generates a model `ring` struct for a single cloud drawn from a thick disk-like structure with parameters defined by the input arguments similar to Pancoast+ 2011 and 2014.

# Arguments: 
- `μ::Float64=500.`: Mean radius of model (in terms of ``r_s``)
- `β::Float64=1.0`: Shape parameter for radial distribution
- `F::Float64=0.5`: Beginning radius in units of μ where clouds can be placed. 
- `ϕ₀::Float64=0.0`: Initial azimuthal angle of cloud (rad) 
- `i::Float64=π/4`: Inclination angle (rad) 
- `rot::Float64=0.0`: Random rotation of cloud about z axis (rad) 
- `rₛ::Float64=1.0`: Scale radius (in terms of ``r_s``)
- `θₒ::Float64=0.0`: Opening angle of cloud (rad) 
- `θₒSystem::Float64=0.0`: Maximum opening angle of the system (rad) 
- `ξ::Float64=1.0`: Fraction of clouds in back side that have not been moved to the front (when ξ = 1.0 clouds equally distributed front - back and when ξ = 0.0 all clouds are on the front side) 
- `I::Union{Function,Float64}=IsotropicIntensity`: Intensity function/value
- `v::Union{Function,Float64}=vCircularCloud`: Velocity function/value
- `kwargs...`: Extra keyword arguments for `I` and `v` functions (see examples)

# Returns:
A model `ring` struct representing the properties of the cloud. 
"""
function drawCloud(;μ::Float64=500.,β::Float64=1.0,F::Float64=0.5,ϕ₀::Float64=0.0,i::Float64=π/4,rot::Float64=0.0,rₛ::Float64=1.0,θₒ::Float64=0.0,θₒSystem::Float64=0.0,
    I::Union{Function,Float64}=IsotropicIntensity,v::Union{Function,Float64}=vCircularCloud,ξ::Float64=1.0,rng::AbstractRNG=Random.GLOBAL_RNG,kwargs...)
    #assuming e = 0.0, vCircular, IsotropicIntensity for now
    # γ = getGamma(μ=μ,β=β,F=F)
    # r = getR(rₛ,γ)
    g = getG(β)
    r = getR(rₛ,μ,β,F,g,rng)
    xyzSys = rotate3D(r,ϕ₀,i,rot,θₒ) #system coordinates xyz, #flip i to match convention of +z being up, relic
    reflect = (xyzSys[3] > midPlaneXZ(xyzSys[1],i)) && (rand(rng) > ξ) #reflect particles from back of disk across disk midplane to front, leaving ξ fraction of particles in back
    if reflect
        xyzSys = reflect!(xyzSys,i)
    end
    undo_tilt = [sin(i) 0.0 cos(i); 0.0 1.0 0.0; -cos(i) 0.0 sin(i)]
    xyzPlane = undo_tilt*xyzSys
    ϕ = atan(xyzPlane[2],-xyzPlane[1]) #ϕ after rotation, measured from +x in disk plane, - sign relic of how rotation matrix was implemented + desire to have ϕ=0 at +x
    η = response(r;kwargs...)
    return ring(r=r,i=i,v=v,I=I,ϕ=ϕ,ΔA=1.0,Δr=1.0,Δϕ=1.0,scale=nothing,rot=rot,θₒ=θₒ,ϕ₀=ϕ₀,reflect=reflect,rng=rng,η=η;kwargs...)
end