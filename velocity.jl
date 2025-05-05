#!/usr/bin/env julia
using Random

vCirc(r::Float64, rₛ::Float64=1.0) = √(rₛ/(2*r))

function vCircularDisk(;r::Union{Float64,Vector{Float64}}, i::Float64, ϕ::Union{Vector{Float64},Float64}, θₒ::Union{Vector{Float64},Float64}, rₛ=1.0, _...) 
    """calculate line of sight velocity for circular orbit at radius r from central mass and inclined at angle i (rad) over grid of azimuthal angles ϕ (rad)
    params: 
        r: radius from central mass {Float64} (rₛ)
        i: inclination angle (rad) {Float64}
        ϕ: list of azimuthal angles (rad) {Vector{Float64}}
        rₛ: Schwarzschild radius {Float64} (optional, to convert to physical units, defaults to 1)
        _: extra kwargs, ignored
    returns:
        line of sight velocity {Vector{Float64}}
    """
    return @. vCirc(r,$rₛ)*$sin(i)*sin(ϕ) #circular velocity where sides of disk are at ±π/2
end

function vCircularRadialDisk(;r::Union{Float64,Vector{Float64}}, i::Float64, ϕ::Union{Vector{Float64},Float64}, vᵣFrac::Union{Vector{Float64},Float64}=0.0, inflow::Union{Vector{Bool},Bool}=true, rₛ=1.0, _...) 
    """
    calculate line of sight velocity for circular orbit at radius r from central mass and inclined at angle i (rad) over grid of azimuthal angles ϕ (rad) with radial inflow/outflow
    params: 
        r: radius from central mass {Float64} (rₛ)
        i: inclination angle (rad) {Float64}
        ϕ: list of azimuthal angles (rad) {Vector{Float64}}
        vᵣFrac: fraction of radial velocity {Union{Vector{Float64},Float64}} (defaults to 0.0)
        inflow: whether the radial velocity is inflowing {Union{Vector{Bool},Bool}} (defaults to true)
        rₛ: Schwarzschild radius {Float64} (optional, to convert to physical units, defaults to 1)
        _: extra kwargs, ignored
    returns:
        line of sight velocity Union{Vector{Float64},Float64}
    """
    vsini = @. vCirc(r,$rₛ)*$sin(i)
    inflow = @. inflow ? 1.0 : -1.0
    vᵣ = @. vsini*cos(ϕ)*vᵣFrac*inflow
    vϕ = @. vsini*sin(ϕ)*(1-vᵣFrac)
    return vᵣ + vϕ
end

function vCircularCloud(;r::Float64, ϕₒ::Float64, i::Float64, rot::Float64, θₒ::Float64, rₛ::Float64=1.0, reflect::Bool=false, _...)
    """calculate line of sight velocity for cloud in 3D space
    params:
        r: radius from central mass (in terms of rₛ) {Float64}
        ϕₒ: starting azimuthal angle in ring plane (rad) {Float64}
        i: inclination angle of ring plane (rad) {Float64}
        rot: rotation of system plane about z axis (rad) {Float64}
        θₒ: opening angle of point {Float64}
        rₛ: Schwarzschild radius {Float64} (optional, to convert to physical units, defaults to 1)
        reflect: whether the point is reflected across the midplane of the disk {Bool}
        _: extra kwargs, ignored
    returns:
        line of sight velocity {Float64}
    """
    vₒ = vCirc(r,rₛ)
    vXYZ = [-vₒ*sin(ϕₒ),vₒ*cos(ϕₒ),0.0] #match velocity sign conventions such that left side is coming towards observer
    r3D = get_r3D(i,rot,θₒ)
    vXYZ = r3D*vXYZ
    if reflect
        vXYZ = BLR.reflect!(vXYZ,i)
    end
    return vXYZ[1] #line of sight velocity is x component after rotation (camera is at +x)
end

function vCloudTurbulentEllipticalFlow(;σρᵣ::Float64,σρc::Float64, σΘᵣ::Float64, σΘc::Float64, θₑ::Float64, fEllipse::Float64, fFlow::Float64, σₜ::Float64, 
    r::Float64, i::Float64, rot::Float64, θₒ::Float64, rₛ::Float64=1.0, ϕₒ::Float64=0.0, reflect::Bool=false, rng::AbstractRNG=Random.GLOBAL_RNG, _...) 
    """calculate line of sight velocity for cloud in 3D space with potential for elliptical orbital velocities, in/outflow, and turbulence as in Pancoast+14
    params:
        σρᵣ: radial standard deviation around radial orbits {Float64}
        σρc: radial standard deviation around circular orbits {Float64}
        σΘᵣ: angular standard deviation around radial orbits {Float64}
        σΘc: angular standard deviation around circular orbits {Float64}
        θₑ: angle in vᵩ-vᵣ plane {Float64}
        fEllipse: fraction of elliptical orbits {Float64}
        fFlow: if < 0.5, inflow, otherwise, outflow {Float64}
        σₜ: standard deviation of turbulent velocity {Float64}
        r: radius from central mass (in terms of rₛ) {Float64}
        i: inclination angle of ring plane (rad) {Float64}
        rot: rotation of system plane about z axis (rad) {Float64}
        θₒ: opening angle of point {Float64}
        rₛ: Schwarzschild radius {Float64} (optional, to convert to physical units, defaults to 1)
        ϕₒ: starting azimuthal angle in ring plane (rad) {Float64}
        reflect: whether the point is reflected across the midplane of the disk {Bool}
        _: extra kwargs, ignored
    returns:
        line of sight velocity {Float64}
    """
    vc = vCirc(r,rₛ)
    vₜ = rand(rng,Normal(0.0,σₜ))*vc
    ρ = 0.0; Θ = 0.0
    if rand(rng) < fEllipse #elliptical orbit, distribution deviating from circular by σρc, σΘc, Pancoast 14 2.5.1
        ρ = rand(rng,Normal(vc,σρc))
        Θ = π/2 + rand(rng,Normal(0.0,σΘc))
    else #in/outflowing orbit, distribution deviating from circular by σρᵣ, σΘᵣ, Pancoast14 2.5.2
        ρ = rand(rng,Normal(vc,σρᵣ))
        Θ = fFlow < 0.5 ? rand(rng,Normal(0.0,σΘᵣ)) + (π - θₑ) : rand(rng,Normal(0.0,σΘᵣ)) + θₑ
    end
    vx = √2*ρ*cos(Θ); vy = ρ*sin(Θ) #without any rotation, radial direction is along x and ϕ is along y at ϕ = 0 when left is rotating towards observer
    vXYZ = [vx*cos(ϕₒ)-vy*sin(ϕₒ),vx*sin(ϕₒ)+vy*cos(ϕₒ),0.0] #rotate around z by ϕₒ, match velocity sign conventions (left = towards observer)
    r3D = get_r3D(i,rot,θₒ) #transform initial coordinates to system coordinates
    vXYZ = r3D*vXYZ #rotate into system coordinates
    if reflect
        vXYZ = BLR.reflect!(vXYZ,i)
    end
    return -vXYZ[1]+vₜ #line of sight velocity is x component after rotation (camera is at +x), turbulence only along line of sight (see Pancoast14 2.5.3), negative sign to match disk convention left towards observer
end

function vElliptical(;a::Float64, i::Float64, ϕ::Union{Vector{Float64},Float64}, G::Float64=1.0, M::Float64=1.0, e::Float64=0.0, _...) 
    """calculate line of sight velocity for elliptical orbit with semimajor axis a from central mass and inclined at angle i (rad) over grid of azimuthal angles ϕ (rad)
    params: 
        a: semimajor axis of ellipse {Float64}
        i: inclination angle (rad) {Float64}
        ϕ: list of azimuthal angles (rad) {Vector{Float64}}
        G: gravitational constant {Float64}
        M: central mass {Float64}
        e: eccentricity {Float64}
        _: extra kwargs, ignored
    returns:
        line of sight velocity {Vector{Float64}}, radius {Vector{Float64}}
    """
    v² = @. G*M*((1 .+ 2*e.*cos.(ϕ) .+ e^2)./(a*(1 - e^2)))
    return @. sqrt.(v²).*sin(i).*cos.(ϕ)
end

function rElliptical(;a::Float64, ϕ::Union{Vector{Float64},Float64}, e::Float64=0.0, _...) 
    """calculate distances for elliptical orbit with semimajor axis a from central mass over grid of azimuthal angles ϕ (rad)
    params: 
        a: semimajor axis of ellipse {Float64}
        ϕ: list of azimuthal angles (rad) {Vector{Float64}}
        e: eccentricity {Float64}
        _: extra kwargs, ignored
    returns:
        radius {Vector{Float64}}
    """
    return @. (a*(1 - e^2))./(1 .+ e.*cos.(ϕ))
end