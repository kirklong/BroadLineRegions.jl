#!/usr/bin/env julia

function vCircular(;r::Union{Float64,Vector{Float64}}, i::Float64, ϕ::Union{Vector{Float64},Float64}, rₛ=1.0, _...) 
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
    return @. sqrt(rₛ/(2*r))*$sin(i)*cos(ϕ)
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