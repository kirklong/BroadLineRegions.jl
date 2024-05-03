#!/usr/bin/env julia
function vCircular(;r::Float64, i::Float64, ϕ::Vector{Float64}, G::Float64=1.0, M::Float64=1.0) 
    """calculate line of sight velocity for circular orbit at radius r from central mass and inclined at angle i (rad) over grid of azimuthal angles ϕ (rad)
    params: 
        r: radius from central mass {Float64}
        i: inclination angle (rad) {Float64}
        ϕ: list of azimuthal angles (rad) {Vector{Float64}}
        G: gravitational constant {Float64}
        M: central mass {Float64}
    returns:
        line of sight velocity {Vector{Float64}}
    """
    return sqrt.(G*M/r).*sin(i).*cos.(ϕ)
end

function vElliptical(;r::Float64, i::Float64, ϕ::Vector{Float64}, G::Float64=1.0, M::Float64=1.0, e::Float64=0.0) 
    """calculate line of sight velocity for elliptical orbit at radius r from central mass and inclined at angle i (rad) over grid of azimuthal angles ϕ (rad)
    params: 
        r: semimajor axis of ellipse {Float64}
        i: inclination angle (rad) {Float64}
        ϕ: list of azimuthal angles (rad) {Vector{Float64}}
        G: gravitational constant {Float64}
        M: central mass {Float64}
        e: eccentricity {Float64}
    returns:
        line of sight velocity {Vector{Float64}}, radius {Vector{Float64}}
    """
    v² = G*M*((1 .+ 2*e.*cos.(ϕ) .+ e^2)./(r*(1 - e^2)))
    r = (r*(1 - e^2))./(1 .+ e.*cos.(ϕ)) # r as a function of ϕ, initial r is taken as the semi-major axis
    return sqrt.(v²).*sin(i).*cos.(ϕ), r
end

function DiskWindIntensity(;r::Float64, i::Float64, ϕ::Vector{Float64}, f1::Float64, f2::Float64, f3::Float64, f4::Float64, α::Float64)
"""calculate the intensity of the disk wind at radius r from the central mass and inclined at angle i (rad) over grid of azimuthal angles ϕ (rad)
    Follows the prescription given in Long+ 2023 and 2024, similar to Chiang and Murray 1996, assumes optically thick line emission limit
    params: 
        r: radius from central mass (in terms of rₛ) {Float64}
        i: inclination angle (rad) {Float64}
        ϕ: list of azimuthal angles (rad) {Vector{Float64}}
        f1: strength of radial velocity gradient dvᵣ/dr {Float64}
        f2: strength of Keplerian shear dvϕ/dr {Float64}
        f3: strength of velocity gradient dvθ/dr {Float64}
        f4: strength of velocity gradient dvθ/dθ (or isoptropic emission) {Float64}
        α: power-law index of the source function S(r) ∝ r^(-α) {Float64}
    returns:
        intensity (arbitrary units) {Vector{Float64}}
"""
    pre = sqrt.(1 ./(2 .* r.^3)); cosϕ = cos.(ϕ); sinϕ = sin.(ϕ); sini = sin(i); cosi = cos(i)
    term12 = (3*sini^2).*(cosϕ) .* (√2*f1 .* cosϕ .+ f2/2 .*sinϕ)
    term3 = ((-f3*3*sini*cosi)) .* cosϕ
    term4 = √2*f4*cosi^2 .* ones(length(ϕ))
    dvl_dl = pre .* (term12 .+ term3 .+ term4)
    I = r.^(-α) .* dvl_dl
    return I
end

function IsotropicIntensity(;ϕ::Union{Vector{Float64},Float64},rescale::Float64=1.0)
    """calculate the intensity of the isotropic emission over grid of azimuthal angles ϕ (rad) or at one azimuthal angle
    params: 
        ϕ: list of azimuthal angles or single angle (rad) {Vector{Float64}}
        rescale: rescale factor for the intensity {Float64}
    returns:
        intensity (arbitrary units) {typeof(ϕ)}
    """
    return typeof(ϕ) == Float64 ? rescale : rescale.*ones(length(ϕ))
end

@kwdef struct ring{F<:Float64,V<:Vector{F}} 
    r::Union{V,F}
    i::F
    e::F = 0.0
    v::Union{V,F,Function}
    I::Union{V,F,Function}
    ϕ::Union{V,F}
    function ring(;r,i,ϕ,v=vCircular,I=DiskWindIntensity,e=0.0)
        @assert typeof(r) == Float64 "r must be Float64, got $(typeof(r))"
        @assert typeof(i) == Float64 "i must be Float64, got $(typeof(i))"
        @assert typeof(e) == Float64 "e must be Float64, got $(typeof(e))"

    #     @assert
    end
end
