#!/usr/bin/env julia

function get_rMinMaxDiskWind(r̄::Float64,rFac::Float64,α::Float64) 
    """calculate the minimum and maximum radius of model given the (intensity weighted) mean radius r̄, the radius factor rFac, and the power-law index α
    This is valid for the DiskWind model detailed in Long+ 2023
    params:
        r̄: mean radius of model (in terms of rₛ) {Float64}
        rFac: radius factor {Float64} 
        α: power-law index of the source function S(r) ∝ r^(-α) {Float64} (cannot be 1/2 or 3/2 as this divides by zero)
    returns:
        rMin: minimum radius of model (in terms of rₛ) {Float64}
        rMax: maximum radius of model (in terms of rₛ) {Float64}
    """
    rMin = r̄*(3-2*α)/(1-2*α)*(rFac^(1/2-α)-1)/(rFac^(3/2-α)-1)
    rMax = rMin*rFac
    return rMin,rMax
end

function raytrace(r_c::Float64, ϕ_c::Float64, i::Float64)
    """calculate where ray traced back from camera coordinates r_c, ϕ_c intersects the system
    params:
        r_c: radius from image center (in terms of rₛ) {Float64}
        ϕ_c: image azimuthal angle (rad) {Float64}
        i: inclination angle (rad) {Float64}
    returns:
        r: radius of system ring plane at intersection {Float64}
        ϕ: azimuthal angle of system ring plane at intersection {Float64}
    """
    α = r_c*cos(ϕ_c); β = r_c*sin(ϕ_c) #camera coordinates
    r = √(α^2 + (β/cos(i))^2) 
    ϕ = atan(β/cos(i),α)
    return r,ϕ
end

function photograph(r::Float64, ϕ::Float64, i::Float64)
    """calculate the image coordinates from system coordinates r, ϕ + inclination angle i
    params:
        r: radius from central mass (in terms of rₛ) {Float64}
        ϕ: azimuthal angle in ring plane (rad) {Float64}
        i: inclination angle of ring plane (rad) {Float64}
    returns:
        r_c: radius from image center (in terms of rₛ) {Float64}
        ϕ_c: image azimuthal angle (rad) {Float64}
    """
    x = r*cos(ϕ); y = r*sin(ϕ) #system x and y
    r_c = √(x^2 + (y*cos(i))^2)
    ϕ_c = atan(y*cos(i),x)
    return r_c,ϕ_c
end