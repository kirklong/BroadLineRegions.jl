#!/usr/bin/env julia
using LinearAlgebra: ⋅, *, inv, mul!, norm, cross
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
    return rMin, rMax
end

function getϕ_rMat(iSystem,θₒSystem)
    """calculate the rotation matrix to rotate the system plane normal into the XY plane
    params:
        iSystem: inclination angle of system (rad) {Float64}
        θₒSystem: opening angle of system {Float64}
    returns:
        rMat: rotation matrix {Array{Float64,2}}
    """
    P = [-cos(iSystem)*cos(θₒSystem/2)-sin(iSystem)*sin(θₒSystem/2),0,cos(θₒSystem/2)*sin(iSystem)-cos(iSystem)*sin(θₒSystem/2)] #average plane normal of disk
    XY = [0.,0.,1.] #want to rotate disk plane into XY plane to calculate ϕ
    c = (P⋅XY)/(norm(P)*norm(XY)) #cosine of angle between P and XY
    ax = cross(P,XY)/norm(cross(P,XY)) #axis of rotation
    s = sqrt(1-c^2) #sine of angle between P and XY
    rMat = [ax[1]^2*(1-c)+c ax[1]*ax[2]*(1-c)-ax[3]*s ax[1]*ax[3]*(1-c)+ax[2]*s;
            ax[2]*ax[1]*(1-c)+ax[3]*s ax[2]^2*(1-c)+c ax[2]*ax[3]*(1-c)-ax[1]*s;
            ax[3]*ax[1]*(1-c)-ax[2]*s ax[3]*ax[2]*(1-c)+ax[1]*s ax[3]^2*(1-c)+c] #rotation matrix
    return rMat
end

getα(r::Float64,ϕ::Float64,i::Float64,rot::Float64,θₒ::Float64) = r*sin(ϕ)*cos(rot)-r*cos(ϕ)*cos(θₒ)*sin(rot) #only valid for non-reflecting points
getβ(r::Float64,ϕ::Float64,i::Float64,rot::Float64,θₒ::Float64) = cos(i)*(r*cos(ϕ)*cos(rot)*cos(θₒ)+r*sin(ϕ)*sin(rot))+sin(i)*sin(θₒ)*r*cos(ϕ) #only valid for non-reflecting points                   






