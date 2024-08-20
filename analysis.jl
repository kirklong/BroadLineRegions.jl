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

function raytrace(α::Float64, β::Float64, i::Float64, rot::Float64, θₒPoint::Float64)
    """calculate where ray traced back from camera coordinates r_c, ϕ_c intersects the system (assumes circular geometry)
    params:
        α: image x coordinate (in terms of rₛ) {Float64}
        β: image y coordinate (rad) {Float64}
        i: inclination angle of system (rad) {Float64}
        rot: rotation of current point about z axis (rad) {Float64}
        θₒPoint: opening angle of current point {Float64}
    returns:
        r: radius of system ring plane at intersection {Float64}
        ϕ: azimuthal angle of system ring plane at intersection {Float64}
    """

    ###NOTE: THIS IS STUPID
    #for real raytracing, fix y and z (α and β) and step through x, checking if that matches a system point or not until you go "through" the system. 
    #use rotation matrix at each step to invert to original coordinates? or just check from converted coordinates 
    xRing = (β*cos(rot) - α*cos(i)*sin(rot))/(cos(i)*cos(θₒPoint)+cos(rot)*sin(i)*sin(θₒPoint)) #system x
    yRing = (α*(cos(i)*cos(θₒPoint)+sec(rot)*sin(i)*sin(θₒPoint))+β*cos(θₒPoint)*tan(rot))/(cos(i)*cos(θₒPoint)*sec(rot)+sin(i)*sin(θₒPoint)) #system y
    r = √(xRing^2 + yRing^2)
    ϕₒ = atan(yRing,xRing) #original ϕₒ (no rotation)
    xyzSys = rotate3D(r,ϕₒ,i,rot,θₒPoint) #system coordinates xyz
    ϕ = atan(xyzSys[2],xyzSys[1]) #ϕ after rotation, measured from +x in disk plane
    #rMat = getϕ_rMat(i,θₒSystem) #rotation matrix
    #xyzXY = rMat*xyzSys #rotate system plane into XY plane
    #ϕ = atan(xyzXY[2],xyzXY[1]) #ϕ after rotation, measured from +x in disk plane
    return r, ϕ, ϕₒ
end


#NOTE FOR RAYTRACING: i don't think I need to do this -- just take atan(y,x) of the system coordinates because camera is at +x in 3D space
function raytrace(α::Float64, β::Float64, i::Float64, rot::Float64, θₒPoint::Float64, r3D::Matrix{Float64}, xyz::Vector{Float64}, matBuff::Matrix{Float64}, colBuff::Vector{Float64})
    """performant version of raytrace function -- calculate where ray traced back from camera coordinates r_c, ϕ_c intersects the system (assumes circular geometry)
    params:
        α: image x coordinate (in terms of rₛ) {Float64}
        β: image y coordinate {Float64}
        i: inclination angle of system (rad) {Float64}
        rot: rotation of current point about z axis (rad) {Float64}
        θₒPoint: opening angle of current point {Float64}
        ϕ_rMat: matrix that rotates system plane into XY plane {Matrix{Float64}}
        r3D: 3D matrix that rotates initial system coordinates into XYZ space {Matrix{Float64}}
        xyz: preallocated xyz vector (but not precalculated) {Vector{Float64}}
        matBuff: preallocated buffer matrix for storing result of 3x3 matrix multiplication {Matrix{Float64}}
        colBuff: preallocated buffer vector for storing final matrix multiplication result {Vector{Float64}}
    returns:
        r: radius of system ring plane at intersection {Float64}
        ϕ: azimuthal angle of system ring plane at intersection {Float64}
        ϕₒ: original azimuthal angle in ring plane {Float64}
    """
    cosr = cos(rot); sinr = sin(rot); cosi = cos(i); sini = sin(i); cosθₒ = cos(θₒPoint); sinθₒ = sin(θₒPoint)
    xRing = (β*cosr - α*cosi*sinr)/(cosi*cosθₒ+cosr*sini*sinθₒ) #system x
    yRing = (α*(cosi*cosθₒ+sini/cosr*sinθₒ)+β*cosθₒ*sinr/cosr)/(cosi*cosθₒ/cosr+sini*sinθₒ)
    r = √(xRing^2 + yRing^2)
    ϕₒ = atan(yRing,xRing) #original ϕₒ (no rotation)
    xyz[1] = xRing; xyz[2] = yRing; xyz[3] = 0.0
    # xyz[1] = r*cos(ϕₒ); xyz[2] = r*sin(ϕₒ); xyz[3] = 0.0
    # mul!(matBuff,ϕ_rMat,r3D)
    # mul!(colBuff,matBuff,xyz) #rotate system plane into XY plane
    mul!(colBuff,r3D,xyz)
    undo_tilt = [sini 0.0 cosi; 0.0 1.0 0.0; -cosi 0.0 sini]
    mul!(xyz,undo_tilt,colBuff)
    ϕ = atan(xyz[2],xyz[1]) #ϕ after rotation and being "puffed up", measured from +x in disk plane -- this is fine even for puffed up clouds but note ϕ is measured wrt to disk midplane then
    #really this whole thing is stupid and should be removed this is not raytracing
    return r, ϕ, ϕₒ
end

getα(r::Float64,ϕ::Float64,i::Float64,rot::Float64,θₒ::Float64) = r*sin(ϕ)*cos(rot)-r*cos(ϕ)*cos(θₒ)*sin(rot) #only valid for non-reflecting points
getβ(r::Float64,ϕ::Float64,i::Float64,rot::Float64,θₒ::Float64) = cos(i)*(r*cos(ϕ)*cos(rot)*cos(θₒ)+r*sin(ϕ)*sin(rot))+sin(i)*sin(θₒ)*r*cos(ϕ) #only valid for non-reflecting points

function photograph(r::Float64, ϕₒ::Float64, i::Float64, rot::Float64, θₒ::Float64, reflect::Bool=false)
    """calculate the image coordinates from system coordinates r, ϕ + inclination angle i
    params:
        r: radius from central mass (in terms of rₛ) {Float64}
        ϕₒ: unrotated azimuthal angle in ring plane (rad) {Float64}
        i: inclination angle of ring plane (rad) {Float64}
        rot: rotation of system plane about z axis (rad) {Float64}
        θₒ: ring opening angle {Float64}
        reflect: whether the point is reflected across the midplane of the disk {Bool}
    returns:
        α: image x coordinate (in terms of rₛ) {Float64}
        β: image y coordinate {Float64}
    """
    xyzSys = rotate3D(r,ϕₒ,i,rot,θₒ,reflect)
    #α = getα(r,ϕₒ,i,rot,θₒ)
    #β = getβ(r,ϕₒ,i,rot,θₒ)
    α = xyzSys[2] #camera is at +x, so α is y
    β = xyzSys[3] #and β is z
    return α, β
end