#!/usr/bin/env julia
using LinearAlgebra

"""
    raytrace(α::Float64, β::Float64, i::Float64, rot::Float64, θₒPoint::Float64)

Calculate where ray traced back from camera coordinates `α` and `β` intersects the system (assumes circular geometry).

# Arguments
- `α::Float64`: image x coordinate (in terms of rₛ)
- `β::Float64`: image y coordinate (in terms of rₛ)
- `i::Float64`: inclination angle of system (rad)
- `rot::Float64`: how the point was rotated about z axis (rad)
- `θₒPoint::Float64`: opening angle of current point

# Returns
- `r::Float64`: distance from central mass (in terms of rₛ)
- `ϕ::Float64`: azimuthal angle of system ring plane at intersection
- `ϕ₀::Float64`: original azimuthal angle in ring plane (no rotation)

# Note
This function is *coordinate* raytracing only. To raytrace models and combine intensities, see `raytrace!`. 
"""
function raytrace(α::Float64, β::Float64, i::Float64, rot::Float64, θₒPoint::Float64)
    xRing = (β*cos(rot) - α*cos(i)*sin(rot))/(cos(i)*cos(θₒPoint)+cos(rot)*sin(i)*sin(θₒPoint)) #system x
    yRing = (α*(cos(i)*cos(θₒPoint)+sec(rot)*sin(i)*sin(θₒPoint))+β*cos(θₒPoint)*tan(rot))/(cos(i)*cos(θₒPoint)*sec(rot)+sin(i)*sin(θₒPoint)) #system y
    r = √(xRing^2 + yRing^2)
    ϕ₀ = atan(yRing,xRing) #original ϕ₀ (no rotation)
    xyzSys = rotate3D_scalar(r,ϕ₀,i,rot,θₒPoint) #system coordinates xyz
    ϕ = atan(xyzSys[2],xyzSys[1]) #ϕ after rotation, measured from +x in disk plane
    return r, ϕ, ϕ₀
end

"""
    raytrace(α::Float64, β::Float64, i::Float64, rot::Float64, θₒPoint::Float64, M::Matrix{Float64})

Allocation-free version of `raytrace` — calculate where ray traced back from camera coordinates `α`, `β`
intersects the system (assumes circular geometry), using a single precomputed rotation matrix.

`M` must be `undo_tilt * get_r3D(i, rot, θₒPoint)` where
`undo_tilt = [sin(i) 0.0 -cos(i); 0.0 1.0 0.0; cos(i) 0.0 sin(i)]` — both factors are constant for fixed
`(i, rot, θₒ)` so `M` should be computed once outside any pixel loop.

# Returns
- `r::Float64`: distance from central mass (in terms of rₛ)
- `ϕ::Float64`: azimuthal angle of system ring plane at intersection
- `ϕ₀::Float64`: original azimuthal angle in ring plane

# Note
Supersedes the buffer-based 9-argument method (which is retained for compatibility).
"""
function raytrace(α::Float64, β::Float64, i::Float64, rot::Float64, θₒPoint::Float64, M::Matrix{Float64})
    cosr = cos(rot); sinr = sin(rot); cosi = cos(i); sini = sin(i); cosθₒ = cos(θₒPoint); sinθₒ = sin(θₒPoint) #flip i to match convention that +x is closer
    xRing = -(β*cosr - α*cosi*sinr)/(cosi*cosθₒ+cosr*sini*sinθₒ) #system x, negative because want bottom towards camera
    yRing = (α*(cosi*cosθₒ+sini/cosr*sinθₒ)+β*cosθₒ*sinr/cosr)/(cosi*cosθₒ/cosr+sini*sinθₒ)
    r = √(xRing^2 + yRing^2)
    ϕ₀ = atan(yRing,xRing) #original ϕ₀ (no rotation)
    x = M[1,1]*xRing + M[1,2]*yRing #z component of ring coordinates is 0 so third column never contributes
    y = M[2,1]*xRing + M[2,2]*yRing
    ϕ = atan(y,x) #ϕ after rotation and being "puffed up", measured from +x in disk plane
    return r, ϕ, ϕ₀
end

"""
    raytrace(α::Float64, β::Float64, i::Float64, rot::Float64, θₒPoint::Float64, M::Matrix{Float64}, r3D::Matrix{Float64})

Same as the 6-argument `raytrace(α, β, i, rot, θₒPoint, M)` method but additionally returns the 3D
*system* coordinates of the intersection point (camera at +x), computed from the supplied
`r3D = get_r3D(i, rot, θₒPoint)` matrix. Used at model construction to fill the per-point
coordinate cache (`ring.x/y/z`) without a second rotation pass.

# Returns
- `r, ϕ, ϕ₀, x, y, z` (all `Float64`)
"""
function raytrace(α::Float64, β::Float64, i::Float64, rot::Float64, θₒPoint::Float64, M::Matrix{Float64}, r3D::Matrix{Float64})
    cosr = cos(rot); sinr = sin(rot); cosi = cos(i); sini = sin(i); cosθₒ = cos(θₒPoint); sinθₒ = sin(θₒPoint) #flip i to match convention that +x is closer
    xRing = -(β*cosr - α*cosi*sinr)/(cosi*cosθₒ+cosr*sini*sinθₒ) #system x, negative because want bottom towards camera
    yRing = (α*(cosi*cosθₒ+sini/cosr*sinθₒ)+β*cosθₒ*sinr/cosr)/(cosi*cosθₒ/cosr+sini*sinθₒ)
    r = √(xRing^2 + yRing^2)
    ϕ₀ = atan(yRing,xRing) #original ϕ₀ (no rotation)
    x = r3D[1,1]*xRing + r3D[1,2]*yRing #system coordinates (z component of ring coordinates is 0 so third columns never contribute)
    y = r3D[2,1]*xRing + r3D[2,2]*yRing
    z = r3D[3,1]*xRing + r3D[3,2]*yRing
    xd = M[1,1]*xRing + M[1,2]*yRing #disk-plane ("puffed up" + tilt undone) coordinates for ϕ
    yd = M[2,1]*xRing + M[2,2]*yRing
    ϕ = atan(yd,xd) #ϕ after rotation and being "puffed up", measured from +x in disk plane
    return r, ϕ, ϕ₀, x, y, z
end

"""
    raytrace(α::Float64, β::Float64, i::Float64, rot::Float64, θₒPoint::Float64, 
            r3D::Matrix{Float64}, xyz::Vector{Float64}, matBuff::Matrix{Float64}, 
            colBuff::Vector{Float64})

Performant version of `raytrace` function -- calculate where ray traced back from camera coordinates `α`, `β` intersects the system (assumes circular geometry).

# Arguments
- `α::Float64`: image x coordinate (in terms of rₛ)
- `β::Float64`: image y coordinate (in terms of rₛ)
- `i::Float64`: inclination angle of system (rad)
- `rot::Float64`: rotation of current point about z axis (rad)
- `θₒPoint::Float64`: opening angle of current point
- `r3D::Matrix{Float64}`: matrix that rotates system plane into XY plane
- `xyz::Vector{Float64}`: preallocated xyz vector (but not precalculated)
- `matBuff::Matrix{Float64}`: preallocated buffer matrix for storing result of 3x3 matrix multiplication
- `colBuff::Vector{Float64}`: preallocated buffer vector for storing final matrix multiplication result

# Returns
- `r::Float64`: distance from central mass (in terms of rₛ)
- `ϕ::Float64`: azimuthal angle of system ring plane at intersection
- `ϕ₀::Float64`: original azimuthal angle in ring plane

# Note
This function is *coordinate* raytracing only. To raytrace models and combine intensities, see `raytrace!`.

!!! warning "Deprecated"
    Prefer the 6-argument method `raytrace(α, β, i, rot, θₒPoint, M)` with a single precomputed
    `M = undo_tilt * get_r3D(i, rot, θₒ)` — it is allocation-free and faster. This buffer-based method
    is retained for backwards compatibility only.
"""
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
        ϕ₀: original azimuthal angle in ring plane {Float64}
    """
    cosr = cos(rot); sinr = sin(rot); cosi = cos(i); sini = sin(i); cosθₒ = cos(θₒPoint); sinθₒ = sin(θₒPoint) #flip i to match convention that +x is closer
    xRing = -(β*cosr - α*cosi*sinr)/(cosi*cosθₒ+cosr*sini*sinθₒ) #system x, negative because want bottom towards camera
    yRing = (α*(cosi*cosθₒ+sini/cosr*sinθₒ)+β*cosθₒ*sinr/cosr)/(cosi*cosθₒ/cosr+sini*sinθₒ)
    r = √(xRing^2 + yRing^2)
    ϕ₀ = atan(yRing,xRing) #original ϕ₀ (no rotation)
    xyz[1] = xRing; xyz[2] = yRing; xyz[3] = 0.0
    mul!(colBuff,r3D,xyz)
    undo_tilt = [sini 0.0 -cosi; 0.0 1.0 0.0; cosi 0.0 sini] 
    mul!(xyz,undo_tilt,colBuff)
    ϕ = atan(xyz[2],xyz[1]) #ϕ after rotation and being "puffed up", measured from +x in disk plane 
    return r, ϕ, ϕ₀
end

"""
    photograph(r::Float64, ϕ₀::Float64, i::Float64, rot::Float64, θₒ::Float64, reflect::Bool=false)

Calculate the image coordinates from system coordinates r, ϕ + inclination angle i.

# Arguments
- `r::Float64`: radius from central mass (in terms of rₛ)
- `ϕ₀::Float64`: unrotated azimuthal angle in ring plane (rad)
- `i::Float64`: inclination angle of ring plane (rad)
- `rot::Float64`: rotation of system plane about z axis (rad)
- `θₒ::Float64`: ring opening angle
- `reflect::Bool=false`: whether the point is reflected across the midplane of the disk

# Returns
- `α::Float64`: image x coordinate (in terms of rₛ)
- `β::Float64`: image y coordinate (in terms of rₛ)

# Note
This function is *coordinate* photography only. To visualize models, see `Image`.`
"""
function photograph(r::Float64, ϕ₀::Float64, i::Float64, rot::Float64, θₒ::Float64, reflect::Bool=false)
    xyzSys = rotate3D_scalar(r,ϕ₀,i,rot,θₒ,reflect)
    α = xyzSys[2] #camera is at +x, so α is y
    β = xyzSys[3] #and β is z
    return α, β
end

"""
    zeroDiskObscuredClouds!(m::model; diskCloudIntensityRatio::Float64=1.0, rotate3D::Function=rotate3D_scalar)

Zero out the intensities of clouds that are obscured by the disk.

Performs simple raytracing for an optically thick obscuring disk. The function
modifies the input model by setting the intensity of obscured cloud points to zero
and adjusting the disk intensity according to the specified ratio.

# Arguments
- `m::model`: Model to zero out disk obscured clouds. Should be a combined model consisting of a disk component and a cloud component. 
- `diskCloudIntensityRatio::Float64=1.0`: Ratio of disk to cloud intensity, used to scale 
  the disk intensities after zeroing out clouds
- `rotate3D::Function=rotate3D_scalar`: Function to rotate coordinates in 3D space (must return an indexable `(x,y,z)`)

# Returns
- `m::model`: Model with disk obscured clouds zeroed out

# See also 
- `removeDiskObscuredClouds!`: Function to remove disk obscured clouds instead of zeroing them out
"""
function zeroDiskObscuredClouds!(m::model;diskCloudIntensityRatio::Float64=1.,rotate3D::Function=rotate3D_scalar)
    isCombined = length(m.subModelStartInds) > 1 #check if model is combined
    if !isCombined
        @warn "did not detect combined model so nothing to zero -- returning unaltered input model"
        return m
    end
    # The logic below classifies rings by type (disk = vector ϕ/r, cloud = scalar),
    # not by submodel, so it handles any number of submodels (disk + N cloud groups)
    # in one pass. (The old `length(startInds) > 2` recursion sliced rings without
    # slicing the camera/subModelStartInds and called itself positionally — broken.)
    diskFlagRing = [!(typeof(r.ϕ) == Float64 && typeof(r.r) == Float64) for r in m.rings]
    iDisk = m.rings[diskFlagRing][1].i
    for ring in m.rings[diskFlagRing]
        if ring.i != iDisk
            @warn "detected at least two different inclinations in disk rings"
            @warn "this method assumes all disk cells have the same inclination -- returning unaltered input model"
            return m
        end
    end
    for ring in m.rings[.!diskFlagRing]
        xyzCloud = rotate3D === rotate3D_scalar ? getXYZ(ring) : rotate3D(ring.r,ring.ϕ₀,ring.i,ring.rot,ring.θₒ,ring.reflect) #cached system coordinates xyz unless a custom rotate3D was passed
        zDisk = midPlaneXZ(xyzCloud[1],iDisk) #z value of disk at x value of cloud
        if xyzCloud[3] < zDisk #cloud below disk -- invisible to camera
            ring.I = 0.0
        end
    end
    totalCloudI = sum([ring.I for ring in m.rings[.!diskFlagRing]])
    totalDiskI = sum([sum(ring.I) for ring in m.rings[diskFlagRing]])
    ratio = totalDiskI/totalCloudI
    for ring in m.rings[diskFlagRing]
        ring.I *= diskCloudIntensityRatio/ratio
    end
    reset!(m)
    return m
end

"""
    removeDiskObscuredClouds!(m::model, rotate3D::Function=rotate3D_scalar)

Remove clouds that are obscured by the disk.

Performs simple raytracing for an optically thick obscuring disk. The function
modifies the input model by removing cloud points that are obscured by the disk.
Note that this is a mutating operation and the input model will be modified in place.

# Arguments
- `m::model`: Model to remove disk obscured clouds. Should be a combined model consisting 
  of a disk component and a cloud component.
- `rotate3D::Function=rotate3D_scalar`: Function to rotate coordinates in 3D space (must return an indexable `(x,y,z)`)

# Returns
- `m::model`: Model with disk obscured clouds removed

# See also 
- `zeroDiskObscuredClouds!`: Function to zero out disk obscured clouds instead of removing them
"""
function removeDiskObscuredClouds!(m::model,rotate3D::Function=rotate3D_scalar)
    isCombined = length(m.subModelStartInds) > 1 #check if model is combined
    if !isCombined
        @warn "did not detect combined model so nothing to remove -- returning unaltered input model"
        return m
    end
    # Classifies rings by type (disk vs cloud), not submodel, so any number of
    # submodels works in one pass. (The old `length(startInds) > 2` recursion sliced
    # rings without the camera/subModelStartInds and was logically redundant.)
    removeFlag = falses(length(m.rings))
    diskRingInds = [i for i in 1:length(m.rings) if !(typeof(m.rings[i].ϕ) == Float64 && typeof(m.rings[i].r) == Float64)]
    cloudRingInds = [i for i in 1:length(m.rings) if typeof(m.rings[i].ϕ) == Float64 && typeof(m.rings[i].r) == Float64]
    αβStartInds = zeros(Int64,length(m.rings))
    currentαβInd = 1
    for (i,ring) in enumerate(m.rings)
        if typeof(ring.ϕ) == Float64 && typeof(ring.r) == Float64
            αβStartInds[i] = currentαβInd
            currentαβInd += 1
        else
            αβStartInds[i] = currentαβInd
            currentαβInd += length(ring.ϕ)
        end
    end
    diskFlag = Array{Bool}(undef,length(m.camera.α))
    ind = 1
    for ring in m.rings
        if typeof(ring.ϕ) == Float64 && typeof(ring.r) == Float64
            diskFlag[ind] = false
            ind += 1
        else
            stride = length(ring.ϕ)
            diskFlag[ind:ind+stride-1] .= true
            ind += stride
        end
    end
    diskα = m.camera.α[diskFlag]; diskβ = m.camera.β[diskFlag]
    rDisk = sqrt.(diskα.^2 .+ diskβ.^2)
    rMinDisk = minimum(rDisk); rMaxDisk = maximum(rDisk)
    rUnique = unique(r -> round(r,sigdigits=12), rDisk)
    ϕList = [m.rings[i].ϕ for i in 1:length(m.rings)]
    for i in cloudRingInds #check if r inside rMin/max, then find closest disk cell, then compare xs and flag for removal if xCloud behind xDisk
        xyzSys = rotate3D === rotate3D_scalar ? getXYZ(m.rings[i]) : rotate3D(m.rings[i].r,m.rings[i].ϕ₀,m.rings[i].i,m.rings[i].rot,m.rings[i].θₒ,m.rings[i].reflect) #cached system coordinates xyz unless a custom rotate3D was passed
        xCloud = xyzSys[1]
        α,β = m.camera.α[αβStartInds[i]],m.camera.β[αβStartInds[i]]
        rCam = sqrt(α^2 + β^2)
        ϕCam = atan(β,α)
        if (rCam > rMinDisk) && (rCam < rMaxDisk)
            rDiskPos = argmin(abs.(rCam .- rUnique))
            rDiskInd = diskRingInds[rDiskPos]
            ϕDiskInd = argmin(abs.(ϕCam .- ϕList[rDiskInd]))
            xDisk = rotate3D === rotate3D_scalar ? getXYZ(m.rings[rDiskInd])[1][ϕDiskInd] : rotate3D(m.rings[rDiskInd].r[ϕDiskInd],m.rings[rDiskInd].ϕ₀[ϕDiskInd],m.rings[rDiskInd].i,m.rings[rDiskInd].rot,m.rings[rDiskInd].θₒ,m.rings[rDiskInd].reflect)[1] #system coordinates xyz (cached x vector indexed at this point unless a custom rotate3D was passed)
            if xCloud < xDisk #cloud behind disk
                removeFlag[i] = true
            end
        end
    end
    newα = deepcopy(m.camera.α); newβ = deepcopy(m.camera.β)
    for (αβi, flag) in zip(αβStartInds,removeFlag)
        if flag
            newα[αβi] = NaN
            newβ[αβi] = NaN
        end
    end
    newα = filter!(x -> !isnan(x),newα)
    newβ = filter!(x -> !isnan(x),newβ)
    m.camera = camera(newα,newβ,false)
    m.rings = m.rings[.!removeFlag]
    reset!(m)
    return m
end

"""
    raytrace!(m::model; IRatios::Union{Float64,Array{Float64,}}=1.0, 
            τCutOff::Float64=1.0, raytraceFreeClouds::Bool=false)

Perform raytracing for a model, combining overlapping components along line of sight.

!!! warning "Slow"
    This function not very performant and can take a long time to combine large models.
    Consider using [`removeDiskObscuredClouds!`](@ref BLR.removeDiskObscuredClouds!) for
    simple disk obscuration removal if you do not need full raytracing.

This function should be called after combining all relevant models (i.e. `mCombined = m1 + m2 + m3...`).
It performs raytracing in discrete steps (no absorption, only adding intensity in chunks along 
the line of sight until maximum optical depth `τ` is reached) and generates a new model object 
with extraneous points removed. Note that this function will mutate the input model objects.

# Arguments
- `m::model`: Model to raytrace
- `IRatios::Union{Float64,Array{Float64,}}=1.0`: Global emissivity weights for each submodel
  - If `Float64`, applies to all submodels equally
  - If array, applies to each submodel individually (must match number of submodels)
  - Used when combining models with different intensity functions if they aren't properly normalized
- `τCutOff::Float64=1.0`: Maximum optical depth to raytrace to (stops when `τ > τCutOff`)
- `raytraceFreeClouds::Bool=false`: Whether to raytrace free clouds (cloud-cloud raytracing)
  - If `false`, clouds are only raytraced if they overlap with a continuous model
  - If `true`, clouds will be checked for overlap with other clouds and raytraced accordingly

# Returns
- `m::model`: Model with raytraced points
"""

struct _RaytracePoint
    submodel::Int
    ring::Int
    col::Int
    α::Float64
    β::Float64
    rCam::Float64
    ϕCam::Float64
    x::Float64
    r::Float64
    ϕ::Float64
    ϕ₀::Float64
    i::Float64
    rot::Float64
    θₒ::Float64
    v::Float64
    I::Float64
    ΔA::Float64
    τ::Float64
    η::Float64
    reflect::Bool
    discrete::Bool
    τ_Δv::Float64
end

struct _RaytracePixel
    ring::Int
    col::Int
    α::Float64
    β::Float64
    rCam::Float64
    ϕCam::Float64
    ΔA::Float64
end

struct _RaytraceGrid
    submodel::Int
    ringinds::Vector{Int}
    outinds::Vector{Int}
    rMin::Vector{Float64}
    rMax::Vector{Float64}
    Δϕ::Float64
    nϕ::Int
    localToPixel::Dict{Tuple{Int,Int},Int}
end

_rt_point_length(r::ring) = r.I isa AbstractVector ? length(r.I) : 1
_rt_point_value(x, col::Int) = x isa AbstractVector ? x[col] : x
_rt_wrap_π(ϕ::Float64) = atan(sin(ϕ), cos(ϕ))

function _rt_iratio_vector(IRatios::Float64, nSubmodels::Int)
    return fill(IRatios, nSubmodels)
end

function _rt_iratio_vector(IRatios::Array{Float64}, nSubmodels::Int)
    if ndims(IRatios) != 1 || length(IRatios) != nSubmodels
        throw(ArgumentError("IRatios must be a Float64 or a vector with one entry per submodel; got size $(size(IRatios)) for $nSubmodels submodels"))
    end
    return vec(IRatios)
end

function _rt_submodel_ring_range(m::model, submodel::Int)
    start = m.subModelStartInds[submodel]
    stop = submodel == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[submodel+1]-1
    return start:stop
end

function _rt_submodel_camera_range(camStartInds::Vector{Int}, m::model, submodel::Int)
    start = camStartInds[submodel]
    stop = submodel == length(camStartInds) ? length(m.camera.α) : camStartInds[submodel+1]-1
    return start:stop
end

function _rt_is_discrete(m::model, submodel::Int)
    return m.rings[m.subModelStartInds[submodel]].ϕ isa Float64
end

function _rt_submodel_rspan(m::model, camStartInds::Vector{Int}, submodel::Int)
    cr = _rt_submodel_camera_range(camStartInds, m, submodel)
    rCam = sqrt.(m.camera.α[cr].^2 .+ m.camera.β[cr].^2)
    return minimum(rCam), maximum(rCam)
end

function _rt_ring_edges(ring::ring, rCam::Float64)
    ΔrUp = ring.scale == :log ? rCam*(exp(ring.Δr)-1) : ring.Δr
    ΔrDown = ring.scale == :log ? min(rCam, rCam*(1-1/exp(ring.Δr))) : min(rCam, ring.Δr)
    return rCam - ΔrDown/2, rCam + ΔrUp/2
end

function _rt_ϕ_col(ϕ::Float64, Δϕ::Float64, nϕ::Int)
    shifted = mod(_rt_wrap_π(ϕ) - (-π - Δϕ/2), 2π)
    col = floor(Int, shifted/Δϕ) + 1
    return col > nϕ ? 1 : col
end

function _rt_make_vector!(r::ring, field::Symbol, n::Int)
    val = getfield(r, field)
    if !(val isa AbstractVector)
        setfield!(r, field, fill(val, n))
    end
end

function _rt_prepare_output_ring!(r::ring)
    n = _rt_point_length(r)
    for field in (:r, :i, :rot, :θₒ, :v, :I, :ϕ, :ϕ₀, :ΔA, :reflect, :τ, :η)
        _rt_make_vector!(r, field, n)
    end
    if isnothing(r.x) || !(r.x isa AbstractVector)
        r.x = [rotate3D_scalar(r.r[j], r.ϕ₀[j], r.i[j], r.rot[j], r.θₒ[j], r.reflect[j])[1] for j in 1:n]
    end
    if isnothing(r.y) || !(r.y isa AbstractVector)
        r.y = [rotate3D_scalar(r.r[j], r.ϕ₀[j], r.i[j], r.rot[j], r.θₒ[j], r.reflect[j])[2] for j in 1:n]
    end
    if isnothing(r.z) || !(r.z isa AbstractVector)
        r.z = [rotate3D_scalar(r.r[j], r.ϕ₀[j], r.i[j], r.rot[j], r.θₒ[j], r.reflect[j])[3] for j in 1:n]
    end
    r.I .= NaN
    return r
end

function _rt_set_point!(r::ring, col::Int, res)
    r.I[col] = res.I
    r.v[col] = res.v
    r.r[col] = res.r
    r.ϕ[col] = res.ϕ
    r.ϕ₀[col] = res.ϕ₀
    r.i[col] = res.i
    r.rot[col] = res.rot
    r.θₒ[col] = res.θₒ
    r.τ[col] = res.τ
    r.η[col] = res.η
    r.reflect[col] = res.reflect
    r.x[col] = res.x
    r.y[col] = res.y
    r.z[col] = res.z
end

function _rt_copy_cloud_point(p::_RaytracePoint, Ieff::Float64)
    return ring(r=p.r, i=p.i, rot=p.rot, θₒ=p.θₒ, v=p.v, I=Ieff, ϕ=p.ϕ, ϕ₀=p.ϕ₀,
        ΔA=p.ΔA, reflect=p.reflect, τ=p.τ, η=p.η, Δr=1.0, Δϕ=1.0, scale=nothing,
        x=p.x, y=p.α, z=p.β)
end

function _rt_mask_ring!(r::ring, mask::AbstractVector{Bool})
    for field in (:r, :i, :rot, :θₒ, :v, :I, :ϕ, :ϕ₀, :ΔA, :reflect, :τ, :η, :x, :y, :z)
        val = getfield(r, field)
        if val isa AbstractVector
            setfield!(r, field, val[mask])
        end
    end
    return r
end

function _rt_rebuild_camera(rings::Vector{ring})
    α = Float64[]
    β = Float64[]
    subStarts = [1]
    offset = 1
    while offset <= length(rings)
        nPer = _rt_point_length(rings[offset])
        stop = offset
        while stop < length(rings) && _rt_point_length(rings[stop+1]) == nPer
            stop += 1
        end
        for col in 1:nPer
            for ringInd in offset:stop
                r = rings[ringInd]
                push!(α, Float64(_rt_point_value(r.y, col)))
                push!(β, Float64(_rt_point_value(r.z, col)))
            end
        end
        offset = stop + 1
        offset <= length(rings) && push!(subStarts, offset)
    end
    return α, β, subStarts
end

function _rt_compact_rings(rings::Vector{ring})
    compact = ring[]
    for r in rings
        if r.I isa AbstractVector
            mask = isfinite.(r.I)
            any(mask) && push!(compact, _rt_mask_ring!(r, mask))
        elseif isfinite(r.I)
            push!(compact, r)
        end
    end
    return compact
end

function _rt_flatten_points(m::model, camStartInds::Vector{Int})
    points = _RaytracePoint[]
    sizehint!(points, length(m.camera.α))
    for s in 1:length(m.subModelStartInds)
        rr = _rt_submodel_ring_range(m, s)
        cr = _rt_submodel_camera_range(camStartInds, m, s)
        discrete = _rt_is_discrete(m, s)
        nRings = length(rr)
        nPer = _rt_point_length(m.rings[first(rr)])
        camOffset = first(cr) - 1
        for col in 1:nPer
            for (localRing, ringInd) in enumerate(rr)
                ring = m.rings[ringInd]
                camInd = camOffset + (col-1)*nRings + localRing
                α = Float64(m.camera.α[camInd])
                β = Float64(m.camera.β[camInd])
                rVal = Float64(_rt_point_value(ring.r, col))
                ϕ₀Val = Float64(_rt_point_value(ring.ϕ₀, col))
                iVal = Float64(_rt_point_value(ring.i, col))
                rotVal = Float64(_rt_point_value(ring.rot, col))
                θₒVal = Float64(_rt_point_value(ring.θₒ, col))
                reflectVal = Bool(_rt_point_value(ring.reflect, col))
                xVal = if !isnothing(ring.x)
                    Float64(_rt_point_value(ring.x, col))
                else
                    rotate3D_scalar(rVal, ϕ₀Val, iVal, rotVal, θₒVal, reflectVal)[1]
                end
                vVal = Float64(_rt_point_value(ring.v, col))
                τVal = Float64(_rt_point_value(ring.τ, col))
                τ_Δv = if ring.τ isa AbstractVector && nPer > 1
                    v2 = col == 1 ? Float64(_rt_point_value(ring.v, min(col+1, nPer))) : Float64(_rt_point_value(ring.v, col-1))
                    abs((vVal + v2)/2)
                else
                    Inf
                end
                push!(points, _RaytracePoint(
                    s, ringInd, col, α, β, hypot(α, β), atan(β, α), xVal, rVal,
                    Float64(_rt_point_value(ring.ϕ, col)), ϕ₀Val, iVal, rotVal, θₒVal, vVal,
                    Float64(_rt_point_value(ring.I, col)), Float64(_rt_point_value(ring.ΔA, col)),
                    τVal, Float64(_rt_point_value(ring.η, col)), reflectVal, discrete, τ_Δv))
            end
        end
    end
    return points
end

function _rt_build_output(m::model, camStartInds::Vector{Int})
    continuous = [s for s in 1:length(m.subModelStartInds) if !_rt_is_discrete(m, s)]
    sort!(continuous, by=s -> _rt_submodel_rspan(m, camStartInds, s)[1])
    grids = _RaytraceGrid[]
    pixels = _RaytracePixel[]
    outRings = ring[]
    αout = Float64[]
    βout = Float64[]
    subStarts = Int[]
    covered = Tuple{Vector{Float64},Vector{Float64}}[]

    for s in continuous
        rr = collect(_rt_submodel_ring_range(m, s))
        nRings = length(rr)
        nϕ = _rt_point_length(m.rings[rr[1]])
        camOffset = first(_rt_submodel_camera_range(camStartInds, m, s)) - 1
        selected = Int[]
        selectedLocal = Int[]
        rMinList = Float64[]
        rMaxList = Float64[]
        for (localRing, ringInd) in enumerate(rr)
            camInd = camOffset + localRing
            rCam = hypot(m.camera.α[camInd], m.camera.β[camInd])
            rMin, rMax = _rt_ring_edges(m.rings[ringInd], rCam)
            alreadyCovered = any(any((rCam .>= cov[1]) .& (rCam .<= cov[2])) for cov in covered)
            if !alreadyCovered
                push!(selected, ringInd)
                push!(selectedLocal, localRing)
                push!(rMinList, rMin)
                push!(rMaxList, rMax)
            end
        end
        isempty(selected) && continue
        push!(subStarts, length(outRings)+1)
        outinds = Int[]
        for ringInd in selected
            push!(outRings, _rt_prepare_output_ring!(deepcopy(m.rings[ringInd])))
            push!(outinds, length(outRings))
        end
        localToPixel = Dict{Tuple{Int,Int},Int}()
        for col in 1:nϕ
            for (idx, localRing) in enumerate(selectedLocal)
                ringInd = selected[idx]
                outRing = outinds[idx]
                camInd = camOffset + (col-1)*nRings + localRing
                push!(αout, Float64(m.camera.α[camInd]))
                push!(βout, Float64(m.camera.β[camInd]))
                pixelInd = length(pixels) + 1
                localToPixel[(idx, col)] = pixelInd
                push!(pixels, _RaytracePixel(outRing, col, αout[end], βout[end], hypot(αout[end], βout[end]),
                    atan(βout[end], αout[end]), Float64(_rt_point_value(m.rings[ringInd].ΔA, col))))
            end
        end
        push!(covered, (rMinList, rMaxList))
        push!(grids, _RaytraceGrid(s, selected, outinds, rMinList, rMaxList, m.rings[selected[1]].Δϕ, nϕ, localToPixel))
    end
    return grids, pixels, outRings, αout, βout, subStarts
end

function _rt_find_pixel(p::_RaytracePoint, grids::Vector{_RaytraceGrid})
    for grid in grids
        ringPos = findfirst(k -> p.rCam >= grid.rMin[k] && p.rCam <= grid.rMax[k], eachindex(grid.rMin))
        isnothing(ringPos) && continue
        col = _rt_ϕ_col(p.ϕCam, grid.Δϕ, grid.nϕ)
        pix = get(grid.localToPixel, (ringPos, col), 0)
        pix != 0 && return pix
    end
    return 0
end

function _rt_scan_bucket(points::Vector{_RaytracePoint}, inds::Vector{Int}, IRatios::Vector{Float64}, outputΔA::Float64, τCutOff::Float64)
    finiteInds = [idx for idx in inds if isfinite(points[idx].I) && isfinite(points[idx].x)]
    isempty(finiteInds) && return nothing
    order = sort(finiteInds, by=idx -> points[idx].x, rev=true)
    τvals = [points[idx].τ for idx in order]
    τ_Δv = [points[idx].τ_Δv for idx in order]
    if !(all(isinf, τ_Δv) || length(unique(τvals)) <= length(order))
        error("velocity-dependent optical depths not yet implemented -- pass τ as a float when creating models if you want to use raytracing")
    end
    weights = [points[idx].I * IRatios[points[idx].submodel] * points[idx].ΔA / outputΔA for idx in order]
    areas = [points[idx].ΔA for idx in order]
    firstPoint = points[order[1]]
    new_τ = firstPoint.τ
    new_I = weights[1]
    new_v = firstPoint.v * new_I
    new_r = firstPoint.r * new_I
    new_ϕ = firstPoint.ϕ * new_I
    new_ϕ₀ = firstPoint.ϕ₀ * new_I
    new_i = firstPoint.i * new_I
    new_rot = firstPoint.rot * new_I
    new_θₒ = firstPoint.θₒ * new_I
    new_η = firstPoint.η * new_I
    new_x = firstPoint.x * new_I
    new_y = firstPoint.α * new_I
    new_z = firstPoint.β * new_I
    new_reflect = firstPoint.reflect
    j = 2
    while j <= length(order)
        # fraction of the rear tile obscured by the largest tile in front of it:
        # obscuredFrac = A_rear / A_front_max. When >1 the rear tile is larger than
        # everything in front, so only A_front_max/A_rear = 1/obscuredFrac of it is
        # covered (attenuated); the rest peeks around at full brightness. Recomputed
        # per tile (not once at j=2).
        obscuredFrac = areas[j]/maximum(@view areas[1:j-1])
        (new_τ < τCutOff || obscuredFrac > 1.0) || break
        p = points[order[j]]
        if obscuredFrac > 1.0
            obscured_I = weights[j]/obscuredFrac          # covered part, attenuated
            unobscured_I = weights[j] - obscured_I        # peeking part, full strength
            tmp_I = exp(-new_τ)*obscured_I + unobscured_I
        else
            tmp_I = exp(-new_τ)*weights[j]
        end
        new_v += p.v * tmp_I
        new_r += p.r * tmp_I
        new_ϕ += p.ϕ * tmp_I
        new_ϕ₀ += p.ϕ₀ * tmp_I
        new_i += p.i * tmp_I
        new_rot += p.rot * tmp_I
        new_θₒ += p.θₒ * tmp_I
        new_η += p.η * tmp_I
        new_x += p.x * tmp_I
        new_y += p.α * tmp_I
        new_z += p.β * tmp_I
        new_I += tmp_I
        new_τ += p.τ
        j += 1
    end
    den = new_I == 0.0 ? 1.0 : new_I
    return (I=new_I, v=new_v/den, r=new_r/den, ϕ=new_ϕ/den, ϕ₀=new_ϕ₀/den,
        i=new_i/den, rot=new_rot/den, θₒ=new_θₒ/den, τ=new_τ, η=new_η/den,
        x=new_x/den, y=new_y/den, z=new_z/den, reflect=new_reflect)
end

function _rt_attenuate_free_clouds(points::Vector{_RaytracePoint}, freeInds::Vector{Int}, IRatios::Vector{Float64}, τCutOff::Float64)
    isempty(freeInds) && return ring[]
    radii = [sqrt(points[idx].ΔA/π) for idx in freeInds]
    cellSize = 2*maximum(radii)
    cellSize = cellSize == 0.0 ? 1.0 : cellSize
    cells = Dict{Tuple{Int,Int},Vector{Int}}()
    for (pos, idx) in enumerate(freeInds)
        p = points[idx]
        key = (floor(Int, p.α/cellSize), floor(Int, p.β/cellSize))
        push!(get!(cells, key, Int[]), pos)
    end
    out = ring[]
    for (pos, idx) in enumerate(freeInds)
        p = points[idx]
        key = (floor(Int, p.α/cellSize), floor(Int, p.β/cellSize))
        τfront = 0.0
        for di in -1:1, dj in -1:1
            for otherPos in get(cells, (key[1]+di, key[2]+dj), Int[])
                otherPos == pos && continue
                q = points[freeInds[otherPos]]
                dist = hypot(p.α - q.α, p.β - q.β)
                if dist < radii[pos] + radii[otherPos] && (q.x > p.x || (q.x == p.x && otherPos < pos))
                    τfront += q.τ
                end
            end
        end
        τfront > τCutOff && continue
        push!(out, _rt_copy_cloud_point(p, p.I * IRatios[p.submodel] * exp(-τfront)))
    end
    return out
end

function raytrace!(m::model;IRatios::Union{Float64,Array{Float64,}}=1.0,τCutOff::Float64=1.0,raytraceFreeClouds=false)
    if m.subModelStartInds == [1]
        @warn "raytrace! called on a model with no submodels -- maybe you already raytraced? Returning unaltered model."
        return m
    elseif m.camera.raytraced
        @warn "raytrace! called on a model that has already been raytraced -- returning unaltered model."
        return m
    end

    @info "raytracing model with $(length(m.subModelStartInds)) submodels"
    nSubmodels = length(m.subModelStartInds)
    IR = _rt_iratio_vector(IRatios, nSubmodels)
    camStartInds = getFlattenedCameraIndices(m)
    points = _rt_flatten_points(m, camStartInds)
    grids, pixels, outRings, αout, βout, subStarts = _rt_build_output(m, camStartInds)

    buckets = [Int[] for _ in pixels]
    freeInds = Int[]
    for (idx, p) in enumerate(points)
        pix = _rt_find_pixel(p, grids)
        if pix == 0
            p.discrete && isfinite(p.I) && push!(freeInds, idx)
        else
            push!(buckets[pix], idx)
        end
    end

    for (pixInd, bucket) in enumerate(buckets)
        res = _rt_scan_bucket(points, bucket, IR, pixels[pixInd].ΔA, τCutOff)
        isnothing(res) && continue
        pix = pixels[pixInd]
        _rt_set_point!(outRings[pix.ring], pix.col, res)
        outRings[pix.ring].y[pix.col] = pix.α
        outRings[pix.ring].z[pix.col] = pix.β
    end

    freeRings = if raytraceFreeClouds
        _rt_attenuate_free_clouds(points, freeInds, IR, τCutOff)
    else
        [_rt_copy_cloud_point(points[idx], points[idx].I * IR[points[idx].submodel]) for idx in freeInds]
    end
    if !isempty(freeRings)
        push!(subStarts, length(outRings)+1)
        for r in freeRings
            push!(outRings, r)
            xyz = getXYZ(r)
            push!(αout, xyz[2])
            push!(βout, xyz[3])
        end
    end
    outRings = _rt_compact_rings(outRings)
    if isempty(outRings)
        error("raytrace! produced an empty model")
    end
    αout, βout, subStarts = _rt_rebuild_camera(outRings)
    out = model(outRings, Dict{Symbol,profile}(), camera(αout, βout, true), subStarts)
    out.cache = Dict{Any,Array}()
    return out
end
