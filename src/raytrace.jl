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
    xyzSys = rotate3D(r,ϕ₀,i,rot,θₒPoint) #system coordinates xyz
    ϕ = atan(xyzSys[2],-xyzSys[1]) #ϕ after rotation, measured from +x in disk plane, -x because of how rotation matrix was implemented and desire to have ϕ=0 at +x
    return r, ϕ, ϕ₀
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
    cosr = cos(rot); sinr = sin(rot); cosi = cos(-i); sini = sin(-i); cosθₒ = cos(θₒPoint); sinθₒ = sin(θₒPoint) #flip i to match convention in rotate3D, +x bottom of disk
    xRing = (β*cosr - α*cosi*sinr)/(cosi*cosθₒ+cosr*sini*sinθₒ) #system x
    yRing = (α*(cosi*cosθₒ+sini/cosr*sinθₒ)+β*cosθₒ*sinr/cosr)/(cosi*cosθₒ/cosr+sini*sinθₒ)
    r = √(xRing^2 + yRing^2)
    ϕ₀ = atan(yRing,xRing) #original ϕ₀ (no rotation)
    xyz[1] = xRing; xyz[2] = yRing; xyz[3] = 0.0
    mul!(colBuff,r3D,xyz)
    undo_tilt = [sini 0.0 cosi; 0.0 1.0 0.0; -cosi 0.0 sini]
    mul!(xyz,undo_tilt,colBuff)
    ϕ = atan(xyz[2],-xyz[1]) #ϕ after rotation and being "puffed up", measured from +x in disk plane -- this is fine even for puffed up clouds but note ϕ is measured wrt to disk midplane then. -x because of how rotation matrix was implemented...
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
    xyzSys = rotate3D(r,ϕ₀,i,rot,θₒ,reflect)
    α = xyzSys[2] #camera is at +x, so α is y
    β = xyzSys[3] #and β is z
    return α, β
end

"""
    zeroDiskObscuredClouds!(m::model; diskCloudIntensityRatio::Float64=1.0, rotate3D::Function=rotate3D)

Zero out the intensities of clouds that are obscured by the disk.

Performs simple raytracing for an optically thick obscuring disk. The function
modifies the input model by setting the intensity of obscured cloud points to zero
and adjusting the disk intensity according to the specified ratio.

# Arguments
- `m::model`: Model to zero out disk obscured clouds. Should be a combined model consisting of a disk component and a cloud component. 
- `diskCloudIntensityRatio::Float64=1.0`: Ratio of disk to cloud intensity, used to scale 
  the disk intensities after zeroing out clouds
- `rotate3D::Function=rotate3D`: Function to rotate coordinates in 3D space

# Returns
- `m::model`: Model with disk obscured clouds zeroed out

# See also 
- `removeDiskObscuredClouds!`: Function to remove disk obscured clouds instead of zeroing them out
"""
function zeroDiskObscuredClouds!(m::model;diskCloudIntensityRatio::Float64=1.,rotate3D::Function=rotate3D)
    isCombined = length(m.subModelStartInds) > 1 #check if model is combined
    startInds = m.subModelStartInds #start indices of submodels
    if !isCombined
        @warn "did not detect combined model so nothing to zero -- returning unaltered input model"
        return m
    end
    if length(startInds) > 2
        #zero disk obscured clouds recursively combining two at a time
        mList = [deepcopy(m) for i=1:length(startInds)]
        for (i,m) in enumerate(mList)
            s = startInds[i]; e = i == length(startInds) ? length(m.rings) : startInds[i+1]-1
            m.rings = m.rings[s:e]
        end
        mFinal = mList[1]
        for i=2:length(mList)
            mFinal = mFinal + mList[i]
            mFinal = zeroDiskObscuredClouds!(mFinal,rotate3D)
        end
        m = mFinal
        return mFinal
    end
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
        xyzCloud = rotate3D(ring.r,ring.ϕ₀,ring.i,ring.rot,ring.θₒ,ring.reflect) #system coordinates xyz
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
    removeDiskObscuredClouds!(m::model, rotate3D::Function=rotate3D)

Remove clouds that are obscured by the disk.

Performs simple raytracing for an optically thick obscuring disk. The function
modifies the input model by removing cloud points that are obscured by the disk.
Note that this is a mutating operation and the input model will be modified in place.

# Arguments
- `m::model`: Model to remove disk obscured clouds. Should be a combined model consisting 
  of a disk component and a cloud component.
- `rotate3D::Function=rotate3D`: Function to rotate coordinates in 3D space

# Returns
- `m::model`: Model with disk obscured clouds removed

# See also 
- `zeroDiskObscuredClouds!`: Function to zero out disk obscured clouds instead of removing them
"""
function removeDiskObscuredClouds!(m::model,rotate3D::Function=rotate3D)
    isCombined = length(m.subModelStartInds) > 1 #check if model is combined
    startInds = m.subModelStartInds #start indices of submodels
    if !isCombined
        @warn "did not detect combined model so nothing to remove -- returning unaltered input model"
        return m
    end
    if length(startInds) > 2
        #remove disk obscured clouds recursively combining two at a time
        mList = [deepcopy(m) for i=1:length(startInds)]
        for (i,m) in enumerate(mList)
            s = startInds[i]; e = i == length(startInds) ? length(m.rings) : startInds[i+1]-1
            m.rings = m.rings[s:e]
        end
        mFinal = mList[1]
        for i=2:length(mList)
            mFinal = mFinal + mList[i]
            mFinal = removeDiskObscuredClouds!(mFinal,rotate3D)
        end
        return mFinal
    end
    removeFlag = falses(length(m.rings))
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
            diskFlag[ind:ind+stride] .= true
            ind += stride
        end
    end
    diskα = m.camera.α[diskFlag]; diskβ = m.camera.β[diskFlag]
    rDisk = sqrt.(diskα.^2 .+ diskβ.^2)
    rMinDisk = minimum(rDisk); rMaxDisk = maximum(rDisk)
    rUnique = unique(r -> round(r,sigdigits=12),sqrt.(m.camera.α.^2 .+ m.camera.β.^2))
    ϕList = [m.rings[i].ϕ for i in 1:length(m.rings)]
    cloudRingInds = [i for i in 1:length(m.rings) if typeof(m.rings[i].ϕ) == Float64 && typeof(m.rings[i].r) == Float64]
    for i in cloudRingInds #check if r inside rMin/max, then find closest disk cell, then compare xs and flag for removal if xCloud behind xDisk
        xyzSys = rotate3D(m.rings[i].r,m.rings[i].ϕ₀,m.rings[i].i,m.rings[i].rot,m.rings[i].θₒ,m.rings[i].reflect) #system coordinates xyz
        xCloud = xyzSys[1]
        α,β = m.camera.α[αβStartInds[i]],m.camera.β[αβStartInds[i]]
        rCam = sqrt(α^2 + β^2)
        ϕCam = atan(β,α)
        if (rCam > rMinDisk) && (rCam < rMaxDisk)
            rDiskInd = argmin(abs.(rCam .- rUnique))
            ϕDiskInd = argmin(abs.(ϕCam .- ϕList[rDiskInd]))
            xDisk = rotate3D(m.rings[rDiskInd].r[ϕDiskInd],m.rings[rDiskInd].ϕ₀[ϕDiskInd],m.rings[rDiskInd].i,m.rings[rDiskInd].rot,m.rings[rDiskInd].θₒ,m.rings[rDiskInd].reflect)[1] #system coordinates xyz
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
- `IRatios::Union{Float64,Array{Float64,}}=1.0`: Intensity ratios for each submodel
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
function raytrace!(m::model;IRatios::Union{Float64,Array{Float64,}}=1.0,τCutOff::Float64=1.0,raytraceFreeClouds=false)
    if m.subModelStartInds == [1]
        @warn "raytrace! called on a model with no submodels -- maybe you already raytraced? Returning unaltered model."
        return m
    elseif m.camera.raytraced
        @warn "raytrace! called on a model that has already been raytraced -- returning unaltered model."
        return m
    else
        @info "raytracing model with $(length(m.subModelStartInds)) submodels"
        #1 model with coarsest grid becomes the new "base" model -- if tie pick the one with most points? 
        #2 any models with smaller grid are added to corresponding larger grid model, if they are outside of the larger grid keep original small grid in those places and put them to back of array + set startInds to still show combined
        #3 work from +x camera backwards until τ ~ 1, adding intensity to each pixel as needed
        #4 if τ is fixed assume it applies at all wavelengths, else if τ is a function of wavelength check for spatial and wavelength overlap
        #5 remove extra points from final model if they are obscured 

        #1 detect grid sizes, ranges, types
        gridInfo = Array{Tuple{Float64,Float64,Int}}(undef,length(m.subModelStartInds)) #contains (rMin,rMax,gridType) for each submodel, where gridType is 0 for continuous (i.e. disk) and 1 for discrete (i.e. cloud)
        camStartInds = getFlattenedCameraIndices(m) #flattened camera indices for each submodel
        for i in 1:length(m.subModelStartInds) 
            discrete = typeof(m.rings[m.subModelStartInds[i]].ϕ) == Float64 #if cloud/discrete model, each ring contains just a single point, otherwise each ring contains a grid of points
            camStartInd = camStartInds[i] #flattened camera index for this submodel
            camEndInd = i == length(m.subModelStartInds) ? length(m.camera.α) : camStartInds[i+1]-1 #flattened camera index for the next submodel
            αSegment = m.camera.α[camStartInd:camEndInd]
            βSegment = m.camera.β[camStartInd:camEndInd]

            rCamSegment = sqrt.(αSegment.^2 .+ βSegment.^2) #camera coordinates
            rMin = minimum(rCamSegment)
            rMax = maximum(rCamSegment)
            if discrete
                gridInfo[i] = (rMin,rMax,1)
            else
                gridInfo[i] = (rMin,rMax,0)
            end
        end
        contCounter = sum([gridInfo[i][3]!=1 for i in 1:length(m.subModelStartInds)]) #number of continuous models
        overlaps = []
        baseModel = nothing #initialize
        if contCounter >= 1
            for i in 1:length(m.subModelStartInds)
                if gridInfo[i][3] == 0
                    #check if this model is inside any of the continuous models
                    for j in 1:length(m.subModelStartInds)
                        if gridInfo[j][3] == 0 && i != j
                            if (gridInfo[i][1] >= gridInfo[j][1] && gridInfo[i][1] <= gridInfo[j][2]) || (gridInfo[i][2] <= gridInfo[j][2] && gridInfo[i][2] >= gridInfo[j][1])
                                #continuous models have some overlap -- whose grid do we use?
                                if (j,i) ∉ overlaps #only add unique pairs
                                    push!(overlaps, (i,j))
                                end
                            end
                        end
                    end
                end
            end
            #2 combine models if necessary
            baseModelInd = 1
            if length(overlaps) > 0 # set base model ind to be the one with smallest rMin
                innerInd = 0 
                rMin = Inf
                for (i,j) in overlaps
                    if gridInfo[i][1] < gridInfo[j][1] && gridInfo[i][1] < rMin
                        innerInd = i
                        rMin = gridInfo[i][1]
                    elseif gridInfo[j][1] < gridInfo[i][1] && gridInfo[j][1] < rMin
                        innerInd = j
                        rMin = gridInfo[j][1]
                    end
                end
                baseModelInd = innerInd
            end
            baseModelRings = deepcopy(m.rings[m.subModelStartInds[baseModelInd]:m.subModelStartInds[baseModelInd+1]-1])
            camStartInd = camStartInds[baseModelInd] #flattened camera index for this submodel
            camEndInd = baseModelInd == length(m.subModelStartInds) ? length(m.camera.α) : camStartInds[baseModelInd+1]-1 #flattened camera index for the next submodel
            baseModelα = deepcopy(m.camera.α[camStartInd:camEndInd])
            baseModelβ = deepcopy(m.camera.β[camStartInd:camEndInd])
            baseModel_rCam = sqrt.(baseModelα.^2 .+ baseModelβ.^2) #camera coordinates
            baseModel_rUnique = unique(r -> round(r,sigdigits=9),baseModel_rCam) #unique radii in camera coordinates
            baseModel_ϕCam = atan.(baseModelβ,baseModelα) #camera coordinates
            
            baseCam = camera(baseModelα,baseModelβ,false) #camera object for base model
            baseModel = model(baseModelRings,Dict{Symbol,profile}(),baseCam,[1]) #initialize new model object with base model rings
            base_r = getVariable(baseModel,:r,flatten=true)
            base_i = getVariable(baseModel,:i,flatten=true)
            base_ϕ = getVariable(baseModel,:ϕ,flatten=true)
            base_ϕ₀ = getVariable(baseModel,:ϕ₀,flatten=true)    
            base_v = getVariable(baseModel,:v,flatten=true)
            base_I = getVariable(baseModel,:I,flatten=true)
            base_τ = getVariable(baseModel,:τ,flatten=true)
            base_η = getVariable(baseModel,:η,flatten=true)
            
            #plan: go through each ring of the base model and check all submodels for their overlaps, push to list
            #then raytrace from +x backwards until τ ~ 1, adding intensity to each pixel as needed in base model, with intensity ratios if supplied 
            #for indices that get added to base model, mark them for removal by setting I = 0
            #after going through base model, check if there are any other continuous models that overlap with each other (not base model) and repeat this process
            #when all overlaps have been checked, remove all rings and camera points that are marked for removal
            #create new model with combined rings and camera points, going from smallest to largest rMin 
            for i in 1:length(baseModel_rCam)
                #check if any other submodels overlap with our base model between this ring and the next
                #rCam and ϕCam are the centers of camera pixels, use Δr and Δϕ to calculate edges
                rCam = baseModel_rCam[i]
                ϕCam = baseModel_ϕCam[i]
                Δr_base = getVariable(baseModel,:Δr)
                Δϕ_base = getVariable(baseModel,:Δϕ)
                ring,col = getRingFromFlattenedInd(baseModel,i) #get ring and column for this pixel     
                ΔrUp = baseModelRings[ring].scale == :log ? rCam*(exp(Δr_base[ring])-1) : Δr_base[ring] #log scale
                ΔrDown = m.rings[ring].scale == :log ? minimum([rCam,rCam*(1-1/exp(Δr_base[ring]))]) : minimum([rCam,Δr_base[ring]]) #if first ring, then use scale to calculate what "should be" the next inner ring and take minimum of that and the current radius (so don't go through zero)
                Δϕ = Δϕ_base[ring] #always linear in ϕ
                rMin = rCam - ΔrDown/2
                rMax = rCam + ΔrUp/2
                ϕMin = atan(sin(ϕCam - Δϕ/2), cos(ϕCam - Δϕ/2)) #azimuthal angle of lower boundary of pixel
                ϕMax = atan(sin(ϕCam + Δϕ/2), cos(ϕCam + Δϕ/2)) #correct for periodicity / make sure angular coordinates are the same   
                i_i = length(base_i) < length(base_I) ? baseModelRings[ring].i : base_i[i] #inclination angle of this point
                θₒi = baseModelRings[ring].θₒ #opening angle of this point
                reflecti = baseModelRings[ring].reflect #reflection flag of this point
                roti = (length(baseModelRings[ring].rot) < length(baseModelRings[ring].I) || length(baseModelRings[ring].I) == 1) ? baseModelRings[ring].rot : baseModelRings[ring].rot[i] #rotation of this point
                xList = [rotate3D(base_r[i],base_ϕ₀[i],i_i,roti,θₒi,reflecti)[1]] #x coordinate in physical space of pixel
                IList = [base_I[i]]; vList = [base_v[i]]; ϕList = [base_ϕ[i]]; rList = [base_r[i]]; ΔAList = [baseModelRings[ring].ΔA[col]]; #initialize lists for intensities, velocities, azimuthal angles, radii, and areas
                iList = length(base_i) < length(base_I) ? [baseModelRings[ring].i] : [base_i[i]]; 
                τList = length(base_τ) < length(base_I) ? [baseModelRings[ring].τ] : [base_τ[i]]; ηList = length(base_η) < length(base_I) ? [baseModelRings[ring].η] : [base_η[i]]; #initialize lists for intensities, velocities, azimuthal angles, radii, inclinations, and optical depths
                Δv = i == 1 ? (base_v[i]+base_v[i+1])/2 : (base_v[i]+base_v[i-1])/2 #calculate Δv as average of neighboring points
                if length(base_τ) < length(base_I)
                    τ_ΔvList = [Inf] 
                else
                    τ_ΔvList = [abs(Δv)]
                end
                #check if any other submodels have points within this pixel
                for j in 1:length(m.subModelStartInds) 
                    if j != baseModelInd
                        endInd = j == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[j+1]-1 #end index of submodel
                        subModel = model(m.rings[m.subModelStartInds[j]:endInd],nothing,nothing,[j]) #submodel to check
                        camStartInd = camStartInds[j] #flattened camera index for this submodel
                        camEndInd = j == length(m.subModelStartInds) ? length(m.camera.α) : camStartInds[j+1]-1 #flattened camera index for the next submodel
                        αSegment = m.camera.α[camStartInd:camEndInd]
                        βSegment = m.camera.β[camStartInd:camEndInd]
                        rCamSegment = sqrt.(αSegment.^2 .+ βSegment.^2) #camera coordinates
                        ϕCamSegment = atan.(βSegment,αSegment) #camera coordinates
                        sub_r = getVariable(subModel,:r,flatten=true)
                        sub_ϕ = getVariable(subModel,:ϕ,flatten=true)
                        sub_v = getVariable(subModel,:v,flatten=true)
                        sub_I = getVariable(subModel,:I,flatten=true)
                        sub_θₒ = getVariable(subModel,:θₒ,flatten=true)
                        sub_reflect = getVariable(subModel,:reflect,flatten=true)
                        sub_ϕ₀ = getVariable(subModel,:ϕ₀,flatten=true)
                        sub_i = getVariable(subModel,:i,flatten=true)
                        sub_rot = getVariable(subModel,:rot,flatten=true)
                        sub_τ = getVariable(subModel,:τ,flatten=true)
                        sub_η = getVariable(subModel,:η,flatten=true)
                        #check if any of the camera points are within this pixel                        
                        for k in 1:length(rCamSegment)
                            ϕMin_k = copy(ϕMin); ϕMax_k = copy(ϕMax) #copy to avoid modifying original
                            rCam_k = rCamSegment[k]
                            ϕCam_k = ϕCamSegment[k]
                            #check if this camera point is within the pixel
                            if (αSegment[k] < 0.0) && (sign(ϕMin)!=sign(ϕMax)) && (baseModelα[i] < 0.0) #at ±π boundary in both submodel and base model
                                if sign(ϕCam_k) == 1 #top left qudrant, need to change sign of ϕMax 
                                    ϕMax_k = ϕMax + 2π #now ϕMax and ϕMin are both positive, with ϕMax extended across the ±π boundary
                                else #bottom right quadrant, need to change sign of ϕMin
                                    ϕMin_k = ϕMin - 2π #now ϕMax and ϕMin are both negative, with ϕMin extended across the ±π boundary
                                end
                            end
                            if (rCam_k >= rMin && rCam_k <= rMax) && (ϕCam_k >= ϕMin_k && ϕCam_k <= ϕMax_k) 
                                subRing, subCol = getRingFromFlattenedInd(subModel,k)
                                #check if this camera point is within the ring
                                discrete = gridInfo[j][3] == 1 #if cloud/discrete model, each ring contains just a single point, otherwise each ring contains a grid of points
                                sub_ik = length(sub_i) < length(sub_I) ? subModel.rings[subRing].i : sub_i[k] #inclination angle of this point
                                sub_θₒk = length(sub_θₒ) < length(sub_I) ? subModel.rings[subRing].θₒ : sub_θₒ[k] #opening angle of this point
                                sub_reflectk = typeof(sub_reflect[subRing]) == Float64 ? Bool(sub_reflect[subRing]) : sub_reflect[subRing] #reflection flag of this point
                                sub_rotk = length(sub_rot) < length(sub_I) ? subModel.rings[subRing].rot : sub_rot[k] #rotation of this point
                                x = rotate3D(sub_r[k],sub_ϕ₀[k],sub_ik,sub_rotk,sub_θₒk,sub_reflectk)[1]  #x coordinate in physical space of pixel
                                push!(xList,x) #add to list of x coordinates
                                ΔARatio = subModel.rings[subRing].ΔA[subCol] / baseModelRings[ring].ΔA[col] #ratio of area of this point to area of base model point
                                Itmp = typeof(IRatios) == Float64 ? sub_I[k]*IRatios*ΔARatio : sub_I[k]*IRatios[j]*ΔARatio #apply intensity ratio if supplied, regardless apply ratio based on area change
                                push!(IList,Itmp) #add to list of intensities
                                push!(vList,sub_v[k]) #add to list of velocities
                                push!(ϕList,sub_ϕ[k]) #add to list of azimuthal angles
                                push!(rList,sub_r[k]) #add to list of radii
                                push!(ΔAList,subModel.rings[subRing].ΔA[subCol]) #add to list of areas

                                #i, τ, and η are at minimum defined for each ring, but not always defined for each point in the ring, so check and add accordingly
                                if length(sub_i) < length(sub_I)
                                    push!(iList,subModel.rings[subRing].i) #add to list of inclinations
                                else
                                    push!(iList,sub_i[k]) #add to list of inclinations
                                end
                                if length(sub_τ) < length(sub_I)
                                    push!(τList,subModel.rings[subRing].τ) #add to list of optical depths
                                    push!(τ_ΔvList,Inf)
                                else
                                    push!(τList,sub_τ[k]) #add to list of optical depths
                                    Δv = k == 1 ? (sub_v[k]+sub_v[k+1])/2 : (sub_v[k]+sub_v[k-1])/2 #calculate Δv as average of neighboring points
                                    push!(τ_ΔvList,abs(Δv))
                                end
                                if length(sub_η) < length(sub_I)
                                    push!(ηList,subModel.rings[subRing].η) #add to list of η values
                                else
                                    push!(ηList,sub_η[k]) #add to list of η values
                                end
                                if typeof(m.rings[m.subModelStartInds[j]+subRing-1].I) == Float64 #single point, can't index
                                    m.rings[m.subModelStartInds[j]+subRing-1].I = NaN 
                                else
                                    m.rings[m.subModelStartInds[j]+subRing-1].I[subCol] = NaN #mark for removal -- all I = 0.0 points removed after pass
                                end                            
                            end
                        end
                    end
                end
                #lists now contain all points in submodels that overlap with basemodel pixel -- combine them!
                if length(IList) > 1 #any extra points?
                    tmpxList = [isnan(xList[i]) ? -Inf : xList[i] for i in 1:length(xList)] #replace NaNs with -Inf so they are sorted to the end
                    order = sortperm(tmpxList,rev=true) #sort by x coordinate, +x is closest to camera
                    NaNList = isnan.(IList) #NaN flag = do nothing, marked for removal
                    start = findfirst(!,NaNList) #find first point where we need to do something
                    if !isnothing(start) #if there are any points we need to combine
                        new_τ = τList[order][1]
                        new_I = IList[order][1]; new_v = vList[order][1]*new_I; new_ϕ = ϕList[order][1]*new_I; new_r = rList[order][1]*new_I; new_i = iList[order][1]*new_I; new_η = ηList[order][1]*new_I #take first point as the "base" point 
                        j = 2
                        #number of unique possible τ values is number of continuous models + number of clouds 
                        nPossibleUnique_τ = sum(gridInfo[i][3] == 0 for i in 1:length(m.subModelStartInds))
                        for k in 1:length(m.subModelStartInds)
                            if gridInfo[k][3] == 1
                                endInd = k == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[k+1]-1 #end index of submodel
                                nPossibleUnique_τ += k == length(m.subModelStartInds) ? length(m.rings)-m.subModelStartInds[k] : m.subModelStartInds[k+1]-m.subModelStartInds[k] #add number of discrete points in this submodel
                            end
                        end
                        if all(τ_ΔvList .== Inf) || length(unique(τList)) <= nPossibleUnique_τ #if all τ are the same or τ is a single value, we can just add points until τ ~ 1
                            obscuredFrac = ΔAList[order][j]/maximum(ΔAList[order][1:j-1]) #fraction of area of this tile obscured by those in front of it, i.e. how much of the point is visible -- if < 1 fully covered by previous tile(s)
                            while (new_τ < τCutOff || obscuredFrac > 1.0) && j <= length(order) #3: add points until τ ~ 1 
                                if !isnan(IList[order][j]) && !isnan(xList[order][j]) #skip points that are marked for removal or unphysical points (r = NaN)
                                    if obscuredFrac > 1.0
                                        obscured_I = IList[order][j]/obscuredFrac
                                        unobscured_I = IList[order][j]-obscured_I
                                        tmp_I = exp(-new_τ)*obscured_I #add the part that is obscured
                                        tmp_I += IList[order][j]*unobscured_I #add unboscured intensity from previous point without attenuation 
                                    else
                                        tmp_I = exp(-new_τ)*IList[order][j] #assume new intensity is added to previous one in chunks, no absorption in between 
                                    end
                                    new_v += vList[order][j]*tmp_I 
                                    new_r += rList[order][j]*tmp_I 
                                    new_ϕ += ϕList[order][j]*tmp_I
                                    new_i += iList[order][j]*tmp_I 
                                    new_η += ηList[order][j]*tmp_I
                                    new_I += tmp_I #add intensity
                                    new_τ += τList[order][j] #new optical depth after adding this point 
                                end
                                j += 1
                            end
                            den = new_I == 0.0 ? 1.0 : new_I #avoid division by zero
                            new_v /= den #average velocity
                            new_r /= den #average radius
                            new_ϕ /= den #average azimuthal angle
                            new_i /= den #average inclination
                            new_η /= den #average η
                        else #4: think about how to do this properly...
                            error("velocity-dependent optical depths not yet implemented -- pass τ as a float when creating models if you want to use raytracing")
                        end      
                        baseModel.rings[ring].I[col] = new_I #set intensity
                        baseModel.rings[ring].v[col] = new_v #set velocity
                        baseModel.rings[ring].ϕ[col] = new_ϕ #set azimuthal angle
                        baseModel.rings[ring].r[col] = new_r #set radius
                        if typeof(baseModel.rings[ring].i) == Float64 #if i is a single value, need to fill all rings with vectors then replace point with new_i
                            for r in 1:length(baseModel.rings)
                                baseModel.rings[r].i = fill(baseModel.rings[r].i,length(baseModel.rings[r].I)) #fill with vector of same length as I
                            end
                            base_i = getVariable(baseModel,:i,flatten=true) #get new base_i after filling
                        end
                        baseModel.rings[ring].i[col] = new_i #set inclination
                        if typeof(baseModel.rings[ring].τ) == Float64 #if τ is a single value, fill with vector then replace point with new_τ
                            for r in 1:length(baseModel.rings)
                                baseModel.rings[r].τ = fill(baseModel.rings[r].τ,length(baseModel.rings[r].I)) #fill with vector of same length as I
                            end
                            base_τ = getVariable(baseModel,:τ,flatten=true) #get new base_τ after filling
                        end
                        baseModel.rings[ring].τ[col] = new_τ #set optical depth
                        if typeof(baseModel.rings[ring].η) == Float64 #if η is a single value, fill with vector then replace point with new_η
                            for r in 1:length(baseModel.rings)
                                baseModel.rings[r].η = fill(baseModel.rings[r].η,length(baseModel.rings[r].I)) #fill with vector of same length as I
                            end
                            base_η = getVariable(baseModel,:η,flatten=true) #get new base_η after filling
                        end
                        baseModel.rings[ring].η[col] = new_η #set η
                    end
                end
            end
            #set all base model rings in original model to NaN so they can be removed later
            endInd = baseModelInd == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[baseModelInd+1]-1 #ring index for the next submodel
            for ring in m.rings[m.subModelStartInds[baseModelInd]:endInd]
                ring.I .= NaN #mark for removal
            end
            #call raytrace again with new "model" that is sum of all models with overlaps but not overlapping base model
            if length(overlaps) > 1 
                for (i,j) in overlaps #both continuous models 
                    #check if i or j is the base model, if so skip it
                    if i != baseModelInd && j != baseModelInd
                    #i is new "base" model, check for all other models that overlap with it 
                        endInd = i = length(m.subModelStartInds) ? length(m.camera.α) : camStartInds[i+1]-1 #flattened camera index for the next submodel
                        camera_i = camera(m.camera.α[camStartInds[i]:endInd],m.camera.β[camStartInds[i]:endInd],false) #camera α coordinates for submodel i
                        endInd = j == length(m.subModelStartInds) ? length(m.camera.α) : camStartInds[j+1]-1 #flattened camera index for the next submodel
                        camera_j = camera(m.camera.α[camStartInds[j]:endInd],m.camera.β[camStartInds[j]:endInd],false) #camera α coordinates for submodel j
                        endInd = i == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[i+1]-1 #ring index for the next submodel
                        subModelI = model(m.rings[m.subModelStartInds[i]:endInd],nothing,camera_i,[1]) #submodel i to check
                        endInd = j == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[j+1]-1 #ring index for the next submodel
                        subModelI += model(m.rings[m.subModelStartInds[j]:endInd],nothing,camera_j,[1]) #add j to check
                        IRatiosTmp = deepcopy(IRatios[i,j])
                        for ii in 1:length(m.subModelStartInds) #check all submodels for overlaps with submodel i
                            if ii != i && ii != j #not same pair
                                if ii != baseModelInd #not base model
                                    endInd = ii = length(m.subModelStartInds) ? length(m.camera.α) : camStartInds[ii+1]-1 #flattened camera index for the next submodel
                                    camera_ii = camera(m.camera.α[camStartInds[ii]:endInd],m.camera.β[camStartInds[ii]:endInd],false) #camera α coordinates for submodel ii
                                    endInd = ii == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[ii+1]-1 #ring index for the next submodel
                                    if gridInfo[ii][3] == 1 #discrete model
                                        rMin = minimum([gridInfo[i][1],gridInfo[j][1]]) #minimum radius of continuous models i and j
                                        rMax = maximum([gridInfo[i][2],gridInfo[j][2]]) #maximum radius of continuous models i and j
                                        subModel_ii = model(m.rings[m.subModelStartInds[ii]:endInd],nothing,camera_ii,[1]) #submodel ii to check
                                        for (jj,ring) in enumerate(subModel_ii.rings)
                                            rtmp = sqrt(subModel_ii.camera.α[jj]^2 + subModel_ii.camera.β[jj]^2) #camera coordinates
                                            if rtmp <= rMin || rtmp >= rMax #free floating cloud, mark for "deletion" (still in original model m, but don't want to add duplicates in this step)
                                                subModel_ii.rings[jj].I .= NaN #mark for removal
                                            else
                                                original_model_index = m.subModelStartInds[ii] + jj - 1 #original model index
                                                m.rings[original_model_index].I .= NaN #mark original model index for removal, because it will be added to subModelI and then the new base model
                                            end
                                        end
                                    end
                                    subModelI += subModel_ii #add submodel ii to submodel i
                                    push!(IRatios,IRatiosTmp[jj]) #add multiplier for jj
                                end
                            end
                        end
                        subModelI = raytrace!(subModelI;IRatios=IRatiosTmp,τCutOff=τCutOff) #raytrace submodel i
                        #remove NaN points from submodel i
                        subModelI = removeNaN!(subModelI) 
                        #add submodel i to base model
                        baseModel += subModelI #baseModel now has every continuous model + overlapping clouds, just need to add back in non-overlapping clouds
                    end
                end
            end
        end

        #raytrace discrete models that don't overlap with any continuous models (last thing left to check)
        nLeft = 0; firstSubModel = 0
        for i in 1:length(m.subModelStartInds)
            if gridInfo[i][3] == 1 
                if firstSubModel == 0 #first discrete model found, set firstSubModel to i
                    firstSubModel = i
                end
                nLeft += 1 #count number of discrete models left
            end
        end
        if nLeft > 0
            endInd = firstSubModel == length(m.subModelStartInds) ? length(m.camera.α) : camStartInds[firstSubModel+1]-1 #flattened camera index for the next submodel
            camtmp = camera(m.camera.α[camStartInds[firstSubModel]:endInd],m.camera.β[camStartInds[firstSubModel]:endInd],false) #camera α coordinates for first discrete model
            endInd = firstSubModel == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[firstSubModel+1]-1 #ring index for the next submodel
            subModel = model(m.rings[m.subModelStartInds[firstSubModel]:endInd],nothing,camtmp,[1]) #submodel to check
            for i in firstSubModel:length(m.subModelStartInds) 
                if (gridInfo[i][3] == 1) && (i > firstSubModel) #not first discrete model, check for overlaps
                    endInd = i == length(m.subModelStartInds) ? length(m.camera.α) : camStartInds[i+1]-1 #flattened camera index for the next submodel
                    camtmp = camera(m.camera.α[camStartInds[i]:endInd],m.camera.β[camStartInds[i]:endInd],false) #camera α coordinates for submodel i
                    endInd = i == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[i+1]-1 #ring index for the next submodel
                    subModel += model(m.rings[m.subModelStartInds[i]:endInd],nothing,camtmp,[1]) #add to list of discrete models 
                end
            end
            if raytraceFreeClouds
                αSubModel = subModel.camera.α
                βSubModel = subModel.camera.β
                ΔASubModel = getVariable(subModel,:ΔA)
                ΔrSubModel = sqrt.(ΔASubModel./π) #assume circular clouds on camera
                sub_r = getVariable(subModel,:r)
                sub_ϕ = getVariable(subModel,:ϕ)
                sub_v = getVariable(subModel,:v)
                sub_I = getVariable(subModel,:I)
                sub_θₒ = getVariable(subModel,:θₒ)
                sub_reflect = getVariable(subModel,:reflect)
                sub_ϕ₀ = getVariable(subModel,:ϕ₀)
                sub_i = getVariable(subModel,:i)
                sub_rot = getVariable(subModel,:rot)
                sub_τ = getVariable(subModel,:τ)
                sub_η = getVariable(subModel,:η)

                for i in 1:length(αSubModel)
                    αi = αSubModel[i]
                    βi = βSubModel[i]
                    Δri = ΔrSubModel[i]
                    xList = [rotate3D(sub_r[i],sub_ϕ₀[i],sub_i[i],sub_rot[i],sub_θₒ[i],sub_reflect[i])[1]] #x coordinate in physical space of pixel
                    IList = [sub_I[i]]; vList = [sub_v[i]]; ϕList = [sub_ϕ[i]]; rList = [sub_r[i]]; iList = [sub_i[i]]; τList = [sub_τ[i]]; ηList = [sub_η[i]] #initialize lists for intensities, velocities, azimuthal angles, radii, inclinations, and optical depths
                    for j in i:length(αSubModel)
                        αj = αSubModel[j]
                        βj = βSubModel[j]
                        Δrj = ΔrSubModel[j]
                        dist = sqrt((αi-αj)^2 + (βi-βj)^2) #distance between points
                        if dist < (Δri + Δrj) #are they touching?
                            #calculate new α, β, and I
                            push!(xList,rotate3D(sub_r[j],sub_ϕ₀[j],sub_i[j],sub_rot[j],sub_θₒ[j],sub_reflect[j])[1]) #add x coordinate of j to list
                            ΔARatio = subModel.rings[j].ΔA[1] / subModel.rings[i].ΔA[1] #ratio of area of this point to area of base model point
                            Itmp = sub_I[j]*IRatios*ΔARatio #apply intensity ratio if supplied, regardless apply ratio based on area change
                            push!(IList,Itmp) #add to list of intensities
                            push!(vList,sub_v[j]) #add velocity of j to list
                            push!(ϕList,sub_ϕ[j]) #add azimuthal angle of j to list
                            push!(rList,sub_r[j]) #add radius of j to list
                            push!(iList,sub_i[j]) #add inclination of j to list
                            push!(τList,sub_τ[j]) #add optical depth of j to list
                            push!(ηList,sub_η[j]) #add η of j to list
                            #mark for deletion
                            subModel.rings[i].I[1] = NaN #mark for removal -- all I = 0.0 points removed after pass
                        end
                    end
                    if length(xList) > 1 
                        order = sortperm(xList,rev=true) #sort by x coordinate, +x is closest to camera
                        NaNList = isnan.(IList) #NaN flag = do nothing, marked for removal
                        start = findfirst(!,NaNList) #find first point where we need to do something
                        if !isnothing(start) 
                            new_τ = τList[order][1]
                            new_I = IList[order][1]; new_v = vList[order][1]*new_I; new_ϕ = ϕList[order][1]*new_I; new_r = rList[order][1]*new_I; new_i = iList[order][1]*new_I; new_η = ηList[order][1]*new_I #take first point as the "base" point 
                            j = 2
                            while new_τ < τCutOff && j <= length(order) #3: add points until τ ~ 1
                                if !isnan(IList[order][j]) && !isnan(xList[order][j]) #skip points that are marked for removal
                                    tmp_I = exp(-new_τ)*IList[order][j] #assume new intensity is added to previous one in chunks, no absorption 
                                    new_v += vList[order][j]*tmp_I 
                                    new_r += rList[order][j]*tmp_I 
                                    new_ϕ += ϕList[order][j]*tmp_I
                                    new_i += iList[order][j]*tmp_I 
                                    new_η += ηList[order][j]*tmp_I
                                    new_I += tmp_I #add intensity
                                    new_τ += τlist[order][j] #new optical depth after adding this point 
                                end
                                j += 1
                            end
                            new_v /= new_I #average velocity
                            new_r /= new_I #average radius
                            new_ϕ /= new_I #average azimuthal angle
                            new_i /= new_I
                            new_η /= new_I #average η
                            subModel.rings[i].I = new_I #set intensity
                            subModel.rings[i].v = new_v #set velocity
                            subModel.rings[i].ϕ = new_ϕ #set azimuthal angle
                            subModel.rings[i].r = new_r #set radius
                            subModel.rings[i].i = new_i #set inclination
                            subModel.rings[i].τ = new_τ #set optical depth
                            subModel.rings[i].η = new_η #set η
                            m.rings[m.subModelStartInds[firstSubModel]+i-1].I = NaN #mark original model index for removal because any leftover point will be combined in last step
                        end
                    end
                end
            else
                #mark original model indices for removal, because they will be added to base model in final step
                for i in 1:length(subModel.rings)
                    subModel.rings[i] = deepcopy(subModel.rings[i]) #no modifications, so need to preserve original model state (will be overwritten in final step)
                    m.rings[m.subModelStartInds[firstSubModel]+i-1].I = NaN 
                end
            end
            if isnothing(baseModel)
                baseModel = subModel #if no base model = only cloud models so set it to subModel
            else
                baseModel += subModel #add free-floating clouds to base model
            end
        end
        for i in 1:length(m.subModelStartInds) #5: add back in non-overlapping points
            camEndInd = i == length(m.subModelStartInds) ? length(m.camera.α) : camStartInds[i+1]-1 #flattened camera index for the next submodel
            camtmp = camera(m.camera.α[camStartInds[i]:camEndInd],m.camera.β[camStartInds[i]:camEndInd],true) #camera α coordinates for submodel i
            endInd = i == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[i+1]-1 #end index of submodel
            baseModel += model(m.rings[m.subModelStartInds[i]:endInd],nothing,camtmp,[1]) #add submodel i to base model
        end
        #5: cleanup redundant points and return
        return removeNaN!(baseModel) 
    end
end
            
            