#!/usr/bin/env julia
using LinearAlgebra

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
    ϕ₀ = atan(yRing,xRing) #original ϕ₀ (no rotation)
    xyzSys = rotate3D(r,ϕ₀,i,rot,θₒPoint) #system coordinates xyz
    ϕ = atan(xyzSys[2],-xyzSys[1]) #ϕ after rotation, measured from +x in disk plane, -x because of how rotation matrix was implemented and desire to have ϕ=0 at +x
    #rMat = getϕ_rMat(i,θₒSystem) #rotation matrix
    #xyzXY = rMat*xyzSys #rotate system plane into XY plane
    #ϕ = atan(xyzXY[2],xyzXY[1]) #ϕ after rotation, measured from +x in disk plane
    return r, ϕ, ϕ₀
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
        ϕ₀: original azimuthal angle in ring plane {Float64}
    """
    cosr = cos(rot); sinr = sin(rot); cosi = cos(i); sini = sin(i); cosθₒ = cos(θₒPoint); sinθₒ = sin(θₒPoint)
    xRing = (β*cosr - α*cosi*sinr)/(cosi*cosθₒ+cosr*sini*sinθₒ) #system x
    yRing = (α*(cosi*cosθₒ+sini/cosr*sinθₒ)+β*cosθₒ*sinr/cosr)/(cosi*cosθₒ/cosr+sini*sinθₒ)
    r = √(xRing^2 + yRing^2)
    ϕ₀ = atan(yRing,xRing) #original ϕ₀ (no rotation)
    xyz[1] = xRing; xyz[2] = yRing; xyz[3] = 0.0
    # xyz[1] = r*cos(ϕ₀); xyz[2] = r*sin(ϕ₀); xyz[3] = 0.0
    # mul!(matBuff,ϕ_rMat,r3D)
    # mul!(colBuff,matBuff,xyz) #rotate system plane into XY plane
    mul!(colBuff,r3D,xyz)
    undo_tilt = [sini 0.0 cosi; 0.0 1.0 0.0; -cosi 0.0 sini]
    mul!(xyz,undo_tilt,colBuff)
    ϕ = atan(xyz[2],-xyz[1]) #ϕ after rotation and being "puffed up", measured from +x in disk plane -- this is fine even for puffed up clouds but note ϕ is measured wrt to disk midplane then. -x because of how rotation matrix was implemented...
    #really this whole thing is stupid and should be removed this is not raytracing
    return r, ϕ, ϕ₀
end

function photograph(r::Float64, ϕ₀::Float64, i::Float64, rot::Float64, θₒ::Float64, reflect::Bool=false)
    """calculate the image coordinates from system coordinates r, ϕ + inclination angle i
    params:
        r: radius from central mass (in terms of rₛ) {Float64}
        ϕ₀: unrotated azimuthal angle in ring plane (rad) {Float64}
        i: inclination angle of ring plane (rad) {Float64}
        rot: rotation of system plane about z axis (rad) {Float64}
        θₒ: ring opening angle {Float64}
        reflect: whether the point is reflected across the midplane of the disk {Bool}
    returns:
        α: image x coordinate (in terms of rₛ) {Float64}
        β: image y coordinate {Float64}
    """
    xyzSys = rotate3D(r,ϕ₀,i,rot,θₒ,reflect)
    α = xyzSys[2] #camera is at +x, so α is y
    β = xyzSys[3] #and β is z
    return α, β
end

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
    #IDisk0/ICloud0 = ratio --> diskCloudIntensityRatio = IDiskFinal/ICloud0 = IDisk0/ICloud0 * (IDiskFinal/IDisk0) = ratio * (IDiskFinal/IDisk0) --> IDIskFinal = IDisk0 * diskCloudIntensityRatio/ratio
    ratio = totalDiskI/totalCloudI
    for ring in m.rings[diskFlagRing]
        ring.I *= diskCloudIntensityRatio/ratio
    end
    reset!(m)
    return m
end

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
    m.camera = camera(newα,newβ,nothing)
    m.rings = m.rings[.!removeFlag]
    reset!(m)
    return m
end

function raytrace!(m::model;IRatios::Union{Float64,Array{Float64,}}=1.0,τCutOff::Float64=1.0,raytraceFreeClouds=false)
    """perform raytracing for a model
    params:
        m: model to raytrace {model}
    returns:
        m: model with raytraced points {model}
    """
    if m.subModelStartInds == [1]
        @warn "raytrace! called on a model with no submodels -- maybe you already raytraced? Returning unaltered model."
        return m
    else
        @info "raytracing model with $(length(m.subModelStartInds)) submodels"
        #1 model with coarsest grid becomes the new "base" model -- if tie pick the one with most points? 
        #2 any models with smaller grid are added to corresponding larger grid model, if they are outside of the larger grid keep original small grid in those places and put them to back of array + set startInds to still show combined
        #3 work from +x camera backwards until τ ~ 1, adding intensity to each pixel as needed
        #4 if τ is fixed assume it applies at all wavelengths, else if τ is a function of wavelength check for spatial and wavelength overlap
        #5 remove extra points from final model if they are obscured 

        #1 detect grid sizes, ranges, types
        # NEED TO REDO THIS AS WELL USING ALPHA AND BETA COORDINATES
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
        contCounter = sum([gridInfo[i][3] for i in 1:length(m.subModelStartInds)]) #number of continuous models
        overlaps = []
        if contCounter > 1
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
        
        baseCam = camera(baseModelα,baseModelβ,nothing) #camera object for base model
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
        #create new model with combined rings and camera points, going from smallest to largest rMin? 
        for i in 1:length(baseModel_rCam)
            println("\n\ni = $i/$(length(baseModel_rCam))\n\n")
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
            #TO DO: these calculations in log space are wrong -- causing uneccessary overlaps. need to calculate a Δr (upper) and Δr (lower) for each pixel (right now upper boundary is correct but lower boundary then goes too far down.)
            rMin = rCam - ΔrDown/2
            rMax = rCam + ΔrUp/2
            ϕMin = ϕCam - Δϕ/2
            ϕMax = ϕCam + Δϕ/2
            i_i = length(base_i) < length(base_I) ? baseModelRings[ring].i : base_i[i] #inclination angle of this point
            θₒi = baseModelRings[ring].θₒ #opening angle of this point
            reflecti = baseModelRings[ring].reflect #reflection flag of this point
            roti = length(baseModelRings[ring].rot) < length(baseModelRings[ring].I) ? baseModelRings[ring].rot : baseModelRings[ring].rot[i] #rotation of this point
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
                    println("\nj = $j/$(length(m.subModelStartInds))")
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
                        println("\nk = $k/$(length(rCamSegment))")
                        rCam = rCamSegment[k]
                        ϕCam = ϕCamSegment[k]
                        #check if this camera point is within the pixel
                        if (rCam >= rMin && rCam <= rMax) && (ϕCam >= ϕMin && ϕCam <= ϕMax)
                            subRing, subCol = getRingFromFlattenedInd(subModel,k)
                            println("submodel = ")
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
                                                
                            m.rings[m.subModelStartInds[j]+subRing-1].I[subCol] = NaN #mark for removal -- all I = 0.0 points removed after pass
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
                    #order = order[start:end] 
                    println("xList = $xList")
                    println("order = $order")
                    println("IList = $IList")
                    println("vList = $vList")
                    println("ϕList = $ϕList")
                    println("rList = $rList")
                    println("iList = $iList")
                    println("τList = $τList")
                    println("ηList = $ηList")
                    new_τ = τList[order][1]
                    new_I = IList[order][1]; new_v = vList[order][1]*new_I; new_ϕ = ϕList[order][1]*new_I; new_r = rList[order][1]*new_I; new_i = iList[order][1]*new_I; new_η = ηList[order][1]*new_I #take first point as the "base" point 
                    j = 2
                    println("new_I (pre loop) = $new_I; new_τ = $new_τ")
                    if all(τ_ΔvList .== Inf) || length(unique(τList)) == 1 #if all τ are the same or τ is a single value, we can just add points until τ ~ 1
                        obscuredFrac = ΔAList[order][j]/maximum(ΔAList[order][1:j-1]) #fraction of area of this tile obscured by those in front of it, i.e. how much of the point is visible -- if < 1 fully covered by previous tile(s)
                        while (new_τ < τCutOff || obscuredFrac > 1.0) && j <= length(order) #add points until τ ~ 1 --TO DO: add obscuration component to this, i.e. if what if thing in front has τ ~ 1 but doens't fully obscure ring behind? 
                            if !isnan(IList[order][j]) && !isnan(xList[order][j]) #skip points that are marked for removal or unphysical points (r = NaN)
                                tmp_I = exp(-new_τ)*IList[order][j] #assume new intensity is added to previous one in chunks, no absorption 
                                if obscuredFrac > 1.0
                                    tmp_I = exp(-new_τ)*IList[order][j]*obscuredFrac #if obscuredFrac > 1, then this point is not fully obscured by previous points, so obscured intensity needs to be reduced
                                    tmp_I += IList[order][j]*(obscuredFrac-1) #add unboscured intensity from previous point
                                else
                                end
                                new_v += vList[order][j]*tmp_I 
                                new_r += rList[order][j]*tmp_I 
                                new_ϕ += ϕList[order][j]*tmp_I
                                new_i += iList[order][j]*tmp_I 
                                new_η += ηList[order][j]*tmp_I
                                new_I += tmp_I #add intensity
                                new_τ = τList[order][j] #new optical depth after adding this point 
                            end
                            j += 1
                        end
                        den = new_I == 0.0 ? 1.0 : new_I #avoid division by zero
                        new_v /= den #average velocity
                        new_r /= den #average radius
                        new_ϕ /= den #average azimuthal angle
                        new_i /= den #average inclination
                        new_η /= den #average η
                    else #think about how to do this properly...
                        error("velocity-dependent optical depths not yet implemented -- pass τ as a float when creating models if you want to use raytracing")
                    end      
                    println("new_I (post loop) = $new_I; new_τ = $new_τ")
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
        #call removal function on all submodels to remove points with I = 0.0
        #m = removeZeroedPoints!(m) #remove points with I = 0.0 from starting model
        #call raytrace again with new "model" that is sum of all models with overlaps but not overlapping base model
        if length(overlaps) > 1 
            for (i,j) in overlaps #both continuous models 
                #check if i or j is the base model, if so skip it
                if i != baseModelInd && j != baseModelInd
                #i is new "base" model, check for all other models that overlap with it 
                    endInd = i = length(m.subModelStartInds) ? length(m.camera.α) : camStartInds[i+1]-1 #flattened camera index for the next submodel
                    camera_i = camera(m.camera.α[camStartInds[i]:endInd],m.camera.β[camStartInds[i]:endInd],nothing) #camera α coordinates for submodel i
                    endInd = j == length(m.subModelStartInds) ? length(m.camera.α) : camStartInds[j+1]-1 #flattened camera index for the next submodel
                    camera_j = camera(m.camera.α[camStartInd])
                    endInd = i == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[i+1]-1 #ring index for the next submodel
                    subModelI = model(m.rings[m.subModelStartInds[i]:endInd],nothing,camera_i,[1]) #submodel i to check
                    endInd = j == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[j+1]-1 #ring index for the next submodel
                    subModelI += model(m.rings[m.subModelStartInds[j]:endInd],nothing,camera_j,[1]) #add j to check
                    IRatiosTmp = deepcopy(IRatios[i,j])
                    for ii in 1:length(m.subModelStartInds) #check all submodels for overlaps with submodel i
                        if ii != i && ii != j #not same pair
                            if ii != baseModelInd #not base model
                                endInd = ii = length(m.subModelStartInds) ? length(m.camera.α) : camStartInds[ii+1]-1 #flattened camera index for the next submodel
                                camera_ii = camera(m.camera.α[camStartInds[ii]:endInd],m.camera.β[camStartInds[ii]:endInd],nothing) #camera α coordinates for submodel ii
                                #camera_ii = camera(m.camera.α[m.[ii]:m.subModelStartInds[ii+1]-1],m.camera.β[m.subModelStartInds[ii]:m.subModelStartInds[ii+1]-1],nothing) #camera α coordinates for submodel ii
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
                                #subModelI += model(m.rings[m.subModelStartInds[jj]:endInd],nothing,camera_ii,[1]) #add jj to check
                                push!(IRatios,IRatiosTmp[jj]) #add multiplier for jj
                            end
                        end
                    end
                    subModelI = raytrace!(subModelI;IRatios=IRatiosTmp,τCutOff=τCutOff) #raytrace submodel i
                    #remove NaN points from submodel i
                    subModelI = removeNaN!(subModelI) 
                    #m = removeZeroedPoints!(m) #remove overlapping clouds that were added to base model at this stage
                    #call removal function to remove points with I = 0.0
                    #add submodel i to base model
                    baseModel += subModelI #baseModel now has every continuous model + overlapping clouds, just need to add back in non-overlapping clouds
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
            camtmp = camera(m.camera.α[camStartInds[firstSubModel]:endInd],m.camera.β[camStartInds[firstSubModel]:endInd],nothing) #camera α coordinates for first discrete model
            endInd = firstSubModel == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[firstSubModel+1]-1 #ring index for the next submodel
            subModel = model(m.rings[m.subModelStartInds[firstSubModel]:endInd],nothing,camtmp,[1]) #submodel to check
            for i in firstSubModel:length(m.subModelStartInds) 
                if gridInfo[i][3] == 1 
                    endInd = i == length(m.subModelStartInds) ? length(m.camera.α) : camStartInds[i+1]-1 #flattened camera index for the next submodel
                    camtmp = camera(m.camera.α[camStartInds[i]:endInd],m.camera.β[camStartInds[i]:endInd],nothing) #camera α coordinates for submodel i
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
                            #order = order[start:end]
                            new_τ = τList[order][1]
                            new_I = IList[order][1]; new_v = vList[order][1]*new_I; new_ϕ = ϕList[order][1]*new_I; new_r = rList[order][1]*new_I; new_i = iList[order][1]*new_I; new_η = ηList[order][1]*new_I #take first point as the "base" point 
                            j = 2
                            while new_τ < τCutOff && j <= length(order) #add points until τ ~ 1
                                if !isnan(IList[order][j]) && !isnan(xList[order][j]) #skip points that are marked for removal
                                    tmp_I = exp(-new_τ)*IList[order][j] #assume new intensity is added to previous one in chunks, no absorption 
                                    new_v += vList[order][j]*tmp_I 
                                    new_r += rList[order][j]*tmp_I 
                                    new_ϕ += ϕList[order][j]*tmp_I
                                    new_i += iList[order][j]*tmp_I 
                                    new_η += ηList[order][j]*tmp_I
                                    new_I += tmp_I #add intensity
                                    new_τ = τlist[order][j] #new optical depth after adding this point 
                                end
                                j += 1
                            end
                            new_v /= new_I #average velocity
                            new_r /= new_I #average radius
                            new_ϕ /= new_I #average azimuthal angle
                            new_i /= new_I
                            new_η /= new_I #average η
                            subModel.rings[i].I[1] = new_I #set intensity
                            subModel.rings[i].v[1] = new_v #set velocity
                            subModel.rings[i].ϕ[1] = new_ϕ #set azimuthal angle
                            subModel.rings[i].r[1] = new_r #set radius
                            if typeof(subModel.rings[ring].i) == Float64 #if i is a single value, fill with vector then replace point with new_i
                                for r in 1:length(subModel.rings)
                                    subModel.rings[r].i = fill(subModel.rings[r].i,length(subModel.rings[r].I)) #fill with vector of same length as I
                                end
                                sub_i = getVariable(subModel,:i,flatten=true) #get new sub_i after filling
                            end
                            subModel.rings[ring].i[col] = new_i #set inclination
                            if typeof(subModel.rings[ring].τ) == Float64 #if τ is a single value, fill with vector then replace point with new_τ
                                for r in 1:length(subModel.rings)
                                    subModel.rings[r].τ = fill(subModel.rings[r].τ,length(subModel.rings[r].I)) #fill with vector of same length as I
                                end
                                sub_τ = getVariable(subModel,:τ,flatten=true) #get new sub_τ after filling
                            end
                            subModel.rings[ring].τ[col] = new_τ #set optical depth
                            if typeof(subModel.rings[ring].η) == Float64 #if η is a single value, fill with vector then replace point with new_η
                                for r in 1:length(subModel.rings)
                                    subModel.rings[r].η = fill(subModel.rings[r].η,length(subModel.rings[r].I)) #fill with vector of same length as I
                                end
                                sub_η = getVariable(subModel,:η,flatten=true) #get new sub_η after filling
                            end
                            subModel.rings[ring].η[col] = new_η #set η
                            m.rings[m.subModelStartInds[firstSubModel]+i-1].I = NaN #mark original model index for removal because any leftover point will be combined in last step
                        end
                    end
                end
            end
            baseModel += subModel #add free-floating clouds to base model
            baseModel.subModelStartInds = [1] #reset submodel start indices to 1, since we only have one submodel now
        end
        for i in 1:length(m.subModelStartInds)
            camEndInd = i == length(m.subModelStartInds) ? length(m.camera.α) : camStartInds[i+1]-1 #flattened camera index for the next submodel
            camtmp = camera(m.camera.α[camStartInds[i]:camEndInd],m.camera.β[camStartInds[i]:camEndInd],nothing) #camera α coordinates for submodel i
            endInd = i == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[i+1]-1 #end index of submodel
            baseModel += model(m.rings[m.subModelStartInds[i]:endInd],nothing,camtmp,[1]) #add submodel i to base model
        end
        #baseModel += m #add original model to base model so that any points that didn't need to be raytraced are still included                 
        #cleanup redundant points and return
        return baseModel #debug
        return removeNaN!(baseModel) 
    end
end


# function raytrace!(m::model;DiskRayFrac::Union{Float64,Nothing}=nothing,nRays::Union{Float64,Int,Nothing}=nothing,τMax::Float64=1.0,trackRays::Bool=false,rotate3D::Function=rotate3D)
#     #new approach -- brute force raytracing
#     isCombined, startInds, diskFlag = detectCombinedModel(m)
#     if !isCombined
#         @warn "did not detect combined model so nothing to raytrace -- returning unaltered input model"
#         return m
#     end
#     if length(startInds) > 2
#         #raytrace recursively combining two at a time
#         mList = [deepcopy(m) for i=1:length(startInds)]
#         for (i,m) in enumerate(mList)
#             s = startInds[i]; e = i == length(startInds) ? length(m.rings) : startInds[i+1]-1
#             m.rings = m.rings[s:e]
#         end
#         mFinal = mList[1]
#         for i=2:length(mList)
#             mFinal = mFinal + mList[i]
#             mFinal = raytrace!(mFinal,DiskRayFrac,nRays,camera,τMax)
#         end
#         return mFinal
#     end

#     ΔA = vcat([getVariable(m,:ΔA)...])
#     diskFlag = Array{Bool}(undef,length(ΔA))
#     ind = 1
#     for ring in m.rings
#         if typeof(ring.ϕ) == Float64 && typeof(ring.r) == Float64
#             diskFlag[ind] = false
#             ind += 1
#         else
#             stride = length(ring.ϕ)
#             diskFlag[ind:ind+stride] .= true
#             ind += stride
#         end
#     end
#     if sum(diskFlag) == 0
#         @warn "no disk cells detected -- simply returning combined model with no raytracing"
#         return m
#     end

#     if isnothing(DiskRayFrac)
#         DiskRayFrac = sum(ΔA[diskFlag])/sum(ΔA[.!diskFlag]) #ratio of disk to cloud area on camera
#     end
#     @debug "DiskRayFrac: $DiskRayFrac"
#     if typeof(nRays) == Float64
#         nRays = ceil(Int,nRays)
#     end
#     nClouds = sum(.!diskFlag) #number of clouds
#     if isnothing(nRays)
#         nRays = ceil(Int,nClouds*(1+DiskRayFrac)) #total number of rays to trace
#     end
#     @debug "nRays: $nRays"
#     sampleFrac = 1.0
#     if nRays < nClouds*(1+DiskRayFrac)
#         sampleFrac = nRays/(nClouds*(1+DiskRayFrac))
#     end
#     @debug "Sampling fraction: $sampleFrac"
#     if sampleFrac < 0.1
#         @warn "Number of rays is less than 10% of recommended value ($(nClouds*(1+DiskRayFrac)))-- results may be unreliable."
#     end
#     αₒ = vcat([m.camera.α...]); βₒ = vcat([m.camera.β...])
#     diskα = αₒ[diskFlag]; diskβ = βₒ[diskFlag] 
#     cloudα = αₒ[.!diskFlag]; cloudβ = βₒ[.!diskFlag]
#     diskr = sqrt.(diskα.^2 .+ diskβ.^2)
#     cloudr = sqrt.(cloudα.^2 .+ cloudβ.^2)
#     minDiskr = minimum(diskr); maxDiskr = maximum(diskr)
#     cloudMaskAlone = [(r<minDiskr || r>maxDiskr) for r in cloudr] #mask for clouds that are not in disk area
#     # if isnothing(camera)
#     #     camera = m.camera
#     #     α = vcat(diskα,cloudα[cloudMask]); β = vcat(diskβ,cloudβ[cloudMask])
#     #     camera = camera(α,β)
#     #     m.camera = camera
#     # else -- should allow user to specify new camera?
#     #     m.camera = camera
#     # end
#     #not allowing new camera specification, instead just recalculating from existing camera and removing cloud cells that are inside of disk cells
#     α = vcat(diskα,cloudα[cloudMaskAlone]); β = vcat(diskβ,cloudβ[cloudMaskAlone])
#     cam = camera(α,β,nothing)
#     m.camera = cam
#     #add up intensities (with raytracing) and divide by total number of rays in each cell -- then the intensity in each cell is ∝ 1/ΔA and can weight by ΔA in binModel
#     #new camera will be a "disk" camera, so every cell at the end will have ΔA = ΔADisk. This is fine even if a cloud is alone because if there are enough rays most of them will miss in said cell and the intensity will be low
#     #also avoids problem of neighboring cells etc. just calculate what we hit in each cell and divide by total number of rays -- if a cloud is across 2 cells it will get hit with equal probability in both of those cells! 
#     #start with cloud ray paths, then randomly sample new rays from range of r,ϕ allowed by camera. Note that this may mess up statistics if there are very few rays? Need to also have an "empty area" fraction in calculating nRays!
#     #or....set ΔA in cloud only cells to sum(ΔAcloud)/sum(ΔA)/nClouds(?) (average cloud contribution to total area) -- do we need to fill in rest of "ring" at each cloud radius with corresponding number of empty ϕ cells? use sparse arrays?

#     rCam = sqrt.(cam.α.^2 .+ cam.β.^2); ϕCam = atan.(cam.β,cam.α)
#     rUnique = unique(r -> round(r,sigdigits=12),rCam) #note that at 8 sigfigs some clouds were actually identical and this was causing an error??
#     ϕUnique = unique(ϕ -> round(ϕ,sigdigits=12),ϕCam) #note really want a list of ϕ disk values? don't care about cloud ϕ because we use this to index into rings with multiple ϕ values
#     rMin = minimum(rCam); rMax = maximum(rCam)
#     ϕMin = minimum(ϕCam); ϕMax = maximum(ϕCam)

#     #initialize new rings
#     diskFlagRing = [!(typeof(r.ϕ) == Float64 && typeof(r.r) == Float64) for r in m.rings] #true if disk ring, false if cloud ring
#     ϕDisk = [m.rings[diskFlagRing][i].ϕ for i=1:sum(diskFlagRing)]
#     newRingsLength = length(rUnique)
#     newRings = Array{ring}(undef,newRingsLength)
#     for i=1:sum(diskFlagRing) #do disk first to match rCam, ϕCam
#         newRings[i] = deepcopy(m.rings[diskFlagRing][i])
#     end
#     for i=sum(diskFlagRing)+1:newRingsLength #do clouds next
#         newRings[i] = deepcopy(m.rings[.!diskFlagRing][cloudMaskAlone][i-sum(diskFlagRing)])
#         if rand() > sampleFrac
#             newRings[i].I = 0.0 #zero out intensity of fraction of cloud cells that won't be sampled
#         end
#         #don't need to zero out inensity in general here because no actual raytracing will happen in this cell (no disk to pass through, I just stays I_cloud)
#     end
#     old_I = [deepcopy(r.I) for r in newRings[1:sum(diskFlagRing)]] #disk intensities before raytracing
#     rayCounts = [deepcopy(r.I) for r in newRings] #counter for number of rays in each disk cell, copy I just for shape
#     rays = trackRays ? Array{ray}(undef,nRays) : nothing
#     for i=1:sum(diskFlagRing)
#         newRings[i].I .= 0.0 #zero out intensity before raytracing
#         rayCounts[i] .= 0.0 #zero out ray counters
#     end

#     #clouds first
#     cloudInds = findall(.!diskFlag)
#     cloudRingInds = findall(.!diskFlagRing)
#     raysTraced = 0
#     for i in cloudInds #raytrace cloud points first
#         if rand() < sampleFrac
#             @debug "in cloud ray loop -- total progress = $(round(raysTraced/nRays*100,sigdigits=2))%\r"
#             α = αₒ[i]; β = βₒ[i]
#             r = sqrt(α^2 + β^2); ϕ = atan(β,α)
#             cloud_dists = sqrt.((cloudα.-α).^2 .+ (cloudβ.-β).^2) #how far is each cloud center from this ray?
#             cloud_radii = sqrt.(ΔA[cloudInds]./π) #projected radius of each cloud, assuming spherical clouds
#             cloudMaskRay = cloud_dists .< cloud_radii #mask for clouds that are hit by this ray
#             newRingsInd = argmin(abs.(r.-rUnique)) #index of closest r value in camera, which should also match newRings index
#             if r > minDiskr && r < maxDiskr #cloud is in disk area
#                 I = 0.0; τ = 0.0; x = 0.0
#                 if trackRays
#                     I = zeros(1+sum(cloudMaskRay)) #initialize intensity array for ray struct
#                     τ = zeros(1+sum(cloudMaskRay)) #initialize optical depth array for ray struct
#                     x = zeros(1+sum(cloudMaskRay)) #initialize system x array for ray struct
#                 end

#                 #raytrace through disk area
#                 #need to find closest disk cell to cloud cell and raytrace through that cell
#                 newRingsϕ = argmin(abs.(ϕ.-ϕDisk[newRingsInd])) #closest ϕ value in camera, should match order of ϕ in newRings[r].ϕ
#                 disk_I = old_I[newRingsInd][newRingsϕ] #intensity of disk cell without raytracing
#                 τDisk = typeof(newRings[newRingsInd].τ) == Float64 ? newRings[newRingsInd].τ : newRings[newRingsInd].τ[newRingsϕ]
#                 xDisk = rotate3D(newRings[newRingsInd].r[newRingsϕ],newRings[newRingsInd].ϕ₀[newRingsϕ],newRings[newRingsInd].i,newRings[newRingsInd].rot,newRings[newRingsInd].θₒ)[1] #system coordinates xyz
#                 #xCloud = rotate3D(m.rings[i].r,m.rings[i].ϕ₀,m.rings[i].i,m.rings[i].rot,m.rings[i].θₒ,m.rings[i].reflect)[1] #system coordinates xyz
#                 #xClouds = rotate3D.(m.rings[cloudInds][cloudMask].r,m.rings[cloudInds][cloudMask].ϕ₀,m.rings[cloudInds][cloudMask].i,m.rings[cloudInds][cloudMask].rot,m.rings[cloudInds][cloudMask].θₒ,m.rings[cloudInds][cloudMask].reflect) #system coordinates xyz
#                 xClouds = [rotate3D(r.ϕ,r.ϕ₀,r.i,r.rot,r.θₒ,r.reflect)[1] for r in m.rings[cloudRingInds][cloudMaskRay]]
                
#                 IList = [disk_I,[r.I for r in m.rings[cloudRingInds][cloudMaskRay]]...]
#                 τList = [τDisk,[r.τ for r in m.rings[cloudRingInds][cloudMaskRay]]...]
#                 xList = [xDisk,xClouds...]
                
#                 xOrder = sortperm(xList,rev=true)
#                 ILocal = IList[xOrder][1]; τLocal = τList[xOrder][1]; xLocal = xList[xOrder][1]
#                 if trackRays
#                     I[1] = ILocal; τ[1] = τLocal; x[1] = xLocal
#                 end
#                 counter = 1; stop = (τLocal > τMax) || (counter+1 > length(xList))
#                 while !stop
#                     counter += 1
#                     xCloud = xList[xOrder][counter]
#                     ILocal += IList[xOrder][counter]*exp(-τLocal) #add intensity 
#                     τLocal += τList[xOrder][counter] #update optical depth after passing through 
#                     if trackRays
#                         I[counter] = deepcopy(ILocal); τ[counter] = deepcopy(τLocal); x[counter] = deepcopy(xCloud)
#                     end
#                     stop = (τLocal > τMax) || (counter+1 > length(xList))
#                 end
#                 newRings[newRingsInd].I[newRingsϕ] += ILocal #add intensity to disk cell
#                 rayCounts[newRingsInd][newRingsϕ] += 1 #increment ray counter
#                 raysTraced += 1
#                 if trackRays
#                     rays[raysTraced] = ray(r,ϕ,α,β,τ[1:counter],x[1:counter],I[1:counter],1)
#                 end
#             else 
#                 if sum(cloudMaskRay) > 1 #ray through cloud overlaps with other clouds
#                     I = 0.0; τ = 0.0; x = 0.0
#                     if trackRays
#                         I = zeros(sum(cloudMaskRay)) #initialize intensity array for ray struct
#                         τ = zeros(sum(cloudMaskRay)) #initialize optical depth array for ray struct
#                         x = zeros(sum(cloudMaskRay)) #initialize system x array for ray struct
#                     end
#                     xClouds = [rotate3D(r.ϕ,r.ϕ₀,r.i,r.rot,r.θₒ,r.reflect)[1] for r in m.rings[cloudInds][cloudMaskRay]]
#                     IList = [r.I for r in m.rings[cloudInds][cloudMaskRay]]
#                     τList = [r.τ for r in m.rings[cloudInds][cloudMaskRay]]
#                     xOrder = sortperm(xClouds,rev=true)
#                     ILocal = IList[xOrder][1]; τLocal = τList[xOrder][1]; xLocal = xClouds[xOrder][1]
#                     if trackRays
#                         I[1] = ILocal; τ[1] = τLocal; x[1] = xLocal
#                     end
#                     counter = 1; stop = (τLocal > τMax) || (counter+1 > length(xClouds))
#                     while !stop
#                         counter += 1
#                         xLocal = xClouds[xOrder][counter]
#                         ILocal += IList[xOrder]*exp(-τLocal) #add intensity 
#                         τLocal += τList[xOrder][counter] #update optical depth after passing through 
#                         if trackRays
#                             I[counter] = deepcopy(ILocal); τ[counter] = deepcopy(τLocal); x[counter] = deepcopy(xLocal)
#                         end
#                         stop = (τLocal > τMax) || (counter+1 > length(xClouds))
#                     end
#                     newRings[newRingsInd].I += ILocal #add intensity to disk cell
#                     rayCounts[newRingsInd] += 1 #increment ray counter
#                     raysTraced += 1
#                     if trackRays
#                         rays[raysTraced] = ray(r,ϕ,α,β,τ[1:counter],x[1:counter],I[1:counter],2)
#                     end
#                 else #cloud alone
#                     # newRings[newRingsInd].I += m.rings[i].I #copy cloud cell intensity to new rings struct
#                     # xCloud = rotate3D(m.rings[i].r,m.rings[i].ϕ₀,m.rings[i].i,m.rings[i].rot,m.rings[i].θₒ,m.rings[i].reflect)[1] #system coordinates xyz
#                     # rayCounts[newRingsInd] += 1 #increment ray counter
#                     # raysTraced += 1
#                     # if trackRays
#                     #     rays[raysTraced] = ray(r,ϕ,α,β,[[m.rings[i].τ]],[xCloud],[m.rings[i].I])
#                     # end

#                     if length(m.rings[cloudRingInds][cloudMaskRay]) == 1
#                         if rayCounts[newRingsInd] == 0
#                             println("newRings[newRingsInd] = $(newRings[newRingsInd])")
#                             println("m.rings[cloudRingInds][cloudMaskRay][1] = $(m.rings[cloudRingInds][cloudMaskRay][1])")
#                             newRings[newRingsInd] = deepcopy(m.rings[cloudRingInds][cloudMaskRay][1]) #copy cloud cell to new rings struct
#                         else
#                             newRings[newRingsInd].I += m.rings[cloudRingInds][cloudMaskRay][1].I #copy cloud cell intensity to new rings struct
#                             newRings[newRingsInd].I += m.rings[cloudRingInds][cloudMaskRay][1].I #copy cloud cell intensity to new rings struct
#                         end
#                         rayCounts[newRingsInd] += 1 #increment ray counter
#                         xCloud = rotate3D(m.rings[cloudRingInds][cloudMaskRay][1].r,m.rings[cloudRingInds][cloudMaskRay][1].ϕ₀,m.rings[cloudRingInds][cloudMaskRay][1].i,m.rings[cloudRingInds][cloudMaskRay][1].rot,m.rings[cloudRingInds][cloudMaskRay][1].θₒ,m.rings[cloudRingInds][cloudMaskRay][1].reflect)[1] #system coordinates xyz
#                         raysTraced += 1
#                         if trackRays
#                             rays[raysTraced] = ray(r,ϕ,α,β,[m.rings[cloudRingInds][cloudMaskRay][1].τ],[xCloud],[m.rings[cloudRingInds][cloudMaskRay][1].I],3)
#                         end
# \                    else
#                         println("CLOUD ALONE: i = $i")
#                         println("typeof(newRings[newRingsInd].I) = $(typeof(newRings[newRingsInd].I))")
#                         println("typeof(m.rings[cloudRingInds][cloudMaskRay][1].I) = $(typeof(m.rings[cloudRingInds][cloudMaskRay][1].I))")
#                         println("r = $r -- (minDiskr,maxDiskr) = ($minDiskr,$maxDiskr)")
#                         println("newRings[newRingsInd] = $(newRings[newRingsInd])")
#                         println("m.rings[cloudRingInds][cloudMaskRay] = $(m.rings[cloudRingInds][cloudMaskRay])")
#                         error("cloud alone mask returned more than one ring struct -- should not happen")
#                     end
#                 end
#             end
#         end
#     end
#     #rest of system randomly draw rays
#     minα = minimum(αₒ); maxα = maximum(αₒ)
#     minβ = minimum(βₒ); maxβ = maximum(βₒ)
#     αD = Uniform(minα,maxα); βD = Uniform(minβ,maxβ)
#     for rayi=raysTraced+1:nRays #fill in rest of rays with random rays from camera
#         @debug "in random ray loop -- total progress = $(round(rayi/nRays*100,sigdigits=2))%\r"
#         α = rand(αD)
#         β = rand(βD)
#         r = sqrt(α^2 + β^2); ϕ = atan(β,α)
#         cloud_dists = sqrt.((cloudα.-α).^2 .+ (cloudβ.-β).^2) #how far is each cloud center from this ray?
#         cloud_radii = sqrt.(ΔA[cloudInds]./π) #projected radius of each cloud, assuming spherical clouds
#         cloudMaskRay = cloud_dists .< cloud_radii #mask for clouds that are hit by this ray
#         allowed_draw = (r > minDiskr && r < maxDiskr) || sum(cloudMaskRay) > 0 #no rays through empty space
#         while !allowed_draw
#             α = rand(Uniform(minimum(αₒ),maximum(αₒ)))
#             β = rand(Uniform(minimum(βₒ),maximum(βₒ)))
#             r = sqrt(α^2 + β^2); ϕ = atan(β,α)
#             cloud_dists = sqrt.((cloudα.-α).^2 .+ (cloudβ.-β).^2) #
#             cloudMaskRay = cloud_dists .< cloud_radii 
#             allowed_draw = (r > minimum(diskr) && r < maximum(diskr)) || sum(cloudMaskRay) > 0 #no rays through empty space
#         end
#         newRingsInd = argmin(abs.(r.-rUnique)) #index of closest r value in camera, which should also match newRings index
#         if r > minDiskr && r < maxDiskr #random ray in disk area
#             newRingsϕ = argmin(abs.(ϕ.-ϕDisk[newRingsInd])) #closest ϕ value in camera, should match order of ϕ in newRings[r].ϕ
#             if sum(cloudMaskRay) > 0
#                 I = 0.0; τ = 0.0; x = 0.0
#                 if trackRays
#                     I = zeros(1+sum(cloudMaskRay)) #initialize intensity array for ray struct
#                     τ = zeros(1+sum(cloudMaskRay)) #initialize optical depth array for ray struct
#                     x = zeros(1+sum(cloudMaskRay)) #initialize system x array for ray struct
#                 end
#                 disk_I = old_I[newRingsInd][newRingsϕ] #intensity of disk cell without raytracing
#                 xDisk = rotate3D(newRings[newRingsInd].r[newRingsϕ],newRings[newRingsInd].ϕ₀[newRingsϕ],newRings[newRingsInd].i,newRings[newRingsInd].rot,newRings[newRingsInd].θₒ)[1] #system coordinates xyz
#                 xClouds = [rotate3D(r.ϕ,r.ϕ₀,r.i,r.rot,r.θₒ,r.reflect)[1] for r in m.rings[cloudRingInds][cloudMaskRay]]
                
#                 IList = [disk_I,[r.I for r in m.rings[cloudRingInds][cloudMaskRay]]...]
#                 τDisk = typeof(newRings[newRingsInd].τ) == Float64 ? newRings[newRingsInd].τ : newRings[newRingsInd].τ[newRingsϕ]
#                 τList = [τDisk,[r.τ for r in m.rings[cloudRingInds][cloudMaskRay]]...]
#                 xList = [xDisk,xClouds...]
                
#                 xOrder = sortperm(xList,rev=true)
#                 ILocal = IList[xOrder][1]; τLocal = τList[xOrder][1]; xLocal = xList[xOrder][1]
#                 if trackRays
#                     I[1] = ILocal; τ[1] = τLocal; x[1] = xLocal
#                 end
#                 counter = 1; stop = (τLocal > τMax) || (counter+1 > length(xList))
#                 while !stop
#                     counter += 1
#                     xLocal = xList[xOrder][counter]
#                     ILocal += IList[xOrder][counter]*exp(-τLocal) #add intensity 
#                     τLocal += τList[xOrder][counter] #update optical depth after passing through 
#                     if trackRays
#                         I[counter] = deepcopy(ILocal); τ[counter] = deepcopy(τLocal); x[counter] = deepcopy(xLocal)
#                     end
#                     stop = (τLocal > τMax) || (counter+1 > length(xList))
#                 end
#                 newRings[newRingsInd].I[newRingsϕ] += ILocal #add intensity to disk cell
#                 rayCounts[newRingsInd][newRingsϕ] += 1 #increment ray counter
#                 raysTraced += 1
#                 if trackRays
#                     rays[raysTraced] = ray(r,ϕ,α,β,τ[1:counter],x[1:counter],I[1:counter],4)
#                 end
#             else
#                 newRings[newRingsInd].I[newRingsϕ] += old_I[newRingsInd][newRingsϕ] #add disk cell intensity to new rings struct, no raytracing needed
#                 rayCounts[newRingsInd][newRingsϕ] += 1 #increment ray counter
#                 xDisk = rotate3D(newRings[newRingsInd].r[newRingsϕ],newRings[newRingsInd].ϕ₀[newRingsϕ],newRings[newRingsInd].i,newRings[newRingsInd].rot,newRings[newRingsInd].θₒ)[1] #system coordinates xyz
#                 τDisk = typeof(newRings[newRingsInd].τ) == Float64 ? newRings[newRingsInd].τ : newRings[newRingsInd].τ[newRingsϕ]
#                 raysTraced += 1
#                 if trackRays
#                     rays[raysTraced] = ray(r,ϕ,α,β,[τDisk],[xDisk],[newRings[newRingsInd].I[newRingsϕ]],5)
#                 end
#             end
#         else #ray through at least one cloud, not in disk area
#             if sum(cloudMaskRay) > 1 #ray through cloud overlaps with other clouds
#                 I = 0.0; τ = 0.0; x = 0.0
#                 if trackRays
#                     I = zeros(sum(cloudMaskRay)) #initialize intensity array for ray struct
#                     τ = zeros(sum(cloudMaskRay)) #initialize optical depth array for ray struct
#                     x = zeros(sum(cloudMaskRay)) #initialize system x array for ray struct
#                 end
#                 xClouds = [rotate3D(r.ϕ,r.ϕ₀,r.i,r.rot,r.θₒ,r.reflect)[1] for r in m.rings[cloudRingInds][cloudMaskRay]]
#                 IList = [r.I for r in m.rings[cloudRingInds][cloudMaskRay]]
#                 τList = [r.τ for r in m.rings[cloudRingInds][cloudMaskRay]]
#                 xOrder = sortperm(xClouds,rev=true)
#                 ILocal = IList[xOrder][1]; τLocal = τList[xOrder][1]; xLocal = xClouds[xOrder][1]
#                 if trackRays
#                     I[1] = ILocal; τ[1] = τLocal; x[1] = xLocal
#                 end
#                 counter = 1; stop = (τLocal > τMax) || (counter+1 > length(xClouds))
#                 while !stop
#                     counter += 1
#                     xLocal = xList[xOrder][counter]
#                     ILocal += IList[xOrder]*exp(-τLocal) #add intensity 
#                     τLocal += τList[xOrder][counter] #update optical depth after passing through 
#                     if trackRays
#                         I[counter] = deepcopy(ILocal); τ[counter] = deepcopy(τLocal); x[counter] = deepcopy(xLocal)
#                     end
#                     stop = (τLocal > τMax) || (counter+1 > length(xList))
#                 end
#                 newRings[newRingsInd].I += ILocal #add intensity to disk cell
#                 rayCounts[newRingsInd] += 1 #increment ray counter
#                 raysTraced += 1
#                 if trackRays
#                     rays[raysTraced] = ray(r,ϕ,α,β,τ[1:counter],x[1:counter],I[1:counter],6)
#                 end
#             else #cloud alone
#                 if length(m.rings[cloudRingInds][cloudMaskRay]) == 1
#                     newRings[newRingsInd].I += m.rings[cloudRingInds][cloudMaskRay][1].I #copy cloud cell intensity to new rings struct
#                     newRings[newRingsInd].I += m.rings[cloudRingInds][cloudMaskRay][1].I #copy cloud cell intensity to new rings struct
#                     rayCounts[newRingsInd] += 1 #increment ray counter
#                     xCloud = rotate3D(m.rings[cloudRingInds][cloudMaskRay][1].r,m.rings[cloudRingInds][cloudMaskRay][1].ϕ₀,m.rings[cloudRingInds][cloudMaskRay][1].i,m.rings[cloudRingInds][cloudMaskRay][1].rot,m.rings[cloudRingInds][cloudMaskRay][1].θₒ,m.rings[cloudRingInds][cloudMaskRay][1].reflect)[1] #system coordinates xyz
#                     raysTraced += 1
#                     if trackRays
#                         rays[raysTraced] = ray(r,ϕ,α,β,[m.rings[cloudRingInds][cloudMaskRay][1].τ],[xCloud],[m.rings[cloudRingInds][cloudMaskRay][1].I],7)
#                     end
#                 else
#                     error("cloud alone mask returned more than one ring struct -- should not happen")
#                 end
#             end
#         end
#     end
#     for i=1:newRingsLength
#         if typeof(newRings[i].I) == Float64
#             newRings[i].I /= rayCounts[i] #cloud, divide by number of rays in each cell
#         else
#             newRings[i].I ./= rayCounts[i] #broadcast, arrays should be same size from setup
#         end
#     end
#     m.rings = newRings
#     reset!(m)
#     if trackRays
#         m.camera.rays = rays
#     end
#     return m
# end



# function raytrace!(m::model;camera::Union{camera,Nothing}=nothing,τMax=1.0)
#     #don't know if this actually works, but general idea is right

#     #NOTES 8/27: 
#     #should initialize every "ray" coming from an associated α and β camera coordinate from system
#     #go backwards, and if you are in front of something else, do calculation
#     #if unobscured, just return I at that point. 

#     #no no -- do it like old disk way, but modify if there are clouds on top!
#     #for disk we have ΔA -- projected area on camera. just need to sum up intensities inside of this area! then integrand still has ΔA when generating profiles
#     #need to make new rings struct where each ring will have ΔA = ΔAdisk if disk is visible, otherwise ΔAcloud (if just cloud and no disk)
#     #set normalization parameter between total cloud/disk flux beforehand, then just do radiative transfer on top of that 
#     #algorithm: 
#     #1: check each ring for whether it is a cloud or disk element (check if ϕ and r are floats or arrays)
#     #2: go through disk rings first -- at each disk ring define ΔA = ΔAdisk and check if there are multiple instances 
#     #  --a: check if any clouds within ΔA
#     #       do this manually by checking if cloud r,ϕ inside of ΔA. find cloud r,ϕ from α,β and compare to bounds of ΔA (calculate r and ϕ range below)
#     #       recreate r,ϕ (camera) with sqrt.(α.^2 .+ β.^2) and atan.(β,α)
#     #       FOR LOG COORDS: nextr = r + (r[i]+r[i+1])/2*Δr; nextϕ = ϕ + Δϕ
#     #                       diff ([r[i+1]-r[i] for i...]) = (r[i]+r[i+1])/2*Δr --> Δr = 2*(r[i+1]-r[i])/(r[i]+r[i+1])
#     #       FOR EITHER: Δϕ = ϕ[i+1]-ϕ[i]
#     #       then r range = r, nextr, and ϕ range = ϕ, ϕ+Δϕ (but make sure to normalize this to be -π to π)
#     #      --i: if no, I = Idisk*ΔADisk
#     #      --ii: if yes, raytrace backwards through each cloud path. dI = I*exp(-τ)*ΔAcloud where τ is the optical depth along path.
#     #            τ stored in each ring struct is the τ after passing through that point, so sum up τ until τMax is reached
#     #            disk contribution along path is ΔAcloud * Idisk * exp(-τ)
#     #            if ΔAcloud smaller than ΔAdisk (which it should be), then add Idisk*(ΔAdisk/ΔAcloud) to I (unobscured portion of disk)
#     #            if ΔAcloud bigger than ΔAdisk, then add 8 closest cells (in α and β) until ΔAcloud < ΔAdisk,
#     #            then repeat calculation described above for "supercell" and set intensity of all cells in supercell to this value/number of cells in supercell
#     #       --iii: set ΔA to 1.0 because when calling binModel, it will multiply by ΔA but we have already done that in weighting the cloud/disk contributions here
#     #   --b: modify velocities at each point to be weighted average (by new intensity) of all cloud/disk contributions as this will be the "observed" velocity
#     #3: after modifying existing disk ring structs, add the rest of the cloud rings unmodified (assuming clouds are small and don't overlap/low optical depth)
#     #4: reset camera α and β to match new ring structs
#     #5: done! can now call binModel to get profiles, etc. and images will be weighted correctly by area and cloud/disk contributions
    
#     #1:
#     diskFlag = [!(typeof(r.ϕ) == Float64 && typeof(r.r) == Float64) for r in m.rings] #true if disk, false if cloud
#     doneFlag = fill(false,length(diskFlag)) #flase if not done, true if done
#     newRings = deepcopy(m.rings)
#     α = m.camera.α; β = m.camera.β
#     ΔA = vcat([getVariable(m,:ΔA)...]) 
#     I = vcat([getVariable(m,:I)...])
#     v = vcat([getVariable(m,:v)...])
#     vNew = deepcopy(v)
#     rCam = sqrt.(α.^2 .+ β.^2); ϕCam = atan.(β,α)
#     rOrder = sortperm(rCam)
#     rList = unique(x->round(x,sigdigits=8),rCam)
#     ϕList = unique(x->round(x,sigdigits=8),ϕCam)
#     #make lists of lists ordered by r:
#     #at each index in r there is a corresponding list of ϕ values, ΔA values, and DiskFlag values -- order these sublists by ϕ
#     #clouds will be identified as an r entry with just one ϕ value, surrounded by either other cloud entries or disk entries
#     #cloud inds to sum up will go from [last disk ind before cloud, first disk ind after cloud(s)]
#     #calculate r and ϕ inds for disk cells that will potentially be affected by cloud inds 
#     #-- if cloud ΔA < disk ΔA (at previous disk r), then raytrace for portion of disk ΔA that is obscured by cloud ΔA and simply add unobscured portion
#     #-- if cloud ΔA > disk ΔA, then sum up 8 closest cells in α and β until cloud ΔA < disk ΔA, then raytrace for supercell and set all cells in supercell to this value
#     #2:
#     ringCount = 1
#     for i in eachindex(rOrder)
#         if diskFlag[rOrder[i]] && !doneFlag[rOrder[i]]
#             ri, ϕi = rCam[rOrder[i]], ϕCam[rOrder[i]]
#             sameCamMask = (rCam .≈ ri) .& (ϕCam .≈ ϕi) .& diskFlag #overlapping disk cells -- should usually be zero
#             if sum(sameCamMask) > 0
#                 error("Overlapping disk cells detected at camera coordinates $(ri), $(ϕi). Raytracing not supported between overlapping disk structures yet.")
#             end
#             overlap = 0
#             while !diskFlag[rOrder[i+overlap+1]]
#                 overlap += 1 
#             end

#             overlap = 0
#             if i<length(rOrder)
#                 nextInd = rOrder[i+1]
#                 while (overlap+i+1 < length(rOrder)) && (ΔA[nextInd] ≈ ΔA[rOrder[i]]) && (rCam[nextInd] ≈ rCam[rOrder[i]]) && (ϕCam[nextInd] ≈ ϕCam[rOrder[i]]) #check if camera is the same -- multiple disk cells (should usually evaluate to 0)
#                     overlap += 1
#                     nextInd = rOrder[i+overlap+1]
#                 end
#                 while !diskFlag[nextInd] && (overlap+i+1 < length(rOrder)) #go through all clouds in ΔA
#                     overlap += 1
#                     nextInd = rOrder[i+overlap+1]
#                 end
#                 nextr = rCam[nextInd]; nextϕ = ϕCam[nextInd] #next filled in disk cell in camera coordinates
#                 if nextInd == rOrder[i+1] #no overlap, just copy old intensity to new ring (but change ΔA)
#                     oldInd = findfirst(isapprox(rCam[rOrder[i]],rList)) #index of ring in m 
#                     newRings[ringCount] = m.rings[oldInd] #copy ring to newRings
#                     newRings[ringCount].I = I[rOrder[i]].*ΔA[rOrder[i]]
#                     newRings[ringCount].ΔA = 1.0
#                     doneFlag[rOrder[i]] = true
#                     ringCount += 1
#                 else #some overlap -- have to raytrace intensity at this cell
#                     #2c:
#                     #2d:
#                     #2e:
#                     #2f:
#                     #2g:
#                     #2h:
#                     #2i:
#                     #2j:
#                     #2k:
#                     #2l:
#                     #2m:
#                     #2n:
#                     #2o:
#                     #2p:
#                     #2q:
#                     #2r:
#                     #2s:
#                     #2t:
#                     #2u:
#                     #2v:
#                     #2w:
#                     #2x:
#                     #2y:
#                     #2z:
#                     doneFlag[rOrder[i]] = true
#                 end

#             end
#             #2a:
#             #2b:
#             doneFlag[order[i]] = true
#         end
#     end
#     #3:
#     for i in eachindex(diskFlag)
#         if !diskFlag[i] && !doneFlag[i]
#             #3a:
#             doneFlag[i] = true
#         end
#     end
#     #4:


#     # if isnothing(camera)
#     #     camera = m.camera
#     # else
#     #     m.camera = camera
#     # end
#     # α = camera.α
#     # β = camera.β
#     # i = getVariable(m,:i)
#     # rot = getVariable(m,:rot)
#     # θₒPoint = getVariable(m,:θₒ)
#     # xRing = @. (β*cos(rot) - α*cos(i)*sin(rot))/(cos(i)*cos(θₒPoint)+cos(rot)*sin(i)*sin(θₒPoint)) #system x
#     # yRing = @. (α*(cos(i)*cos(θₒPoint)+sec(rot)*sin(i)*sin(θₒPoint))+β*cos(θₒPoint)*tan(rot))/(cos(i)*cos(θₒPoint)*sec(rot)+sin(i)*sin(θₒPoint)) #system y
#     # r = @. √(xRing^2 + yRing^2)
#     # ϕ₀ = @. atan(yRing,xRing) #original ϕ₀ (no rotation)
#     # xyzSys = rotate3D.(r,ϕ₀,i,rot,θₒPoint) #system coordinates xyz
#     # xSys = [xyzSys[i][1] for i in 1:length(xyzSys)]
#     # ySys = [xyzSys[i][2] for i in 1:length(xyzSys)]
#     # zSys = [xyzSys[i][3] for i in 1:length(xyzSys)]
#     # img = zeros(length(α)-1,length(β)-1)
#     # I = getVariable(m,:I)
#     # τ = getVariable(m,:τ)
#     # for i in 1:length(α)-1
#     #     for j in 1:length(β)-1
#     #         rayY_l = α[i]; rayZ_b = β[j]
#     #         rayY_r = α[i+1]; rayZ_t = β[j+1]
#     #         sysMask = (xSys .>= rayY_l) .& (xSys .< rayY_r) .& (ySys .>= rayZ_b) .& (ySys .< rayZ_t)
#     #         if sum(sysMask) == 0.0
#     #             img[i,j] = 0.0
#     #         else
#     #             rayτ = 0.0
#     #             x = xSys[sysMask]; y = ySys[sysMask]; z = zSys[sysMask]
#     #             sortOrder = sortperm(x,rev=true) #highest x closest to camera
#     #             counter = 1
#     #             while rayτ < τMax && counter < length(x)
#     #                 #note that this part specificallly probably won't work because the mask is done on the system coordinates which are not neccessarily the same shape as I? check
#     #                 img[i,j] += sum(I[sysMask][sortOrder[counter]].*exp.(-rayτ)) #note -- not doing dx here, assuming they are all small and equal. could add as parameter to rings?
#     #                 rayτ += τ[sysMask][sortOrder] #τ is the additional optical depth gained *after* passing through each point
#     #                 counter += 1
#     #             end
#     #         end
#     #     end
#     # end
#     # m.camera.img = img
#     return m
# end

            
            