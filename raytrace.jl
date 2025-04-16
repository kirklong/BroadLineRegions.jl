#!/usr/bin/env julia
using LinearAlgebra

function zeroDiskObscuredClouds!(m::model;diskCloudIntensityRatio::Float64=1.,rotate3D::Function=rotate3D)
    isCombined, startInds, diskFlag = detectCombinedModel(m)
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
        xyzCloud = rotate3D(ring.r,ring.ϕₒ,ring.i,ring.rot,ring.θₒ,ring.reflect) #system coordinates xyz
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
    isCombined, startInds, diskFlag = detectCombinedModel(m)
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
        xyzSys = rotate3D(m.rings[i].r,m.rings[i].ϕₒ,m.rings[i].i,m.rings[i].rot,m.rings[i].θₒ,m.rings[i].reflect) #system coordinates xyz
        xCloud = xyzSys[1]
        α,β = m.camera.α[αβStartInds[i]],m.camera.β[αβStartInds[i]]
        rCam = sqrt(α^2 + β^2)
        ϕCam = atan(β,α)
        if (rCam > rMinDisk) && (rCam < rMaxDisk)
            rDiskInd = argmin(abs.(rCam .- rUnique))
            ϕDiskInd = argmin(abs.(ϕCam .- ϕList[rDiskInd]))
            xDisk = rotate3D(m.rings[rDiskInd].r[ϕDiskInd],m.rings[rDiskInd].ϕₒ[ϕDiskInd],m.rings[rDiskInd].i,m.rings[rDiskInd].rot,m.rings[rDiskInd].θₒ,m.rings[rDiskInd].reflect)[1] #system coordinates xyz
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


function raytrace!(m::model;DiskRayFrac::Union{Float64,Nothing}=nothing,nRays::Union{Float64,Int,Nothing}=nothing,τMax::Float64=1.0,trackRays::Bool=false,rotate3D::Function=rotate3D)
    #new approach -- brute force raytracing
    isCombined, startInds, diskFlag = detectCombinedModel(m)
    if !isCombined
        @warn "did not detect combined model so nothing to raytrace -- returning unaltered input model"
        return m
    end
    if length(startInds) > 2
        #raytrace recursively combining two at a time
        mList = [deepcopy(m) for i=1:length(startInds)]
        for (i,m) in enumerate(mList)
            s = startInds[i]; e = i == length(startInds) ? length(m.rings) : startInds[i+1]-1
            m.rings = m.rings[s:e]
        end
        mFinal = mList[1]
        for i=2:length(mList)
            mFinal = mFinal + mList[i]
            mFinal = raytrace!(mFinal,DiskRayFrac,nRays,camera,τMax)
        end
        return mFinal
    end

    ΔA = vcat([getVariable(m,:ΔA)...])
    diskFlag = Array{Bool}(undef,length(ΔA))
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
    if sum(diskFlag) == 0
        @warn "no disk cells detected -- simply returning combined model with no raytracing"
        return m
    end

    if isnothing(DiskRayFrac)
        DiskRayFrac = sum(ΔA[diskFlag])/sum(ΔA[.!diskFlag]) #ratio of disk to cloud area on camera
    end
    @debug "DiskRayFrac: $DiskRayFrac"
    if typeof(nRays) == Float64
        nRays = ceil(Int,nRays)
    end
    nClouds = sum(.!diskFlag) #number of clouds
    if isnothing(nRays)
        nRays = ceil(Int,nClouds*(1+DiskRayFrac)) #total number of rays to trace
    end
    @debug "nRays: $nRays"
    sampleFrac = 1.0
    if nRays < nClouds*(1+DiskRayFrac)
        sampleFrac = nRays/(nClouds*(1+DiskRayFrac))
    end
    @debug "Sampling fraction: $sampleFrac"
    if sampleFrac < 0.1
        @warn "Number of rays is less than 10% of recommended value ($(nClouds*(1+DiskRayFrac)))-- results may be unreliable."
    end
    αₒ = vcat([m.camera.α...]); βₒ = vcat([m.camera.β...])
    diskα = αₒ[diskFlag]; diskβ = βₒ[diskFlag] 
    cloudα = αₒ[.!diskFlag]; cloudβ = βₒ[.!diskFlag]
    diskr = sqrt.(diskα.^2 .+ diskβ.^2)
    cloudr = sqrt.(cloudα.^2 .+ cloudβ.^2)
    minDiskr = minimum(diskr); maxDiskr = maximum(diskr)
    cloudMaskAlone = [(r<minDiskr || r>maxDiskr) for r in cloudr] #mask for clouds that are not in disk area
    # if isnothing(camera)
    #     camera = m.camera
    #     α = vcat(diskα,cloudα[cloudMask]); β = vcat(diskβ,cloudβ[cloudMask])
    #     camera = camera(α,β)
    #     m.camera = camera
    # else -- should allow user to specify new camera?
    #     m.camera = camera
    # end
    #not allowing new camera specification, instead just recalculating from existing camera and removing cloud cells that are inside of disk cells
    α = vcat(diskα,cloudα[cloudMaskAlone]); β = vcat(diskβ,cloudβ[cloudMaskAlone])
    cam = camera(α,β,nothing)
    m.camera = cam
    #add up intensities (with raytracing) and divide by total number of rays in each cell -- then the intensity in each cell is ∝ 1/ΔA and can weight by ΔA in binModel
    #new camera will be a "disk" camera, so every cell at the end will have ΔA = ΔADisk. This is fine even if a cloud is alone because if there are enough rays most of them will miss in said cell and the intensity will be low
    #also avoids problem of neighboring cells etc. just calculate what we hit in each cell and divide by total number of rays -- if a cloud is across 2 cells it will get hit with equal probability in both of those cells! 
    #start with cloud ray paths, then randomly sample new rays from range of r,ϕ allowed by camera. Note that this may mess up statistics if there are very few rays? Need to also have an "empty area" fraction in calculating nRays!
    #or....set ΔA in cloud only cells to sum(ΔAcloud)/sum(ΔA)/nClouds(?) (average cloud contribution to total area) -- do we need to fill in rest of "ring" at each cloud radius with corresponding number of empty ϕ cells? use sparse arrays?

    rCam = sqrt.(cam.α.^2 .+ cam.β.^2); ϕCam = atan.(cam.β,cam.α)
    rUnique = unique(r -> round(r,sigdigits=12),rCam) #note that at 8 sigfigs some clouds were actually identical and this was causing an error??
    ϕUnique = unique(ϕ -> round(ϕ,sigdigits=12),ϕCam) #note really want a list of ϕ disk values? don't care about cloud ϕ because we use this to index into rings with multiple ϕ values
    rMin = minimum(rCam); rMax = maximum(rCam)
    ϕMin = minimum(ϕCam); ϕMax = maximum(ϕCam)

    #initialize new rings
    diskFlagRing = [!(typeof(r.ϕ) == Float64 && typeof(r.r) == Float64) for r in m.rings] #true if disk ring, false if cloud ring
    ϕDisk = [m.rings[diskFlagRing][i].ϕ for i=1:sum(diskFlagRing)]
    newRingsLength = length(rUnique)
    newRings = Array{ring}(undef,newRingsLength)
    for i=1:sum(diskFlagRing) #do disk first to match rCam, ϕCam
        newRings[i] = deepcopy(m.rings[diskFlagRing][i])
    end
    for i=sum(diskFlagRing)+1:newRingsLength #do clouds next
        newRings[i] = deepcopy(m.rings[.!diskFlagRing][cloudMaskAlone][i-sum(diskFlagRing)])
        if rand() > sampleFrac
            newRings[i].I = 0.0 #zero out intensity of fraction of cloud cells that won't be sampled
        end
        #don't need to zero out inensity in general here because no actual raytracing will happen in this cell (no disk to pass through, I just stays I_cloud)
    end
    old_I = [deepcopy(r.I) for r in newRings[1:sum(diskFlagRing)]] #disk intensities before raytracing
    rayCounts = [deepcopy(r.I) for r in newRings] #counter for number of rays in each disk cell, copy I just for shape
    rays = trackRays ? Array{ray}(undef,nRays) : nothing
    for i=1:sum(diskFlagRing)
        newRings[i].I .= 0.0 #zero out intensity before raytracing
        rayCounts[i] .= 0.0 #zero out ray counters
    end

    #clouds first
    cloudInds = findall(.!diskFlag)
    cloudRingInds = findall(.!diskFlagRing)
    raysTraced = 0
    for i in cloudInds #raytrace cloud points first
        if rand() < sampleFrac
            @debug "in cloud ray loop -- total progress = $(round(raysTraced/nRays*100,sigdigits=2))%\r"
            α = αₒ[i]; β = βₒ[i]
            r = sqrt(α^2 + β^2); ϕ = atan(β,α)
            cloud_dists = sqrt.((cloudα.-α).^2 .+ (cloudβ.-β).^2) #how far is each cloud center from this ray?
            cloud_radii = sqrt.(ΔA[cloudInds]./π) #projected radius of each cloud, assuming spherical clouds
            cloudMaskRay = cloud_dists .< cloud_radii #mask for clouds that are hit by this ray
            newRingsInd = argmin(abs.(r.-rUnique)) #index of closest r value in camera, which should also match newRings index
            if r > minDiskr && r < maxDiskr #cloud is in disk area
                I = 0.0; τ = 0.0; x = 0.0
                if trackRays
                    I = zeros(1+sum(cloudMaskRay)) #initialize intensity array for ray struct
                    τ = zeros(1+sum(cloudMaskRay)) #initialize optical depth array for ray struct
                    x = zeros(1+sum(cloudMaskRay)) #initialize system x array for ray struct
                end

                #raytrace through disk area
                #need to find closest disk cell to cloud cell and raytrace through that cell
                newRingsϕ = argmin(abs.(ϕ.-ϕDisk[newRingsInd])) #closest ϕ value in camera, should match order of ϕ in newRings[r].ϕ
                disk_I = old_I[newRingsInd][newRingsϕ] #intensity of disk cell without raytracing
                τDisk = typeof(newRings[newRingsInd].τ) == Float64 ? newRings[newRingsInd].τ : newRings[newRingsInd].τ[newRingsϕ]
                xDisk = rotate3D(newRings[newRingsInd].r[newRingsϕ],newRings[newRingsInd].ϕₒ[newRingsϕ],newRings[newRingsInd].i,newRings[newRingsInd].rot,newRings[newRingsInd].θₒ)[1] #system coordinates xyz
                #xCloud = rotate3D(m.rings[i].r,m.rings[i].ϕₒ,m.rings[i].i,m.rings[i].rot,m.rings[i].θₒ,m.rings[i].reflect)[1] #system coordinates xyz
                #xClouds = rotate3D.(m.rings[cloudInds][cloudMask].r,m.rings[cloudInds][cloudMask].ϕₒ,m.rings[cloudInds][cloudMask].i,m.rings[cloudInds][cloudMask].rot,m.rings[cloudInds][cloudMask].θₒ,m.rings[cloudInds][cloudMask].reflect) #system coordinates xyz
                xClouds = [rotate3D(r.ϕ,r.ϕₒ,r.i,r.rot,r.θₒ,r.reflect)[1] for r in m.rings[cloudRingInds][cloudMaskRay]]
                
                IList = [disk_I,[r.I for r in m.rings[cloudRingInds][cloudMaskRay]]...]
                τList = [τDisk,[r.τ for r in m.rings[cloudRingInds][cloudMaskRay]]...]
                xList = [xDisk,xClouds...]
                
                xOrder = sortperm(xList,rev=true)
                ILocal = IList[xOrder][1]; τLocal = τList[xOrder][1]; xLocal = xList[xOrder][1]
                if trackRays
                    I[1] = ILocal; τ[1] = τLocal; x[1] = xLocal
                end
                counter = 1; stop = (τLocal > τMax) || (counter+1 > length(xList))
                while !stop
                    counter += 1
                    xCloud = xList[xOrder][counter]
                    ILocal += IList[xOrder][counter]*exp(-τLocal) #add intensity 
                    τLocal += τList[xOrder][counter] #update optical depth after passing through 
                    if trackRays
                        I[counter] = deepcopy(ILocal); τ[counter] = deepcopy(τLocal); x[counter] = deepcopy(xCloud)
                    end
                    stop = (τLocal > τMax) || (counter+1 > length(xList))
                end
                newRings[newRingsInd].I[newRingsϕ] += ILocal #add intensity to disk cell
                rayCounts[newRingsInd][newRingsϕ] += 1 #increment ray counter
                raysTraced += 1
                if trackRays
                    rays[raysTraced] = ray(r,ϕ,α,β,τ[1:counter],x[1:counter],I[1:counter],1)
                end
            else 
                if sum(cloudMaskRay) > 1 #ray through cloud overlaps with other clouds
                    I = 0.0; τ = 0.0; x = 0.0
                    if trackRays
                        I = zeros(sum(cloudMaskRay)) #initialize intensity array for ray struct
                        τ = zeros(sum(cloudMaskRay)) #initialize optical depth array for ray struct
                        x = zeros(sum(cloudMaskRay)) #initialize system x array for ray struct
                    end
                    xClouds = [rotate3D(r.ϕ,r.ϕₒ,r.i,r.rot,r.θₒ,r.reflect)[1] for r in m.rings[cloudInds][cloudMaskRay]]
                    IList = [r.I for r in m.rings[cloudInds][cloudMaskRay]]
                    τList = [r.τ for r in m.rings[cloudInds][cloudMaskRay]]
                    xOrder = sortperm(xClouds,rev=true)
                    ILocal = IList[xOrder][1]; τLocal = τList[xOrder][1]; xLocal = xClouds[xOrder][1]
                    if trackRays
                        I[1] = ILocal; τ[1] = τLocal; x[1] = xLocal
                    end
                    counter = 1; stop = (τLocal > τMax) || (counter+1 > length(xClouds))
                    while !stop
                        counter += 1
                        xLocal = xClouds[xOrder][counter]
                        ILocal += IList[xOrder]*exp(-τLocal) #add intensity 
                        τLocal += τList[xOrder][counter] #update optical depth after passing through 
                        if trackRays
                            I[counter] = deepcopy(ILocal); τ[counter] = deepcopy(τLocal); x[counter] = deepcopy(xLocal)
                        end
                        stop = (τLocal > τMax) || (counter+1 > length(xClouds))
                    end
                    newRings[newRingsInd].I += ILocal #add intensity to disk cell
                    rayCounts[newRingsInd] += 1 #increment ray counter
                    raysTraced += 1
                    if trackRays
                        rays[raysTraced] = ray(r,ϕ,α,β,τ[1:counter],x[1:counter],I[1:counter],2)
                    end
                else #cloud alone
                    # newRings[newRingsInd].I += m.rings[i].I #copy cloud cell intensity to new rings struct
                    # xCloud = rotate3D(m.rings[i].r,m.rings[i].ϕₒ,m.rings[i].i,m.rings[i].rot,m.rings[i].θₒ,m.rings[i].reflect)[1] #system coordinates xyz
                    # rayCounts[newRingsInd] += 1 #increment ray counter
                    # raysTraced += 1
                    # if trackRays
                    #     rays[raysTraced] = ray(r,ϕ,α,β,[[m.rings[i].τ]],[xCloud],[m.rings[i].I])
                    # end

                    if length(m.rings[cloudRingInds][cloudMaskRay]) == 1
                        if rayCounts[newRingsInd] == 0
                            println("newRings[newRingsInd] = $(newRings[newRingsInd])")
                            println("m.rings[cloudRingInds][cloudMaskRay][1] = $(m.rings[cloudRingInds][cloudMaskRay][1])")
                            newRings[newRingsInd] = deepcopy(m.rings[cloudRingInds][cloudMaskRay][1]) #copy cloud cell to new rings struct
                        else
                            newRings[newRingsInd].I += m.rings[cloudRingInds][cloudMaskRay][1].I #copy cloud cell intensity to new rings struct
                            newRings[newRingsInd].I += m.rings[cloudRingInds][cloudMaskRay][1].I #copy cloud cell intensity to new rings struct
                        end
                        rayCounts[newRingsInd] += 1 #increment ray counter
                        xCloud = rotate3D(m.rings[cloudRingInds][cloudMaskRay][1].r,m.rings[cloudRingInds][cloudMaskRay][1].ϕₒ,m.rings[cloudRingInds][cloudMaskRay][1].i,m.rings[cloudRingInds][cloudMaskRay][1].rot,m.rings[cloudRingInds][cloudMaskRay][1].θₒ,m.rings[cloudRingInds][cloudMaskRay][1].reflect)[1] #system coordinates xyz
                        raysTraced += 1
                        if trackRays
                            rays[raysTraced] = ray(r,ϕ,α,β,[m.rings[cloudRingInds][cloudMaskRay][1].τ],[xCloud],[m.rings[cloudRingInds][cloudMaskRay][1].I],3)
                        end
\                    else
                        println("CLOUD ALONE: i = $i")
                        println("typeof(newRings[newRingsInd].I) = $(typeof(newRings[newRingsInd].I))")
                        println("typeof(m.rings[cloudRingInds][cloudMaskRay][1].I) = $(typeof(m.rings[cloudRingInds][cloudMaskRay][1].I))")
                        println("r = $r -- (minDiskr,maxDiskr) = ($minDiskr,$maxDiskr)")
                        println("newRings[newRingsInd] = $(newRings[newRingsInd])")
                        println("m.rings[cloudRingInds][cloudMaskRay] = $(m.rings[cloudRingInds][cloudMaskRay])")
                        error("cloud alone mask returned more than one ring struct -- should not happen")
                    end
                end
            end
        end
    end
    #rest of system randomly draw rays
    minα = minimum(αₒ); maxα = maximum(αₒ)
    minβ = minimum(βₒ); maxβ = maximum(βₒ)
    αD = Uniform(minα,maxα); βD = Uniform(minβ,maxβ)
    for rayi=raysTraced+1:nRays #fill in rest of rays with random rays from camera
        @debug "in random ray loop -- total progress = $(round(rayi/nRays*100,sigdigits=2))%\r"
        α = rand(αD)
        β = rand(βD)
        r = sqrt(α^2 + β^2); ϕ = atan(β,α)
        cloud_dists = sqrt.((cloudα.-α).^2 .+ (cloudβ.-β).^2) #how far is each cloud center from this ray?
        cloud_radii = sqrt.(ΔA[cloudInds]./π) #projected radius of each cloud, assuming spherical clouds
        cloudMaskRay = cloud_dists .< cloud_radii #mask for clouds that are hit by this ray
        allowed_draw = (r > minDiskr && r < maxDiskr) || sum(cloudMaskRay) > 0 #no rays through empty space
        while !allowed_draw
            α = rand(Uniform(minimum(αₒ),maximum(αₒ)))
            β = rand(Uniform(minimum(βₒ),maximum(βₒ)))
            r = sqrt(α^2 + β^2); ϕ = atan(β,α)
            cloud_dists = sqrt.((cloudα.-α).^2 .+ (cloudβ.-β).^2) #
            cloudMaskRay = cloud_dists .< cloud_radii 
            allowed_draw = (r > minimum(diskr) && r < maximum(diskr)) || sum(cloudMaskRay) > 0 #no rays through empty space
        end
        newRingsInd = argmin(abs.(r.-rUnique)) #index of closest r value in camera, which should also match newRings index
        if r > minDiskr && r < maxDiskr #random ray in disk area
            newRingsϕ = argmin(abs.(ϕ.-ϕDisk[newRingsInd])) #closest ϕ value in camera, should match order of ϕ in newRings[r].ϕ
            if sum(cloudMaskRay) > 0
                I = 0.0; τ = 0.0; x = 0.0
                if trackRays
                    I = zeros(1+sum(cloudMaskRay)) #initialize intensity array for ray struct
                    τ = zeros(1+sum(cloudMaskRay)) #initialize optical depth array for ray struct
                    x = zeros(1+sum(cloudMaskRay)) #initialize system x array for ray struct
                end
                disk_I = old_I[newRingsInd][newRingsϕ] #intensity of disk cell without raytracing
                xDisk = rotate3D(newRings[newRingsInd].r[newRingsϕ],newRings[newRingsInd].ϕₒ[newRingsϕ],newRings[newRingsInd].i,newRings[newRingsInd].rot,newRings[newRingsInd].θₒ)[1] #system coordinates xyz
                xClouds = [rotate3D(r.ϕ,r.ϕₒ,r.i,r.rot,r.θₒ,r.reflect)[1] for r in m.rings[cloudRingInds][cloudMaskRay]]
                
                IList = [disk_I,[r.I for r in m.rings[cloudRingInds][cloudMaskRay]]...]
                τDisk = typeof(newRings[newRingsInd].τ) == Float64 ? newRings[newRingsInd].τ : newRings[newRingsInd].τ[newRingsϕ]
                τList = [τDisk,[r.τ for r in m.rings[cloudRingInds][cloudMaskRay]]...]
                xList = [xDisk,xClouds...]
                
                xOrder = sortperm(xList,rev=true)
                ILocal = IList[xOrder][1]; τLocal = τList[xOrder][1]; xLocal = xList[xOrder][1]
                if trackRays
                    I[1] = ILocal; τ[1] = τLocal; x[1] = xLocal
                end
                counter = 1; stop = (τLocal > τMax) || (counter+1 > length(xList))
                while !stop
                    counter += 1
                    xLocal = xList[xOrder][counter]
                    ILocal += IList[xOrder][counter]*exp(-τLocal) #add intensity 
                    τLocal += τList[xOrder][counter] #update optical depth after passing through 
                    if trackRays
                        I[counter] = deepcopy(ILocal); τ[counter] = deepcopy(τLocal); x[counter] = deepcopy(xLocal)
                    end
                    stop = (τLocal > τMax) || (counter+1 > length(xList))
                end
                newRings[newRingsInd].I[newRingsϕ] += ILocal #add intensity to disk cell
                rayCounts[newRingsInd][newRingsϕ] += 1 #increment ray counter
                raysTraced += 1
                if trackRays
                    rays[raysTraced] = ray(r,ϕ,α,β,τ[1:counter],x[1:counter],I[1:counter],4)
                end
            else
                newRings[newRingsInd].I[newRingsϕ] += old_I[newRingsInd][newRingsϕ] #add disk cell intensity to new rings struct, no raytracing needed
                rayCounts[newRingsInd][newRingsϕ] += 1 #increment ray counter
                xDisk = rotate3D(newRings[newRingsInd].r[newRingsϕ],newRings[newRingsInd].ϕₒ[newRingsϕ],newRings[newRingsInd].i,newRings[newRingsInd].rot,newRings[newRingsInd].θₒ)[1] #system coordinates xyz
                τDisk = typeof(newRings[newRingsInd].τ) == Float64 ? newRings[newRingsInd].τ : newRings[newRingsInd].τ[newRingsϕ]
                raysTraced += 1
                if trackRays
                    rays[raysTraced] = ray(r,ϕ,α,β,[τDisk],[xDisk],[newRings[newRingsInd].I[newRingsϕ]],5)
                end
            end
        else #ray through at least one cloud, not in disk area
            if sum(cloudMaskRay) > 1 #ray through cloud overlaps with other clouds
                I = 0.0; τ = 0.0; x = 0.0
                if trackRays
                    I = zeros(sum(cloudMaskRay)) #initialize intensity array for ray struct
                    τ = zeros(sum(cloudMaskRay)) #initialize optical depth array for ray struct
                    x = zeros(sum(cloudMaskRay)) #initialize system x array for ray struct
                end
                xClouds = [rotate3D(r.ϕ,r.ϕₒ,r.i,r.rot,r.θₒ,r.reflect)[1] for r in m.rings[cloudRingInds][cloudMaskRay]]
                IList = [r.I for r in m.rings[cloudRingInds][cloudMaskRay]]
                τList = [r.τ for r in m.rings[cloudRingInds][cloudMaskRay]]
                xOrder = sortperm(xClouds,rev=true)
                ILocal = IList[xOrder][1]; τLocal = τList[xOrder][1]; xLocal = xClouds[xOrder][1]
                if trackRays
                    I[1] = ILocal; τ[1] = τLocal; x[1] = xLocal
                end
                counter = 1; stop = (τLocal > τMax) || (counter+1 > length(xClouds))
                while !stop
                    counter += 1
                    xLocal = xList[xOrder][counter]
                    ILocal += IList[xOrder]*exp(-τLocal) #add intensity 
                    τLocal += τList[xOrder][counter] #update optical depth after passing through 
                    if trackRays
                        I[counter] = deepcopy(ILocal); τ[counter] = deepcopy(τLocal); x[counter] = deepcopy(xLocal)
                    end
                    stop = (τLocal > τMax) || (counter+1 > length(xList))
                end
                newRings[newRingsInd].I += ILocal #add intensity to disk cell
                rayCounts[newRingsInd] += 1 #increment ray counter
                raysTraced += 1
                if trackRays
                    rays[raysTraced] = ray(r,ϕ,α,β,τ[1:counter],x[1:counter],I[1:counter],6)
                end
            else #cloud alone
                if length(m.rings[cloudRingInds][cloudMaskRay]) == 1
                    newRings[newRingsInd].I += m.rings[cloudRingInds][cloudMaskRay][1].I #copy cloud cell intensity to new rings struct
                    newRings[newRingsInd].I += m.rings[cloudRingInds][cloudMaskRay][1].I #copy cloud cell intensity to new rings struct
                    rayCounts[newRingsInd] += 1 #increment ray counter
                    xCloud = rotate3D(m.rings[cloudRingInds][cloudMaskRay][1].r,m.rings[cloudRingInds][cloudMaskRay][1].ϕₒ,m.rings[cloudRingInds][cloudMaskRay][1].i,m.rings[cloudRingInds][cloudMaskRay][1].rot,m.rings[cloudRingInds][cloudMaskRay][1].θₒ,m.rings[cloudRingInds][cloudMaskRay][1].reflect)[1] #system coordinates xyz
                    raysTraced += 1
                    if trackRays
                        rays[raysTraced] = ray(r,ϕ,α,β,[m.rings[cloudRingInds][cloudMaskRay][1].τ],[xCloud],[m.rings[cloudRingInds][cloudMaskRay][1].I],7)
                    end
                else
                    error("cloud alone mask returned more than one ring struct -- should not happen")
                end
            end
        end
    end
    for i=1:newRingsLength
        if typeof(newRings[i].I) == Float64
            newRings[i].I /= rayCounts[i] #cloud, divide by number of rays in each cell
        else
            newRings[i].I ./= rayCounts[i] #broadcast, arrays should be same size from setup
        end
    end
    m.rings = newRings
    reset!(m)
    if trackRays
        m.camera.rays = rays
    end
    return m
end



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
#     # ϕₒ = @. atan(yRing,xRing) #original ϕₒ (no rotation)
#     # xyzSys = rotate3D.(r,ϕₒ,i,rot,θₒPoint) #system coordinates xyz
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

            
            