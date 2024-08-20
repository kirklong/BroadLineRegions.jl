#!/usr/bin/env julia
using LinearAlgebra

function raytrace!(m::model;camera::Union{camera,nothing}=nothing,τMax=1.0)
    #don't know if this actually works, but general idea is right
    if isnothing(camera)
        camera = m.camera
    end
    α = camera.α
    β = camera.β
    i = getVariable(m,:i)
    rot = getVariable(m,:rot)
    θₒPoint = getVariable(m,:θₒ)
    xRing = @. (β*cos(rot) - α*cos(i)*sin(rot))/(cos(i)*cos(θₒPoint)+cos(rot)*sin(i)*sin(θₒPoint)) #system x
    yRing = @. (α*(cos(i)*cos(θₒPoint)+sec(rot)*sin(i)*sin(θₒPoint))+β*cos(θₒPoint)*tan(rot))/(cos(i)*cos(θₒPoint)*sec(rot)+sin(i)*sin(θₒPoint)) #system y
    r = @. √(xRing^2 + yRing^2)
    ϕₒ = @. atan(yRing,xRing) #original ϕₒ (no rotation)
    xyzSys = rotate3D.(r,ϕₒ,i,rot,θₒPoint) #system coordinates xyz
    xSys = [xyzSys[i][1] for i in 1:length(xyzSys)]
    ySys = [xyzSys[i][2] for i in 1:length(xyzSys)]
    zSys = [xyzSys[i][3] for i in 1:length(xyzSys)]
    img = zeros(size(α))
    I = getVariable(m,:I)
    τ = getVariable(m,:τ)
    for i in 1:size(α)[1]-1
        for j in 1:size(α)[2]-1
            rayY_l = α[i,j]; rayZ_b = β[i,j]
            rayY_r = α[i,j+1]; rayZ_t = β[i+1,j]
            sysMask = (xSys .>= rayY_l) .& (xSys .< rayY_r) .& (ySys .>= rayZ_b) .& (ySys .< rayZ_t)
            if sum(sysMask) == 0.0
                img[i,j] = 0.0
            else
                rayτ = 0.0
                x = xSys[sysMask]; y = ySys[sysMask]; z = zSys[sysMask]
                sortOrder = sortperm(x,rev=true) #highest x closest to camera
                counter = 1
                while rayτ < τMax && counter < length(x)
                    #note that this part specificallly probably won't work because the mask is done on the system coordinates which are not neccessarily the same shape as I? check
                    img[i,j] += sum(I[sysMask][sortOrder[counter]].*exp.(-rayτ))*ΔA[sysMask][sortOrder[counter]] #ΔA ~ dx (should probably explicitly have dx as model param like τ?)
                    rayτ += τ[sysMask][sortOrder] #τ is the additional optical depth gained *after* passing through each point
                    counter += 1
                end
            end
        end
    end
    m.camera.img = img
    return m
end

            
            