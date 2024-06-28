#!/usr/bin/env julia

using Plots

# function rotate3D(r,ϕ,i,rot)
#     x = r*(cos(ϕ)*cos(rot) - sin(i)*sin(ϕ)*sin(rot))
#     y = r*(cos(ϕ)*sin(rot) + sin(i)*sin(ϕ)*cos(rot))
#     z = r*cos(i)*sin(ϕ)
#     return x,y,z
# end

function rotate3D(r,ϕ,i,rot,θ0)
    x = r*(sin(i)*(cos(ϕ)*cos(θ0)*cos(rot) - sin(ϕ)*sin(rot)) + cos(i)*cos(ϕ)*sin(θ0))
    y = r*(cos(ϕ)*cos(θ0)*sin(rot) + sin(ϕ)*cos(rot))
    z = r*(cos(i)*(sin(ϕ)*sin(rot) - cos(ϕ)*cos(θ0)*cos(rot)) + sin(i)*cos(ϕ)*sin(θ0))
    return x,y,z
end

# function rotate3D(xyz,i,rot,θ0)
#     matrix = [

#     ]
#     return matrix*xyz
# end
# function invert(xyz,i,rot,θ0)
#     matrix = [cos(rot)*cos(θ0)*sin(i)+sin(θ0)*cos(i) sin(rot)*sin(i) cos(rot)*sin(θ0)*sin(i)-cos(θ0)*cos(i);
#             -sin(rot)*cos(θ0) cos(rot) -sin(rot)*sin(θ0);
#             cos(rot)*cos(θ0)*cos(i)-sin(θ0)*sin(i) sin(rot)*cos(i) cos(rot)*sin(θ0)*cos(i)+cos(θ0)*sin(i)] 
#     return matrix*xyz
# end
function rotate3DSajal(r,ϕ,i,rot,θ0)
    matrix = [
        (cos(rot)*cos(θ0)*sin(i)-cos(i)*sin(θ0)) sin(i)*sin(rot) (-cos(i)*cos(θ0)-cos(rot)*sin(i)*sin(θ0));
        -cos(θ0)*sin(rot) cos(rot) sin(rot)*sin(θ0);
        (cos(i)*cos(rot)*cos(θ0)+sin(i)*sin(θ0)) cos(i)*sin(rot) (cos(θ0)*sin(i)-cos(i)*cos(rot)*sin(θ0))
    ]
    return matrix*[r*cos(ϕ);r*sin(ϕ);0]
end

function invertSajal(xyz,i,rot,θ0)
    matrix = [
        (cos(rot)*cos(θ0)*sin(i)-cos(i)*sin(θ0)) -cos(θ0)*sin(rot) (cos(i)*cos(rot)*cos(θ0)+sin(i)*sin(θ0));
        sin(i)*sin(rot) cos(rot) cos(i)*sin(rot);
        (-cos(i)*cos(θ0)-cos(rot)*sin(i)*sin(θ0)) sin(rot)*sin(θ0) (cos(θ0)*sin(i)-cos(i)*cos(rot)*sin(θ0))
    ]
    return matrix*xyz
end

function plotRotation!(p,r,i,rot,ϕ = range(0, stop=2π, length=100))
    x,y,z = zeros(length(ϕ)),zeros(length(ϕ)),zeros(length(ϕ))
    for (j,ϕ) in enumerate(ϕ)
        x[j],y[j],z[j] = rotate3D(r,ϕ,i,rot)
    end
    p = plot!(p,x,y,z,color=1,label="")
    p = scatter!(p,[x[1]], [y[1]], [z[1]], label = "rotated ϕ = 0",color=1)
    p = scatter!(p,[r*sin(ϕ[1])],[r*cos(ϕ[1])],[0], label = "unrotated ϕ = 0",color=2)
    p = plot!(p,r.*sin.(ϕ),r.*cos.(ϕ),zeros(length(ϕ)),label = "",color=2)
    return p
end

# function taylorXDisk(α,β,i,rot)
#     return -sign(rot)*β*tan(i)
# end

# function recoverCoords(r,ϕ,i,rot) #this works, but there is numerical noise at exactly ±π/2 for rot and at ±π/2 in ϕ because of cos(ϕ) in denominator
#     #to do: taylor expansion for ϕ, rot near ±π/2? ±π/100 should be fine.
#     err1, err2, err3 = zeros(length(ϕ)), zeros(length(ϕ)), zeros(length(ϕ))
#     αList = zeros(length(ϕ)); βList = zeros(length(ϕ))
#     for (j,ϕ) in enumerate(ϕ)
#         x,y,z = rotate3D(r,ϕ,i,rot)
#         α = x; β = z
#         xDisk = (α+β*tan(i)*sin(rot))/cos(rot)
#         yDisk = β/cos(i)
#         rRecovered1 = √(xDisk^2 + yDisk^2)
#         rRecovered2 = β/(cos(i)*cos(ϕ))
#         ϕRecovered = atan(yDisk,xDisk)
#         err1[j], err2[j], err3[j] = rRecovered1-r, rRecovered2-r, ϕRecovered-ϕ
#         αList[j] = α
#         βList[j] = β
#     end
#     return err1, err2, err3, αList, βList
# end

function roots(i,θ0,rot)
    yRoots = zeros(7)
    xRoots = zeros(3)
    yRoots[1] = atan(cot(θ0)*cos(rot))
    yRoots[2] = atan(cot(i)*cos(rot))
    yRoots[3] = θ0 > atan(cot(i)) ? asec(tan(i)*tan(θ0)) : acos(tan(i)*tan(θ0))
    yRoots[4] = atan(cot(θ0)*sec(rot))
    yRoots[5] = atan(cot(i)*sec(rot))
    yRoots[6] = π/2
    yRoots[7] = -π/2
    xRoots[1] = atan(cot(θ0)*sec(rot))
    xRoots[2] = atan(cot(i)*sec(rot))
    xRoots[3] = θ0 > atan(cot(i)) ? -asec(tan(i)*tan(θ0)) : NaN
    yRoots = [(yRoots[1],yRoots[4]),(yRoots[2],yRoots[5]),(yRoots[6],yRoots[7],yRoots[3])]
    return xRoots,yRoots #(i,θ0,rot) indices
end

function recoverCoords(r,ϕ,i,rot,θ0) #this doesn't work
    err1, err2, err3 = zeros(length(ϕ)), zeros(length(ϕ)), zeros(length(ϕ))
    αList = zeros(length(ϕ)); βList = zeros(length(ϕ))
    # xRoots,yRoots = roots(i,θ0,rot)
    # rotRoots = isnan(xRoots[3]) ? [yRoots[3]...] : [xRoots[3],yRoots[3]...]
    # θ0Roots = [xRoots[2],yRoots[2]...]
    # iRoots = [xRoots[1],yRoots[1]...]
    # if sum(isapprox.(rotRoots,rot,atol=1e-8)) != 0 || sum(isapprox.(θ0Roots,θ0,atol=1e-8)) != 0 || sum(isapprox.(iRoots,i,atol=1e-8)) != 0
    #     @warn("Warning: near coordinate singularity, recovered coordinates may be erroneous")
    # end
    for (j,ϕ) in enumerate(ϕ)
        x,y,z = rotate3D(r,ϕ,i,rot,θ0)
        rRecovered1 = √(x^2 + y^2 + z^2)
        # rotateInv = [sin(θ0+i) 0 -cos(θ0+i);
        #                 0 1 0;
        #                 cos(θ0+i) 0 sin(θ0+i)]
        # x,y,z = rotateInv*[x,y,z]
        # ϕRecoverd = atan(y,x)
        # err1[j], err2[j] = rRecovered - r, ϕRecoverd - ϕ
        α = y; β = z
        xDisk = (β*cos(rot) - α*cos(i)*sin(rot))/(cos(i)*cos(θ0)+cos(rot)*sin(i)*sin(θ0))
        yDisk = (α*(cos(i)*cos(θ0)+sec(rot)*sin(i)*sin(θ0))+β*cos(θ0)*tan(rot))/(cos(i)*cos(θ0)*sec(rot)+sin(i)*sin(θ0))
        rRecovered2 = √(xDisk^2 + yDisk^2)
        ϕRecovered = atan(yDisk,xDisk)
        err1[j], err2[j], err3[j] = rRecovered1-r, rRecovered2-r, ϕRecovered-ϕ
        αList[j] = α
        βList[j] = β
    end
    return err1, err2, err3, αList, βList
end