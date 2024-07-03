#!/usr/bin/env julia
using RecipesBase
"""retrive elements from model object and stack them into matrices for easy manipulation
params:
    m: model
        - model object to extract variables from
    variable: Union{String, Symbol, Function}
        - variable to extract from model (if String will be converted to Symbol)
        - must be a valid attribute of model.rings (e.g. :I, :v, :r, :e, :i, :ϕ) or a function that can be applied to model.rings
            - example: Keplerian disk time delays could be calculated like `t(ring) = ring.r*(1 .+ sin.(ring.ϕ).*ring.i))`
returns:
    y: Matrix{Float64}
        - matrix of extracted variable from model.rings, created by stacking the output variable for each ring
        - for example, if variable given is :I y will have shape (length(r), length(ϕ)) as at each r and ϕ there is a value of I
"""
function getVariable(m::model,variable::String) # method for getting variable if String 
    variable = Symbol(variable)
    if variable ∉ fieldnames(ring)
        throw(ArgumentError("variable must be a valid attribute of model.rings\nvalid attributes: $(fieldnames(ring))"))
    end
    return stack([getfield(ring,variable) for ring in m.rings],dims=1)
end,
function getVariable(m::model,variable::Symbol) # method for getting variable if Symbol
    if variable ∉ fieldnames(ring)
        throw(ArgumentError("variable must be a valid attribute of model.rings\nvalid attributes: $(fieldnames(ring))"))
    end
    return stack([getfield(ring,variable) for ring in m.rings],dims=1)
end,
function getVariable(m::model,variable::Function) # method for getting variable if Function
    y = nothing
    try
        y = [variable(ring) for ring in m.rings]
    catch
        throw(error("error in function call $(variable)(ring)"))
    end
    return stack(y,dims=1)
end

@userplot Image #note that then to call this method use lowercase, i.e. image(m,"I") -- PROBLEM: this doesn't actually loop through each x and y point -- need to collapse them into 1D arrays? 
@recipe function f(img::Image)
    model, variable = nothing, nothing
    if length(img.args) == 2
        model, variable = img.args
    else
        throw(ArgumentError("expected 2 arguments (model, variable), got $(img.args)"))
    end
    model, variable = img.args
    z = vec(getVariable(model,variable))
    seriestype := :scatter
    marker_z := z
    markerstrokewidth --> 0.0
    markersize --> 1.
    x := vec(model.camera.α)
    y := vec(model.camera.β)
    xlabel --> "α [rₛ]"
    ylabel --> "β [rₛ]"
    title --> "Image of $variable"
    aspect_ratio --> :equal
    label --> false
    ()
end

function get_r3D(i,rot,θₒ)
    """calculate rotation matrix to transform from initial XY plane coordinates to 3D space
    params:
        i: inclination angle of ring (rad) {Float64}
        rot: rotation of ring plane about z axis (rad) {Float64}
        θₒ: opening angle of point {Float64}
    returns:
        matrix: 3x3 rotation matrix {Matrix{Float64}}
    """
    matrix = [
        (cos(rot)*cos(θₒ)*sin(i)-cos(i)*sin(θₒ)) sin(i)*sin(rot) (-cos(i)*cos(θₒ)-cos(rot)*sin(i)*sin(θₒ));
        -cos(θₒ)*sin(rot) cos(rot) sin(rot)*sin(θₒ);
        (cos(i)*cos(rot)*cos(θₒ)+sin(i)*sin(θₒ)) cos(i)*sin(rot) (cos(θₒ)*sin(i)-cos(i)*cos(rot)*sin(θₒ))
    ]
    return matrix
end
function reflect!(xyzSys,i)
    midSlope = cot(i)
    den = 1+midSlope^2
    xf = (xyzSys[1]-midSlope^2*xyzSys[1]+2*midSlope*xyzSys[3])/den
    zf = (2*midSlope*xyzSys[1]+(midSlope^2-1)*xyzSys[3])/den
    xyzSys[1] = xf
    xyzSys[3] = zf
    return xyzSys
end
function rotate3D(r,ϕₒ,i,rot,θₒ,reflect=false)
    """go from ring coordinates to 3D coordinates where camera is at +x 
    params:
        r: radius from central mass (in terms of rₛ) {Float64}
        ϕₒ: starting azimuthal angle in ring plane (rad) {Float64}
        i: inclination angle of ring plane (rad) {Float64}
        rot: rotation of system plane about z axis (rad) {Float64}
    returns: [x;y;z]
        x: x coordinate in 3D space {Float64}
        y: y coordinate in 3D space {Float64}
        z: z coordinate in 3D space {Float64}
    """
    matrix = get_r3D(i,rot,θₒ)
    xyzSys = matrix*[r*cos(ϕₒ);r*sin(ϕₒ);0] 
    if reflect
        xyzSys = DiskWind2.reflect!(xyzSys,i)
    end
    return xyzSys
end

midPlaneXZ(x,i) = x*cot(i)

@userplot Plot3d #note that then to call this method use lowercase, i.e. plot3d(m,"I") 
@recipe function f(p::Plot3d)
    xlabel --> "x [rₛ]"
    ylabel --> "y [rₛ]"
    zlabel --> "z [rₛ]"
    aspect_ratio --> :equal
    model, variable, annotatedCamera = nothing, nothing, true
    if length(p.args) == 2
        model, tmp = p.args
        if typeof(tmp) == Bool
            annotatedCamera = tmp
        else
            variable = tmp
        end
    elseif length(p.args) == 1
        model = p.args[1]
    elseif length(p.args) == 3
        model, tmp1, tmp2 = p.args
        if typeof(tmp1) == Bool
            annotatedCamera = tmp1
            variable = tmp2
        else
            variable = tmp1
            annotatedCamera = tmp2
        end
    else
        throw(ArgumentError("expected 1, 2, or 3 arguments, got $(length(p.args))"))
    end
    title --> (isnothing(variable) ? "System geometry visualization" : "System geometry + $variable (color) visualization")
    i = getVariable(model,:i)
    r = getVariable(model,:r)
    ϕₒ = getVariable(model,:ϕₒ)
    xtmp = zeros(size(ϕₒ))
    ytmp = zeros(size(ϕₒ))
    ztmp = zeros(size(ϕₒ))
    if typeof(r) == Matrix{Float64} && typeof(ϕₒ) == Matrix{Float64}
        for ii in 1:size(r)[1]
            for jj in 1:size(r)[2]
                rot = model.rings[ii].rot
                xtmp[ii,jj],ytmp[ii,jj],ztmp[ii,jj] = rotate3D(r[ii,jj],ϕₒ[ii,jj],i[ii],rot,model.rings[ii].θₒ,model.rings[ii].reflect) 
            end
        end
    elseif typeof(r) == Vector{Float64} && typeof(ϕₒ) == Matrix{Float64}
        for ii in 1:size(ϕ)[1]
            for jj in 1:size(ϕ)[2]
                rot = model.rings[ii].rot
                xtmp[ii,jj],ytmp[ii,jj],ztmp[ii,jj] = rotate3D(r[ii],ϕₒ[ii,jj],i[ii],rot,model.rings[ii].θₒ,model.rings[ii].reflect) 
            end
        end
    else #if r is just a vector (with ϕ and i matching)
        for ii in 1:length(r)
            rot = model.rings[ii].rot
            xtmp[ii],ytmp[ii],ztmp[ii] = rotate3D(r[ii],ϕₒ[ii],i[ii],rot,model.rings[ii].θₒ,model.rings[ii].reflect) 
        end
    end
    boxSize = 1.1*maximum([maximum(i for i in xtmp if !isnan(i)),maximum(i for i in ytmp if !isnan(i)),maximum(i for i in ztmp if !isnan(i))])
    ylims --> (-boxSize,boxSize)
    xlims --> (-boxSize,boxSize)
    zlims --> (-boxSize,boxSize)
    foreground_color_legend --> nothing
    @series begin
        subplot := 1
        seriestype := :scatter
        palette --> :magma
        marker_z := (isnothing(variable) ? nothing : vec(getVariable(model,variable)))
        markerstrokewidth --> 0.0
        markersize --> 1.
        label --> ""
        x := vec(xtmp)
        y := vec(ytmp)
        z := vec(ztmp)
        ()
    end
    if annotatedCamera
        r = sqrt(maximum(model.camera.α.^2 .+ model.camera.β.^2))
        @series begin
            subplot := 1
            θ = range(0,stop=2π,length=64)
            seriestype := :path
            linestyle --> :dash
            linecolor --> :dodgerblue
            linewidth --> 1.0
            fillalpha --> 0.1
            label --> "" 
            x := vec(ones(length(θ)).*boxSize)
            z := vec(r.*sin.(θ))
            y := vec(r.*cos.(θ))
            ()
        end
        @series begin
            subplot := 1
            x := [0,boxSize]
            y := [0,0]
            z := [0,0]
            seriestype := :path
            linecolor --> :dodgerblue
            linewidth --> 1.0
            linestyle --> :dash
            arrow := :arrow
            label --> "camera"
            ()
        end
        @series begin
            subplot := 1
            x := [boxSize,boxSize]
            y := [r,-r]
            z := [0,0]
            seriestype := :path
            linestyle --> :dash
            linecolor --> :dodgerblue
            linewidth --> 1.0
            label --> ""
            ()
        end
        @series begin
            subplot := 1
            x := [boxSize,boxSize]
            y := [0,0]
            z := [r,-r]
            seriestype := :path
            linestyle --> :dash
            linecolor --> :dodgerblue
            linewidth --> 1.0
            label --> ""
            ()
        end
        @series begin
            subplot := 1 
            x := [-boxSize,-boxSize,boxSize,boxSize]
            y := [-boxSize,boxSize,boxSize,-boxSize]
            z := [midPlaneXZ(-boxSize,i[1]),midPlaneXZ(-boxSize,i[1]),midPlaneXZ(boxSize,i[1]),midPlaneXZ(boxSize,i[1])]
            seriestype := :path
            color --> :crimson
            fill --> true
            fillalpha --> 0.1
            label --> "midplane"
            ()
        end
    end
end

@userplot Profile 
@recipe function f(p::Profile)
    m, variable = nothing, nothing
    if length(p.args) == 2
        m, variable = p.args
        variable = Symbol(variable)
    else
        if typeof(p.args[1]) == model
            m = p.args[1]
        else
            error("expected arguments (model, variable[optional]), got $(p.args)")
        end
    end
    if length(m.profiles) == 0
        error("no profiles set in model")
    end
    if isnothing(variable)
        variable = collect(keys(m.profiles))
    elseif typeof(variable) == Symbol
        variable = [variable]
    end
    title --> (length(variable) == 1 ? "Profile of $(variable[1])" : "Model profiles")
    ylabel --> (length(variable) == 1 ? "$(variable[1])" : "normalized value")
    xlabel --> "Δv [c]"
    for (i,v) in enumerate(variable)
        norm = length(variable) == 1 ? 1.0 : maximum(i for i in m.profiles[v].binSums if !isnan(i))
        @series begin
            subplot := 1
            seriestype := :path
            x := m.profiles[v].binCenters
            y := m.profiles[v].binSums./norm
            marker --> :circle 
            markerstrokewidth --> 0.0
            markersize --> 2. 
            color --> i
            label --> (length(variable) == 1 ? "" : "$v (max = $(round(norm, sigdigits=3)))")
            ()
        end
    end
end



# function image(m::model,variable::Union{String,Symbol,Function}; plot=true, kwargs...)
#     y = getVariable(m,variable)
#     if plot

#         p = f(m,y;kwargs...)
#         return (y,p)
#     else
#         return y
#     end
# end