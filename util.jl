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

@userplot  Image #note that then to call this method use lowercase, i.e. image(m,"I") -- PROBLEM: this doesn't actually loop through each x and y point -- need to collapse them into 1D arrays? 
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

@userplot Plot3d #note that then to call this method use lowercase, i.e. plot3d(m,"I") #note, something doesn't work here -- check tomorrow
@recipe function f(p::Plot3d)
    xlabel --> "x [rₛ]"
    ylabel --> "y [rₛ]"
    zlabel --> "z [rₛ]"
    aspect_ratio --> :equal
    model, variable = nothing, nothing
    if length(p.args) == 2
        model, variable = p.args
    elseif length(p.args) == 1
        model = p.args[1]
    else
        throw(ArgumentError("expected 1 or 2 arguments, got $(length(p.args))"))
    end
    title --> (isnothing(variable) ? "System geometry visualization" : "System geometry + $variable (color) visualization")
    i = getVariable(model,:i)
    xtmp = deepcopy(model.camera.α)
    ytmp = deepcopy(model.camera.β)
    ztmp = zeros(size(xtmp))
    if typeof(model.camera.α) == Matrix{Float64}
        if length(size(i)) == 1 
            for j in 1:length(i)
                ytmp[j,:]./=cos(i[j])
            end
        else
            ytmp./=cos.(i)
        end
        for ri in 1:size(xtmp)[1]
            for ϕi in 1:size(xtmp)[2]
                xyz = [xtmp[ri,ϕi];ytmp[ri,ϕi];ztmp[ri,ϕi]]
                inc = length(size(i)) == 1 ? i[ri] : i[ri,ϕi]
                rot = [1 0 0; 0 cos(π/2-inc) -sin(π/2-inc); 0 sin(π/2-inc) cos(π/2-inc)] #rotate disk coords CCW about x axis by 90 - i deg such that angle between plot z and disk plane is i (draw it out)
                # this means the "camera" is now staring down the y axis of the plot system
                xtmp[ri,ϕi],ytmp[ri,ϕi],ztmp[ri,ϕi] = rot*xyz #have to do this way because r,ϕ are in disk plane (so can't just do spherical coordinate conversion)
            end
        end
    else
        for j in 1:length(i)
            ytmp[j]/=cos(i[j])
        end
        for j in 1:length(xtmp)
            xyz = [xtmp[j];ytmp[j];ztmp[j]]
            inc = i[j]
            rot = [1 0 0; 0 cos(π/2-inc) -sin(π/2-inc); 0 sin(π/2-inc) cos(π/2-inc)] #rotate disk coords CCW about x axis by 90 - i deg such that angle between plot z and disk plane is i (draw it out)
            # this means the "camera" is now staring down the y axis of the plot system
            xtmp[j],ytmp[j],ztmp[j] = rot*xyz #have to do this way because r,ϕ are in disk plane (so can't just do spherical coordinate conversion)
        end
    end
    @series begin
        r = maximum([maximum(model.camera.α),maximum(model.camera.β)])
        θ = range(0,stop=2π,length=64)
        seriestype := :path
        linecolor --> :dodgerblue
        linewidth --> 1.0
        fillalpha --> 0.1
        label --> "camera" 
        x := vec(r.*cos.(θ))
        z := vec(r.*sin.(θ))
        y := vec(ones(length(θ)).*maximum(vec(ytmp)))
        ()
    end
    seriestype := :scatter
    marker_z := isnothing(variable) ? nothing : vec(getVariable(model,variable))
    markerstrokewidth --> 0.0
    markersize --> 1.
    label --> ""
    x := vec(xtmp)
    y := vec(ytmp)
    z := vec(ztmp)
    ()
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