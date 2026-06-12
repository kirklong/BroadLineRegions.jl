#!/usr/bin/env julia
using Plots, RecipesBase

"""
    reset!(m::model; profiles=true, img=false)

Erase existing profiles/raytrace status and invalidate the model's `getVariable` cache.

Call this after mutating ring fields directly (e.g. `m.rings[1].I .= ...`) so that subsequent
`getVariable`/`getProfile` calls see the new values instead of memoized ones.

# Parameters
- `m::model`: Model object to reset
- `profiles::Bool=true`: If true, reset profiles
- `img::Bool=false`: If true, reset raytracing boolean (does not change existing model but allows model to be raytraced again after combining other new models)
"""
function reset!(m::model;profiles=true,img=false)
    if profiles
        m.profiles = Dict{Symbol,profile}()
    end
    if img
        m.camera.raytraced = false
    end
    if !isnothing(m.cache)
        empty!(m.cache) #invalidate memoized getVariable results
    end
end

"""
    removeNaN!(m::model) 

Remove points with `I = NaN` from model.

# Parameters
- `m::model`: Model to remove points from

# Returns
- `m::model`: Model with NaN points removed
"""
function removeNaN!(m::model)
    if !isnothing(m.cache)
        empty!(m.cache) #defend against a stale cache if rings were mutated without reset! -- gather fresh below
    end
    I = getVariable(m,:I,flatten=true)
    NaNMask = .!isnan.(I)
    #remove camera points with I = 0.0
    α = m.camera.α[NaNMask]; β = m.camera.β[NaNMask]
    m.camera = camera(α,β,m.camera.raytraced) #update camera with new α and β
    #remove rings with I = 0.0
    ringMask = Array{Bool}(undef,length(m.rings))
    for i in 1:length(m.subModelStartInds)
        endInd = i == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[i+1]-1
        subRings = m.rings[m.subModelStartInds[i]:endInd]
        #first go through subrings and remove NaN clouds (typeof(I) == Float64)
        for (j,ring) in enumerate(subRings)
            if typeof(ring.I) == Float64 && isnan(ring.I)
                ringMask[m.subModelStartInds[i]+j-1] = false #if I is cloud and NaN, don't keep this ring
            else
                ringMask[m.subModelStartInds[i]+j-1] = true #if I is not NaN and type of ring not cloud, keep this ring
            end
        end
        subRingMask = ringMask[m.subModelStartInds[i]:endInd] #get mask for subrings
        subRings = subRings[subRingMask] #filter out cloud rings with NaN I
        for ring in subRings #remove NaN points from continous rings
            NaNMask = .!isnan.(ring.I) 
            ring.v = ring.v[NaNMask] #keep only not NaN velocities
            ring.ϕ = ring.ϕ[NaNMask] #keep only not NaN azimuthal angles
            ring.r = ring.r[NaNMask] #keep only not NaN radii
            ring.ΔA = ring.ΔA[NaNMask] #keep only not NaN area elements
            if length(ring.i) == length(ring.I)
                ring.i = ring.i[NaNMask] #keep only not NaN inclinations
            elseif sum(NaNMask) != 0 #some points are NaN but i is not array (yet) -- convert to array copying values in place except for where NaN
                ring.i = fill(ring.i,length(ring.I))
                ring.i = ring.i[NaNMask] #keep only not NaN inclinations
            end
            #if sum(NaNMask) == 0 no points are NaN so don't modify anything as nothing will be removed
            if length(ring.τ) == length(ring.I)
                ring.τ = ring.τ[NaNMask] #keep only not NaN optical depths
            elseif sum(NaNMask) != 0 #some points are NaN but τ is not array (yet)
                ring.τ = fill(ring.τ,length(ring.I))
                ring.τ = ring.τ[NaNMask] #keep only not NaN optical depths
            end
            if length(ring.η) == length(ring.I)
                ring.η = ring.η[NaNMask] #keep only not NaN η values
            elseif sum(NaNMask) != 0 #some points are NaN but η is not array (yet)
                ring.η = fill(ring.η,length(ring.I))
                ring.η = ring.η[NaNMask] #keep only not NaN η values
            end
            if length(ring.ϕ₀) == length(ring.I)
                ring.ϕ₀ = ring.ϕ₀[NaNMask] #keep only not NaN starting azimuthal angles
            elseif sum(NaNMask) != 0 #some points are NaN but ϕ₀ is not array (yet)
                ring.ϕ₀ = fill(ring.ϕ₀,length(ring.I))
                ring.ϕ₀ = ring.ϕ₀[NaNMask] #keep only not NaN starting azimuthal angles
            end
            if length(ring.rot) == length(ring.I)
                ring.rot = ring.rot[NaNMask] #keep only not NaN rotation angles
            elseif sum(NaNMask) != 0 #some points are NaN but rot is not array (yet)
                ring.rot = fill(ring.rot,length(ring.I))
                ring.rot = ring.rot[NaNMask] #keep only not NaN rotation angles
            end
            if length(ring.θₒ) == length(ring.I)
                ring.θₒ = ring.θₒ[NaNMask] #keep only not NaN opening angles
            elseif sum(NaNMask) != 0 #some points are NaN but θₒ is not array (yet)
                ring.θₒ = fill(ring.θₒ,length(ring.I))
                ring.θₒ = ring.θₒ[NaNMask] #keep only not NaN opening angles
            end
            if length(ring.reflect) == length(ring.I)
                ring.reflect = ring.reflect[NaNMask] #keep only not NaN reflect flags
            elseif sum(NaNMask) != 0 #some points are NaN but reflect is not array (yet)
                ring.reflect = fill(ring.reflect,length(ring.I))
                ring.reflect = ring.reflect[NaNMask] #keep only not NaN reflect flags
            end
            #cached coordinates are per-point: mask like the other fields, but if the cache is in an
            #unexpected shape just drop it (getXYZ will lazily recompute from the masked geometry)
            if !isnothing(ring.x) && length(ring.x) == length(ring.I)
                ring.x = ring.x[NaNMask] #keep only not NaN x coordinates
            elseif !isnothing(ring.x) && sum(NaNMask) != 0
                ring.x = nothing
            end
            if !isnothing(ring.y) && length(ring.y) == length(ring.I)
                ring.y = ring.y[NaNMask] #keep only not NaN y coordinates
            elseif !isnothing(ring.y) && sum(NaNMask) != 0
                ring.y = nothing
            end
            if !isnothing(ring.z) && length(ring.z) == length(ring.I)
                ring.z = ring.z[NaNMask] #keep only not NaN z coordinates
            elseif !isnothing(ring.z) && sum(NaNMask) != 0
                ring.z = nothing
            end
            ring.I = ring.I[NaNMask] #keep only not NaN intensities
        end
    end
    m.rings = m.rings[ringMask] #filter out cloud NaN rings
    m.rings =  filter(ring -> length(ring.I) > 0, m.rings) #filter out rings with no points left
    m.cache = Dict{Any,Array}() #invalidate (and re-enable, if disabled) the getVariable cache -- ring contents changed
    oldLength = length(m.rings[1].I) #get length of first ring to use as reference
    m.subModelStartInds = [1] #reset subModelStartInds to only include first ring
    for (i,ring) in enumerate(m.rings[2:end])
        if length(ring.I) != oldLength
            m.subModelStartInds = push!(m.subModelStartInds,i+1) #if the length of the ring is not equal to the first ring, add it to the subModelStartInds
            oldLength = length(ring.I) #update oldLength to the new length
        end
    end
    return m
end

"""
    getFlattenedCameraIndices(m::model)

Get flattened camera indices corresponding to rings in model.

# Arguments
- `m::model`: Model object to extract camera indices from

# Returns
- `camStartInds::Vector{Int}`: Vector of camera starting indices with length equal to `m.subModelStartInds`
"""
function getFlattenedCameraIndices(m::model)
    camStartInds = Array{Int}(undef,length(m.subModelStartInds))
    camStartInds[1] = 1
    for i=2:length(m.subModelStartInds)
        nPerRing = length(m.rings[m.subModelStartInds[i-1]].ϕ) #could also use r, just need to check how many points are in each ring 
        nRings = m.subModelStartInds[i]-m.subModelStartInds[i-1]
        camStartInds[i] = camStartInds[i-1] + nPerRing*nRings
    end
    return camStartInds
end

"""
    getRingFromFlattenedInd(m::model, flattenedInd::Int) -> Tuple{Int, Int}

Retrieve the model ring index and subindex (if applicable) from flattened array index.

# Arguments
- `m::model`: Model object with rings 
- `flattenedInd::Int`: The index in the flattened array we need to work back from

# Returns
- `row::Int`: The ring index in `model.rings` that the flattened index corresponds to
- `column::Int`: The subindex that matches the flattened index passed to this function.
"""
function getRingFromFlattenedInd(m::model,flattenedInd::Int)
    #need to do it per chunk 
    row = 0; column = 0
    if length(m.subModelStartInds) == 1 
        column = floor(Int,flattenedInd/length(m.rings)+1)
        row = floor(Int,flattenedInd%length(m.rings))
        if row == 0
            row = length(m.rings) #when there is no remainder that is the bottom of column
            column -= 1 #no remainder = bottom of column, but flattened/length(m.rings)+1 will give column + 1
        end
    else #if there are multiple submodels, do the process above but for just the specific submodel, then add the offset of all previous submodels 
        subModelFlattenedInds = getFlattenedCameraIndices(m)
        subInd = findfirst(flattenedInd .< subModelFlattenedInds)
        if isnothing(subInd)
            subInd = length(subModelFlattenedInds)
        else
            subInd -= 1 #findfirst returns the first index where the condition is true, so we need to subtract 1 to get the correct subindex
        end
        subRings = subInd == length(subModelFlattenedInds) ? m.rings[m.subModelStartInds[subInd]:end] : m.rings[m.subModelStartInds[subInd]:m.subModelStartInds[subInd+1]-1]
        subDifference = flattenedInd - subModelFlattenedInds[subInd] + 1 #difference between flattened index and the start of the submodel, +1 for 1-based indexing
        column = floor(Int,(subDifference)/length(subRings)+1)
        row = floor(Int,(subDifference)%length(subRings))
        if row == 0
            row = length(subRings) #when there is no remainder that is the bottom of column
            column -= 1 #no remainder = bottom of column, but flattened/length(m.rings)+1 will give column + 1
        end
        row += m.subModelStartInds[subInd] - 1 #add the starting index of the submodel to the ring index
    end
    return row, column #ringInd is row, subInd is column, for clouds/point models column should always be 1
end

"""
    getVariable(m::model, variable::Union{String,Symbol,Function}; flatten=false)

Retrieve elements from model object and stack them into matrices for easy manipulation.

# Arguments
- `m::model`: Model object to extract variables from
- `variable::Union{String,Symbol,Function}`: Variable to extract from model
  - If `String`, will be converted to `Symbol`
  - Must be a valid attribute of `model.rings` (e.g. `:I`, `:v`, `:r`, `:e`, `:i`, `:ϕ`) or a function that can be applied to `model.rings`
  - Example: Keplerian disk time delays could be calculated like `t(ring) = ring.r*(1 .+ sin.(ring.ϕ).*ring.i))`
- `flatten::Bool=false`: If true, flatten the result to a vector

# Returns
- `Array{Float64,}`: Matrix or vector of extracted variable from `model.rings`, created by stacking the output variable for each ring
  - For example, if variable given is `:I`, result will have shape `(length(r), length(ϕ))` as at each `r` and `ϕ` there is a value of `I`
  - If `flatten=true`, result will be a flattened vector
"""
function getVariable(m::model,variable::String;flatten=false) # method for getting variable if String 
    variable = Symbol(variable)
    return getVariable(m,variable;flatten)
end

"""
    getVariable(m::model, variable::Symbol; flatten=false) 

Retrieve model variable when specified as a `Symbol`. See main docstring for details.
"""
function getVariable(m::model,variable::Symbol;flatten=false) # method for getting variable if Symbol
    if variable ∉ fieldnames(ring)
        throw(ArgumentError("variable must be a valid attribute of model.rings\nvalid attributes: $(fieldnames(ring))"))
    end
    if !isnothing(m.cache)
        cached = get(m.cache, (variable, flatten), nothing)
        if !isnothing(cached)
            return cached #treat as read-only -- see cache contract in the model docstring
        end
    end
    isCombined = length(m.subModelStartInds) > 1 #check if model is combined
    res = nothing
    if isCombined
        startInds = m.subModelStartInds
        chunks = []
        l = 0
        for i=1:length(startInds)
            s = startInds[i]; e = i == length(startInds) ? length(m.rings) : startInds[i+1]-1
            chunk = [getfield(ring,variable) for ring in m.rings[s:e]]
            push!(chunks,chunk)
            l+=sum(length,chunk)
        end
        res = Array{Float64}(undef,l) #preallocate array
        for (i,chunk) in enumerate(chunks)
            startInd = i == 1 ? 1 : sum(sum(length,chunks[ii]) for ii in 1:(i-1))+1
            endInd = i == length(chunks) ? l : sum(sum(length,chunks[ii]) for ii in 1:(i))
            if typeof(chunk) == Vector{Vector{Float64}} #if chunk is 2D
                res[startInd:endInd] = vec(stack(chunk,dims=1))
            else
                res[startInd:endInd] = chunk
            end
        end
    else
        res = gatherStack(m.rings, variable)
        if flatten
            res = vec(res) #flatten the matrix to a vector
        end
    end
    if !isnothing(m.cache)
        m.cache[(variable, flatten)] = res
    end
    return res
end

"""
    gatherStack(rings::Vector{ring}, variable::Symbol)

Typed helper (function barrier) for `getVariable`: stacks `getfield(ring, variable)` over all rings
into a `Vector` (one scalar per ring, e.g. cloud models) or a `Matrix` with one row per ring
(vector-valued fields), without the intermediate vector-of-vectors that `stack` requires.
Preserves the element type of the field (e.g. `Bool` for `:reflect`).
"""
function gatherStack(rings::Vector{ring}, variable::Symbol)
    fv = getfield(rings[1], variable)
    if fv isa AbstractVector
        res = Matrix{eltype(fv)}(undef, length(rings), length(fv))
        @inbounds for k in eachindex(rings)
            res[k,:] .= getfield(rings[k], variable)
        end
        return res
    else
        res = Vector{typeof(fv)}(undef, length(rings))
        @inbounds for k in eachindex(rings)
            res[k] = getfield(rings[k], variable)
        end
        return res
    end
end

"""
    getVariable(m::model, variable::Function; flatten=false)

Retrieve model variable when specified as a `Function`. See main docstring for details.

# Note
Results are memoized in `m.cache` only for the package's own pure delay functions
(`t`, `tDisk`, `tCloud`) — user-supplied functions are always re-evaluated, since they may
capture state the cache cannot see.
"""
function getVariable(m::model,variable::Function;flatten=false) # method for getting variable if Function
    if variable == t && length(m.subModelStartInds) == 1 && length(m.rings) > 0
        variable = m.rings[1].θₒ == 0.0 ? tDisk : tCloud #resolve t up front so the cache key is the concrete function
    end
    cacheable = variable === tDisk || variable === tCloud #only package-owned pure functions -- never cache user closures
    if cacheable && !isnothing(m.cache)
        cached = get(m.cache, (variable, flatten), nothing)
        if !isnothing(cached)
            return cached #treat as read-only -- see cache contract in the model docstring
        end
    end
    res = getVariableFunction(m, variable; flatten=flatten)
    if cacheable && !isnothing(m.cache)
        m.cache[(variable, flatten)] = res
    end
    return res
end

function getVariableFunction(m::model,variable::Function;flatten=false) #uncached implementation, see getVariable
    res = nothing
    isCombined = m.subModelStartInds != [1]
    if isCombined
        try
            startInds = m.subModelStartInds
            chunks = []
            l = 0
            for i=1:length(startInds)
                s = startInds[i]; e = i == length(startInds) ? length(m.rings) : startInds[i+1]-1
                if variable == t 
                    variable = m.rings[1].θₒ == 0.0 ? tDisk : tCloud #if variable is t, use tDisk or tCloud depending on whether θₒ is 0
                end
                chunk = [variable(ring) for ring in m.rings[s:e]]
                push!(chunks,chunk)
                l+=sum(length,chunk)
            end
            res = Array{Float64}(undef,l) #preallocate array
            for (i,chunk) in enumerate(chunks)
                startInd = i == 1 ? 1 : sum(length,chunks[i-1])*(i-1)+1
                endInd = i == length(chunks) ? l : sum(length,chunk)*i
                if typeof(chunk) == Vector{Vector{Float64}} #if chunk is 2D
                    res[startInd:endInd] = vec(stack(chunk,dims=1))
                else
                    res[startInd:endInd] = chunk
                end
            end
            return res
        catch
            throw(error("error in function call $(variable)(ring)"))
        end
    else
        if variable == t
            variable = m.rings[1].θₒ == 0.0 ? tDisk : tCloud
        end
        res = stack([variable(ring) for ring in m.rings],dims=1)
        if flatten
            return vec(res) #flatten the matrix to a vector
        else
            return res #return as is
        end
    end
end
"""
    get_rMinMaxDiskWind(r̄::Float64, rFac::Float64, α::Float64)

Calculate the minimum and maximum radius of model given the (intensity weighted) mean radius r̄, 
the radius factor rFac, and the power-law index α following Long+ 2023.

# Parameters
- `r̄::Float64`: mean radius of model (in terms of rₛ)
- `rFac::Float64`: radius scaling factor 
- `α::Float64`: power-law index of the source function ``S(r) \\propto r^{-\\alpha}`` (cannot be 1/2 or 3/2 as this divides by zero)

# Returns
- `rMin::Float64`: minimum radius of model (in terms of rₛ)
- `rMax::Float64`: maximum radius of model (in terms of rₛ)
"""
function get_rMinMaxDiskWind(r̄::Float64,rFac::Float64,α::Float64) 
    rMin = r̄*(3-2*α)/(1-2*α)*(rFac^(1/2-α)-1)/(rFac^(3/2-α)-1)
    rMax = rMin*rFac
    return rMin, rMax
end

"""
    BLR.image(m::model, variable::Union{String,Symbol,Function}, kwargs...)

Generate an image of the model where the color of each point is determined by the variable provided.

# Arguments
- `m::model`: Model object to extract variable from
- `variable::Union{String,Symbol,Function}`: Variable to extract from model
  - If `String`, will be converted to `Symbol`
  - Must be a valid attribute of `model.rings` (e.g. `:I`, `:v`, `:r`, `:e`, `:i`, `:ϕ`) or a function that can be applied to `model.rings`
  - Example: Keplerian disk time delays could be calculated like `t(ring) = ring.r*(1 .+ sin.(ring.ϕ).*ring.i))`

# Keywords
- Additional keyword arguments are passed to `Plots.plot`

# Returns
- `p::Plots.plot`: Plot object representing the generated image
"""
@userplot Image #note that then to call this method use lowercase, i.e. image(m,"I") -- PROBLEM: this doesn't actually loop through each x and y point -- need to collapse them into 1D arrays? 
@recipe function f(img::Image)
    model, variable = nothing, nothing
    if length(img.args) == 2
        model, variable = img.args
    elseif length(img.args) == 1
        model = img.args
        variable = :I #default variable is intensity
    else
        throw(ArgumentError("expected up to 2 arguments (model, [variable]), got $(img.args)"))
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

"""
    addGrid!(m::model, colors=nothing, nϕ::Int=64)

Add a grid to the model image plot - mostly a debugging tool to visualize grid cells of overlapping models.

# Arguments
- `m::model`: Model object to add grid to
- `colors=nothing`: Vector of colors for each submodel (if `nothing`, uses default colors)
- `nϕ::Int=64`: Number of azimuthal angles to use for the grid

# Returns
- `::Plots.plot`: Plot with grid added
"""
function addGrid!(m,colors=nothing,nϕ=64)
    ϕGrid = range(0,stop=2π,length=nϕ)
    rAll = sqrt.(m.camera.α.^2 .+ m.camera.β.^2) #get radii from camera α and β
    ϕAll = atan(m.camera.β,m.camera.α) #to-do: draw grid cells based on phi as well not just r 
    rUnique = unique(rAll) #get unique radii
    if isnothing(colors)
        colors = [i for i in 1:length(m.subModelStartInds)]
    end
    for (ri,r) in enumerate(rUnique)
        i = findfirst(r .== rAll)
        ring,col = getRingFromFlattenedInd(m,i) #get the ring index from the flattened index
        subModelInd = findfirst(ring .< m.subModelStartInds)
        if isnothing(subModelInd)
            subModelInd = length(m.subModelStartInds) #if not found, then it is the last submodel
        else
            subModelInd -= 1 #findfirst returns the first index where the condition is true, so we need to subtract 1 to get the correct subindex
        end
        ΔrUp = m.rings[ring].scale == :log ? r*(exp(m.rings[ring].Δr)-1) : m.rings[ring].Δr #log scale
        ΔrDown = m.rings[ring].scale == :log ? minimum([r,r*(1-1/exp(m.rings[ring].Δr))]) : minimum([r,m.rings[ring].Δr]) #if first ring, then use scale to calculate what "should be" the next inner ring and take minimum of that and the current radius (so don't go through zero)
        # Δϕ = m.rings[ring].Δϕ
        # ϕ = 
        #need to do the "ring -1" check for each subModelStartInd -- if you are the first ring of submodel don't go to the last ring of next submodel to calculate this 
        #also analytically do it based on Δr what the inner ring would have been and then do minimum that and 0 
        lx = (r-ΔrDown/2).* cos.(ϕGrid); ly = (r-ΔrDown/2).* sin.(ϕGrid)
        cellx = vcat(lx,[lx[end],lx[end]])
        ux = (r+ΔrUp/2).* cos.(ϕGrid); uy = (r+ΔrUp/2).* sin.(ϕGrid)
        plot!(vcat(lx,ux),vcat(ly,uy),label="",color=colors[subModelInd],fill=true,fillalpha=0.1,linestyle=:dash,lw=2)
        plot!(r.*cos.(ϕGrid),r.*sin.(ϕGrid),label="",color=colors[subModelInd],linestyle=:solid,lw=2)
    end
end


"""
    get_r3D(i::Float64, rot::Float64, θₒ::Float64) -> Matrix{Float64}

Calculate rotation matrix to transform from initial XY plane coordinates to 3D space.

# Parameters
- `i::Float64`: Inclination angle of ring (rad)
- `rot::Float64`: Rotation of ring plane about z axis (rad)
- `θₒ::Float64`: Opening angle of point (rad)

# Returns
- `matrix::Matrix{Float64}`: 3×3 rotation matrix
"""
function get_r3D(i,rot,θₒ)
    matrix = [
        (cos(θₒ)*cos(rot)*sin(i) - sin(θₒ)*cos(i)) -sin(rot)*sin(i) (sin(θₒ)*cos(rot)*sin(i) + cos(θₒ)*cos(i));
        cos(θₒ)*sin(rot) cos(rot) sin(θₒ)*sin(rot);
        -(cos(θₒ)*(cos(rot))*cos(i) + sin(θₒ)*sin(i)) sin(rot)*cos(i) (cos(θₒ)*sin(i) - sin(θₒ)*cos(rot)*cos(i))
    ]
    return matrix
end
"""
    reflect!(xyzSys::Vector{Float64}, i::Float64) -> Vector{Float64}

Reflect coordinates in 3D space across the ring plane.

# Parameters
- `xyzSys::Vector{Float64}`: `[x;y;z]` coordinates in 3D space
- `i::Float64}`: Inclination angle of ring plane (rad)

# Returns
- `xyzSys::Vector{Float64}`: `[x';y';z']` coordinates in 3D space after reflection
"""
function reflect!(xyzSys,i)
    #reflect across line made with inclination angle z - m*x = 0, where m = -cot(i)
    den = 1 + cot(i)^2
    xf = (xyzSys[1]*(1-cot(i)^2) - 2*xyzSys[3]*cot(i))/den
    zf = (xyzSys[3]*(cot(i)^2-1) - 2*xyzSys[1]*cot(i))/den
    xyzSys[1] = xf
    xyzSys[3] = zf
    return xyzSys
end
"""
    reflect_scalar(x::Float64, y::Float64, z::Float64, i::Float64) -> NTuple{3,Float64}

Allocation-free scalar version of `reflect!` — reflect coordinates in 3D space across the ring plane.
Same math as `reflect!` but takes/returns plain scalars instead of mutating a `Vector`.
"""
function reflect_scalar(x::Float64, y::Float64, z::Float64, i::Float64)
    #reflect across line made with inclination angle z - m*x = 0, where m = -cot(i) (math copied from reflect!)
    cti = cot(i)
    den = 1 + cti^2
    xf = (x*(1-cti^2) - 2*z*cti)/den
    zf = (z*(cti^2-1) - 2*x*cti)/den
    return (xf, y, zf)
end

"""
    rotate3D_vector_scalar(vx::Float64, vy::Float64, vz::Float64, i::Float64, rot::Float64, θₒ::Float64, reflect::Bool=false) -> NTuple{3,Float64}

Allocation-free rotation of an arbitrary vector `(vx,vy,vz)` from initial XY-plane coordinates into 3D space
(camera at +x), with optional reflection across the ring plane. The matrix entries are those of
`get_r3D(i,rot,θₒ)` inlined; used for rotating velocity vectors as well as positions.
"""
function rotate3D_vector_scalar(vx::Float64, vy::Float64, vz::Float64, i::Float64, rot::Float64, θₒ::Float64, reflect::Bool=false)
    sini, cosi = sincos(i)
    sinrot, cosrot = sincos(rot)
    sinθₒ, cosθₒ = sincos(θₒ)
    x = (cosθₒ*cosrot*sini - sinθₒ*cosi)*vx - sinrot*sini*vy + (sinθₒ*cosrot*sini + cosθₒ*cosi)*vz
    y = cosθₒ*sinrot*vx + cosrot*vy + sinθₒ*sinrot*vz
    z = -(cosθₒ*cosrot*cosi + sinθₒ*sini)*vx + sinrot*cosi*vy + (cosθₒ*sini - sinθₒ*cosrot*cosi)*vz
    return reflect ? reflect_scalar(x, y, z, i) : (x, y, z)
end

"""
    rotate3D_scalar(r::Float64, ϕ₀::Float64, i::Float64, rot::Float64, θₒ::Float64, reflect::Bool=false) -> NTuple{3,Float64}

Allocation-free version of `rotate3D`: transform from ring coordinates to 3D coordinates where camera
is at +x, returning an `(x, y, z)` tuple instead of a heap-allocated `Vector`. Prefer this in hot loops.
"""
function rotate3D_scalar(r::Float64, ϕ₀::Float64, i::Float64, rot::Float64, θₒ::Float64, reflect::Bool=false)
    sinϕ₀, cosϕ₀ = sincos(ϕ₀)
    return rotate3D_vector_scalar(r*cosϕ₀, r*sinϕ₀, 0.0, i, rot, θₒ, reflect)
end

"""
    rotate3D(r::Float64, ϕ₀::Float64, i::Float64, rot::Float64, θₒ::Float64, reflect::Bool=false)

Transform from ring coordinates to 3D coordinates where camera is at +x.

# Parameters
- `r::Float64`: Radius from central mass (in terms of rₛ)
- `ϕ₀::Float64`: Starting azimuthal angle in ring plane (rad)
- `i::Float64`: Inclination angle of ring plane (rad)
- `rot::Float64`: Rotation of system plane about z axis (rad)
- `θₒ::Float64`: Opening angle of point (rad)
- `reflect::Bool=false`: Whether to reflect across the ring plane

# Returns
- `Vector{Float64}`: `[x; y; z]` coordinates in 3D space

# Note
This method allocates its result; in performance-critical loops use [`rotate3D_scalar`](@ref) which
returns a stack-allocated tuple instead.
"""
function rotate3D(r::Float64,ϕ₀::Float64,i::Float64,rot::Float64,θₒ::Float64,reflect::Bool=false)
    return collect(rotate3D_scalar(r,ϕ₀,i,rot,θₒ,reflect))
end
"""
    getXYZ(ring::ring)

Return the 3D system coordinates `(x, y, z)` of every point in `ring` (camera at +x).

Uses the cached `ring.x/y/z` fields when present; otherwise computes them with
`rotate3D_scalar` and stores them on the (mutable) ring so the transform is only
ever done once per point.

# Returns
- For a cloud/point ring (scalar fields): `NTuple{3,Float64}`
- For a continuous ring (vector fields): `Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}`

!!! warning "Reflection already applied"
    Cached coordinates are *post-reflection* — if `ring.reflect` is true the stored values already
    include the reflection across the disk midplane. Never apply `reflect!`/`reflect_scalar` to the
    result of `getXYZ` a second time.

!!! warning "Geometry mutation"
    The cache is filled from `r`, `ϕ₀`, `i`, `rot`, `θₒ` at first access. If you mutate any of those
    fields afterwards, set `ring.x = nothing; ring.y = nothing; ring.z = nothing` to force a recompute.
"""
function getXYZ(ring::ring)
    if !isnothing(ring.x) && !isnothing(ring.y) && !isnothing(ring.z)
        return ring.x, ring.y, ring.z
    end
    if typeof(ring.r) == Float64 && typeof(ring.ϕ₀) == Float64 #cloud/point ring
        xyz = rotate3D_scalar(ring.r, ring.ϕ₀, ring.i, ring.rot, ring.θₒ, ring.reflect)
        ring.x = xyz[1]; ring.y = xyz[2]; ring.z = xyz[3]
        return xyz
    else #continuous ring -- per-point fields may be scalars or vectors (see removeNaN! for the canonical pattern)
        n = length(ring.r)
        x = Vector{Float64}(undef,n); y = Vector{Float64}(undef,n); z = Vector{Float64}(undef,n)
        for k in 1:n
            ik = typeof(ring.i) == Float64 ? ring.i : ring.i[k]
            rotk = typeof(ring.rot) == Float64 ? ring.rot : ring.rot[k]
            θₒk = typeof(ring.θₒ) == Float64 ? ring.θₒ : ring.θₒ[k]
            reflectk = typeof(ring.reflect) == Bool ? ring.reflect : ring.reflect[k]
            x[k], y[k], z[k] = rotate3D_scalar(ring.r[k], ring.ϕ₀[k], ik, rotk, θₒk, reflectk)
        end
        ring.x = x; ring.y = y; ring.z = z
        return x, y, z
    end
end

"""
    rotate3D(r::Float64, ϕ₀::Float64, i::Float64, matrix::Matrix{Float64}, reflect::Bool=false)
Transform from ring coordinates to 3D coordinates where camera is at +x using a precomputed rotation matrix.
# Parameters
- `r::Float64`: Radius from central mass (in terms of rₛ)
- `ϕ₀::Float64`: Starting azimuthal angle in ring plane (rad)
- `i::Float64`: Inclination angle of ring plane (rad)
- `matrix::Matrix{Float64}`: Precomputed rotation matrix for the ring
- `reflect::Bool=false`: Whether to reflect across the ring plane
# Returns
- `Vector{Float64}`: `[x; y; z]` coordinates in 3D space
"""
function rotate3D(r::Float64,ϕ₀::Float64,i::Float64,matrix::Matrix{Float64},reflect::Bool=false)
    xyzSys = matrix*[r*cos(ϕ₀);r*sin(ϕ₀);0]
    if reflect
        xyzSys = BLR.reflect!(xyzSys,i) #use the inclination angle from the rotation matrix
    end
    return xyzSys
end


"""
    midPlaneXZ(x::Float64, i::Float64) -> Float64

midPlaneXZ(x,i) = -x*cot(i)
"""
midPlaneXZ(x,i) = -x*cot(i)

function geometry(ring) #for 3d visualzation of model geometry
    return length(ring.ϕ) == 1 ? 1.0 : ones(length(ring.ϕ)) #if no variable is given, just return 1.0 for each point
end
"""
    plot3d(m::model, [variable], [annotate], [kwargs...])

Generate a 3D plot of the model geometry, optionally colored by a variable.

## Parameters
- `m::model`: Model object to plot
- `variable::Union{String,Symbol,Function}=nothing`: Variable to color the points by
  - If `String`, will be converted to `Symbol`
  - Must be a valid attribute of `model.rings` (e.g., `:I`, `:v`, `:r`, `:e`, `:i`, `:ϕ`) or a function that can be applied to `model.rings`
  - Example: Keplerian disk time delays could be calculated like `t(ring) = ring.r*(1 .+ sin.(ring.ϕ).*ring.i))`
  - If not provided, defaults to `nothing` (no coloring)
- `annotate::Bool=true`: Whether to annotate the camera position and model orientation in the plot
- `kwargs...`: Additional keyword arguments passed to `Plots.plot`

## Returns
- `p::Plots.plot`: 3D plot of the model geometry, optionally colored by the variable provided
"""
@userplot Plot3d #note that then to call this method use lowercase, i.e. plot3d(m,"I") -- WIP need to make more general for combined Models -- do in 2 steps if isCombined? just call twice
@recipe function f(p::Plot3d)
    xlabel --> "x [rₛ]"
    ylabel --> "y [rₛ]"
    zlabel --> "z [rₛ]"
    aspect_ratio --> :equal
    model, variable, annotate = nothing, nothing, true
    if length(p.args) == 2
        model, tmp = p.args
        if typeof(tmp) == Bool
            annotate = tmp
        else
            variable = tmp
        end
    elseif length(p.args) == 1
        model = p.args[1]
        variable = BLR.geometry #if no variable is given, just return 1.0 for each point
    elseif length(p.args) == 3
        model, tmp1, tmp2 = p.args
        if typeof(tmp1) == Bool
            annotate = tmp1
            variable = tmp2
        else
            variable = tmp1
            annotate = tmp2
        end
    else
        throw(ArgumentError("expected 1, 2, or 3 arguments, got $(length(p.args))"))
    end
    variable = isa(variable,Function) ? variable : Symbol(variable)
    isCombined = length(model.subModelStartInds) > 1 #check if model is combined
    startInds = model.subModelStartInds
    diskFlags = [typeof(model.rings[i].I) != Float64 for i in startInds] 
    mList = [deepcopy(model[i]) for i=1:length(startInds)]
    title --> (variable == geometry ? "System geometry visualization" : "System geometry + $variable (color) visualization")
    boxSizeGlobal = 0.0
    for (mInd,model) in enumerate(mList)
        i = getVariable(model,:i)
        r = getVariable(model,:r)
        ϕ₀ = getVariable(model,:ϕ₀)
        xtmp = zeros(size(ϕ₀))
        ytmp = zeros(size(ϕ₀))
        ztmp = zeros(size(ϕ₀))
        if typeof(r) == Matrix{Float64} && typeof(ϕ₀) == Matrix{Float64}
            for ii in 1:size(r)[1]
                for jj in 1:size(r)[2]
                    rot = model.rings[ii].rot
                    xtmp[ii,jj],ytmp[ii,jj],ztmp[ii,jj] = rotate3D_scalar(r[ii,jj],ϕ₀[ii,jj],i[ii],rot,model.rings[ii].θₒ,model.rings[ii].reflect)
                end
            end
        elseif typeof(r) == Vector{Float64} && typeof(ϕ₀) == Matrix{Float64}
            for ii in 1:size(ϕ)[1]
                for jj in 1:size(ϕ)[2]
                    rot = model.rings[ii].rot
                    xtmp[ii,jj],ytmp[ii,jj],ztmp[ii,jj] = rotate3D_scalar(r[ii],ϕ₀[ii,jj],i[ii],rot,model.rings[ii].θₒ,model.rings[ii].reflect)
                end
            end
        else #if r is just a vector (with ϕ and i matching) -- cloud/point rings, one point per ring
            for (ii,r) in enumerate(model.rings)
                xtmp[ii],ytmp[ii],ztmp[ii] = getXYZ(r) #cached system coordinates
            end
        end
        boxSize = 1.1*maximum([maximum(i for i in xtmp if !isnan(i)),maximum(i for i in ytmp if !isnan(i)),maximum(i for i in ztmp if !isnan(i))])
        boxSizeGlobal = max(boxSize,boxSizeGlobal)
        @series begin
            subplot := 1
            seriestype := :scatter
            palette --> :magma
            mz = vec(getVariable(model,variable))
            nanMask = isnan.(mz)
            marker_z := mz[.!nanMask]
            markeralpha --> (diskFlags[mInd] ? 0.9 : 0.1)
            markerstrokewidth --> 0.0
            markersize --> 1.
            label --> ""
            x := vec(xtmp)[.!nanMask]
            y := vec(ytmp)[.!nanMask]
            z := vec(ztmp)[.!nanMask]
            ()
        end
        if annotate
            @series begin
                subplot := 1
                x := [-boxSize,boxSize]./1.1
                y := [0.0,0.0]
                z := [midPlaneXZ(-boxSize/1.1,i[1]),midPlaneXZ(boxSize/1.1,i[1])]
                seriestype := :path
                cList = [:crimson,:darkorange]
                color --> (mInd <= 2 ? cList[mInd] : mInd)
                label --> (length(mList) == 1 ? "midplane" : "midplane $mInd ($(diskFlags[mInd] ? "disk" : "cloud"))")
                ()
            end
        end
    end
    boxSize = boxSizeGlobal
    ylims --> (-boxSize,boxSize)
    xlims --> (-boxSize,boxSize)
    zlims --> (-boxSize,boxSize)
    foreground_color_legend --> nothing
    colorbar --> variable == BLR.geometry
    if annotate
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
    end
end

"""
    profile(m::model, [variable], [kwargs...])

Plot all profiles set in the model, normalized to the maximum value of each profile.

# Arguments
- `model`: A model object containing profile data.
- `variable`: Optional. A symbol or string (or list of symbols/strings) specifying which profile to plot. If not provided, all profiles set in model will be plotted.
- `kwargs...`: Additional keyword arguments passed to `Plots.plot`
"""
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
        norm = length(variable) == 1 ? 1.0 : maximum(abs(i) for i in m.profiles[v].binSums if !isnan(i))
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