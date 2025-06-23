#!/usr/bin/env julia
"""
    binnedSum(x::Array{Float64,}, y::Array{Float64, }; bins=100, 
            overflow=false, centered=false, minX=nothing, maxX=nothing)

Bin the x and y variables into a histogram, where each bin is the sum of the y values for the corresponding x values.

# Arguments
- `x::Array{Float64,}`: x variable to bin over
- `y::Array{Float64,}`: y variable to bin
- `bins::Union{Int,Vector{Float64}}=100`: Number of bins or bin edges for binning
  - If `Int`: Number of bins with edges equally spaced between min/max of x
  - If `Vector{Float64}`: Specific bin edges, with number of bins = `length(bins)-1`
  - Left edge inclusive, right edge exclusive (except last bin which is inclusive)
- `overflow::Bool=false`: If `true`, include values outside bin range in the first/last bins
- `centered::Bool=false`: If `true`, shift bin edges to center around middle value
- `minX::Union{Float64,Nothing}=nothing`: Minimum value of x for binning (defaults to `minimum(x)`)
- `maxX::Union{Float64,Nothing}=nothing`: Maximum value of x for binning (defaults to `maximum(x)`)

# Returns
- `Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}`: A tuple containing:
  - `binEdges`: Bin edges for the x variable
  - `binCenters`: Bin centers for the x variable
  - `result`: Binned sums of the y variable
"""
function binnedSum(x::Array{Float64,}, y::Array{Float64, }; bins::Union{Int,Vector{Float64}}=100, overflow::Bool = false, centered::Bool=true, minX::Union{Float64,Nothing}=nothing, maxX::Union{Float64,Nothing}=nothing)
    binEdges = nothing 
    if typeof(bins) == Int
        if isnothing(minX)
            minX = isnan(minimum(x)) ? minimum(i for i in x if !isnan(i)) : minimum(x)
        end
        if isnothing(maxX)
            maxX = isnan(maximum(x)) ? maximum(i for i in x if !isnan(i)) : maximum(x)
        end
        if centered
            middle = (maxX+minX)/2
            Δ = (maxX-minX)/bins
            nBins = 2*ceil(Int,(maxX-(middle+Δ/2))/Δ) + 1
            maxX = (middle + Δ/2) + nBins/2*Δ
            minX = (middle - Δ/2) - nBins/2*Δ
            binEdges = collect(range(minX,stop=maxX,length=nBins+1))
        else
            binEdges = collect(range(minX,stop=maxX,length=bins+1))
        end
    else
        binEdges = bins
    end
    binCenters = @. (binEdges[1:end-1] + binEdges[2:end])/2
    result = zeros(length(binEdges)-1)
    for (x,y) in zip(x,y)
        if isfinite(x) && isfinite(y)
            if x <= binEdges[1]
                if overflow
                    result[1] += y
                end
            elseif x >= binEdges[end]
                if overflow
                    result[end] += y
                end
            else
                bin = searchsortedfirst(binEdges,x)-1
                result[bin] += y
            end
        end
    end
    return (binEdges, binCenters, result)
end

"""
    binModel(bins::Union{Int,Vector{Float64}}=100; m::model, yVariable::Union{String,Symbol,Function}, 
            xVariable::Union{String,Symbol,Function}=:v, kwargs...)
    binModel(bins::Vector{Float64}, dx::Array{Float64,}; m::model, yVariable::Union{String,Symbol,Function}, 
            xVariable::Union{String,Symbol,Function}=:v, kwargs...)

Bin the model into a histogram, where each bin is the integrated value of the yVariable as a function of the xVariable.

# Arguments
- `m::model`: Model object to bin
- `yVariable::Union{String,Symbol,Function}`: Variable to bin
  - Must be valid attribute of `model.rings` (e.g. `:I`, `:v`, `:r`, `:e`, `:i`, `:ϕ`) or a function that can be applied to `model.rings`
  - Example: Keplerian disk time delays could be calculated like `t(ring) = ring.r*(1 .+ sin.(ring.ϕ).*ring.i))`
- `bins::Union{Int,Vector{Float64}}`: Number of bins or bin edges for binning
  - If `Int`: Number of bins with edges equally spaced between min/max of xVariable
  - If `Vector{Float64}`: Specific bin edges, with number of bins = `length(bins)-1`
  - Left edge inclusive, right edge exclusive (except last bin which is inclusive)
- `xVariable::Union{String,Symbol,Function}=:v`: Variable to bin over
  - Must be a valid attribute of `model.rings` or a function that can be applied to `model.rings`
- `dx::Array{Float64,}`: Integration element for each bin
  - If provided, used as the associated integral element
  - Otherwise defaults to `ΔA` in each ring struct

# Returns
- `Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}`: A tuple containing:
  - `binEdges`: Bin edges for the xVariable of the histogram
  - `binCenters`: Bin centers for the xVariable
  - `yBinned`: Binned values of the yVariables
"""
function binModel(bins::Union{Int,Vector{Float64}}=100; m::model, yVariable::Union{String,Symbol,Function}, xVariable::Union{String,Symbol,Function}=:v, kwargs...)
    x = getVariable(m,xVariable)
    y = getVariable(m,yVariable) 
    dx = getVariable(m,:ΔA)
    if size(y) != size(dx)
        if size(dx)[1] == size(y)[1] && typeof(m.rings[1].ΔA) == Float64
            dx = stack([dx for ϕ in 1:size(y)[2]],dims=2)
        else
            throw(ArgumentError("y and dA must be the same size, got $(size(y)) and $(size(dx))\nconsider supplying custom dx = Array{Float64,}(size(y)) (use ones for equal weighting)"))
        end
    end

    return binnedSum(x,y.*dx,bins=bins;kwargs...)
end

function binModel(bins::Vector{Float64}, dx::Array{Float64,}; m::model, yVariable::Union{String,Symbol,Function}, xVariable::Union{String,Symbol,Function}=:v, kwargs...)
    x = getVariable(m,xVariable)
    y = getVariable(m,yVariable) 
    if size(y) != size(dx)
        throw(ArgumentError("y and dx must be the same size, got $(size(y)) and $(size(dx))"))
    end

    return binnedSum(x,y.*dx,bins=bins;kwargs...)
end

"""
    tDisk(ring::ring)

    Calculate time delays for a point in a disk as `` t = \\eta r \\left(1 + \\cos(\\phi) \\sin(i)\\right)``
"""
tDisk(ring::ring) = @. ring.η*ring.r*(1 + cos(ring.ϕ)*sin(ring.i)) #time delays for Keplerian disk [rₛ], or a cloud modelled as a point in a disk

"""
    tCloud(ring::ring)

    Calculate time delays for a cloud with opening angle ``\\theta_o`` as the x-coordinate of the point subtracted from the radial distance of the point ``t = r - x``.
"""
tCloud(ring::ring) = begin
    xyzSys = rotate3D(ring.r,ring.ϕ₀,ring.i,ring.rot,ring.θₒ,ring.reflect) #system coordinates xyz
    return ring.η*(ring.r - xyzSys[1]) #could also calculate new incliation angle based on θₒ, but this is simpler, +x
end
"""
    t(ring::ring)

    Calculate time delays for a point in a disk as `` t = \\eta r \\left(1 + \\cos(\\phi) \\sin(i)\\right)`` or a cloud with opening angle ``\\theta_o``
    as the x-coordinate of the point subtracted from the radial distance of the point ``t = r - x``.
"""
t(ring::ring) = ring.θₒ == 0.0 ? tDisk(ring) : tCloud(ring) #if θₒ = 0.0, then this is a point in a disk, otherwise it is a cloud with opening angle θₒ
"""
    t(ring::ring, subFxn::Function=tDisk)

    Calculate time delays for a point in a disk or cloud using a custom function `subFxn` that takes a `ring` struct.
    It is more peformant to pass the function directly rather than figure it out on the fly if known ahead of time.
"""
t(ring::ring,subFxn::Function=tDisk) = subFxn(ring) #allow for custom time delay function, e.g. tDisk or tCloud

"""
    getG(β::Float64)
"""
    getG(β::Float64)

"""
    phase(m::model; U, V, PA, BLRAng, returnAvg=false, offAxisInds=nothing, kwargs...)

Calculate differential phase signal for a model based on specified baselines, model orientation, and BLR angular size.

# Arguments
- `m::model`: Model object to calculate phase for
- `U::Vector{Float64}`: U component of complex visibility in [Mλ]
- `V::Vector{Float64}`: V component of complex visibility in [Mλ]
- `PA::Float64`: On-sky position angle of the model in radians
- `BLRAng::Float64`: Characteristic size of the BLR model in radians (conversion from ``r_s`` to radians)
- `returnAvg::Bool=false`: If `true`, returns the average phase across all baselines
- `offAxisInds::Union{Nothing,Vector{Int}}=nothing`: If provided, only calculates phase for baselines at specified indices

# Returns
- If `returnAvg=true`: `Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}` containing:
  - Bin edges for velocity
  - Bin centers for velocity
  - Average differential phase (in radians)
- If `returnAvg=false`: `Vector{Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}` containing:
  - For each baseline, a tuple of bin edges, bin centers, and differential phase

The differential phase is calculated by integrating the phase over the model at each velocity bin,
weighted by the intensity and area of each ring element.
"""
function phase(m::model; returnAvg::Bool=false, offAxisInds::Union{Nothing,Vector{Int}}=nothing, 
    U::Vector{Float64}, V::Vector{Float64}, PA::Float64, BLRAng::Float64, kwargs...)    
    if size(U) != size(V)
        throw(ArgumentError("U and V must be the same size, got $(size(U)) and $(size(V))"))
    end
    X = m.camera.α.*BLRAng; Y = m.camera.β.*BLRAng
    function getΔϕ(v::Array{Float64,},I::Array{Float64,},ΔA::Array{Float64,},x::Array{Float64,},y::Array{Float64},U::Float64,V::Float64,PA::Float64;kwargs...)
        U′ = -cos(PA)*U-sin(PA)*V; V′ = sin(PA)*U-cos(PA)*V
        edges,centers,Δϕ = binnedSum(v, -2π*(x.*U′ .+ y.*V′).*I.*ΔA.*(180/π*1e6);kwargs...) 
        edges,centers,LP = binnedSum(v, I.*ΔA;kwargs...)
        return (edges,centers,Δϕ./(1.0.+LP)) #divide out LP but then weight by f/(1+f) to get the differential phase
    end
    v = getVariable(m,:v,flatten=true); I = getVariable(m,:I,flatten=true); ΔA = getVariable(m,:ΔA,flatten=true)
    if isnothing(offAxisInds)
        offAxisInds = collect(1:length(U))
    end
    U = U[offAxisInds]; V = V[offAxisInds]
    ΔϕList = [getΔϕ(v,I,ΔA,X,Y,Ui,Vi,PA;kwargs...) for (Ui,Vi) in zip(U,V)]
    if returnAvg
        edges, centers, _ = ΔϕList[1] #get the edges and centers from the first Δϕ
        ΔϕList = [ϕ[3] for ϕ in ΔϕList]
        ΔϕAvg = mean(ΔϕList,dims=1)[1]
        return (edges,centers,ΔϕAvg)
    else
        return ΔϕList
    end
end

"""
    getProfile(m::model, name::Union{String,Symbol,Function}; 
               bins::Union{Int,Vector{Float64}}=100, 
               dx::Union{Array{Float64,},Nothing}=nothing, kwargs...)

Return a profile for the model based on the specified name.

# Arguments
- `m`: Model object to get the profile from
- `name`: Name of the profile to get. Options include:
  - `:line`: Returns the line profile (integrated intensity as function of velocity)
  - `:delay`: Returns the delay profile (mean delays weighted by intensity/responsivity vs. velocity)
  - `:r`: Returns the mean radius (weighted by intensity) as function of velocity 
  - `:ϕ`: Returns the mean azimuthal angle (weighted by intensity) as function of velocity
  - `:phase`: Returns the phase profile (integrated phase as function of velocity)
    - Requires `U` [Mλ], `V` [Mλ], `PA` [rad], and `BLRAng` [rad] as keyword arguments
  - `Function`: Returns the intensity weighted mean of this function vs. velocity
- `bins`: Number of bins or bin edges for binning
  - If `Int`: Number of bins with edges equally spaced between min/max velocity
  - If `Vector{Float64}`: Specific bin edges, with number of bins = `length(bins)-1`
  - Left edge inclusive, right edge exclusive (except last bin which is inclusive)
- `dx`: Integration element for each ring (defaults to `ΔA` in each ring struct if `nothing`)
- Additional `kwargs` passed to `binModel` include:
  - `overflow=true`: Include overflow bins
  - `centered=true`: Center bins around 0.0
  - `minX`, `maxX`: Set min/max bin boundaries

# Returns
- `profile`: A profile object containing bin edges, bin centers, and binned sums.
  Assign to a model object using `setProfile!` to store the profile in the model.
"""
function getProfile(m::model, name::Union{String,Symbol,Function}; bins::Union{Int,Vector{Float64}}=100, dx::Union{Array{Float64,},Nothing}=nothing, kwargs...)
    n = Symbol(name); p = nothing;
    if n == :line
        p = isnothing(dx) ? binModel(bins,m=m,yVariable=:I,xVariable=:v;kwargs...) : binModel(bins,dx,m=m,yVariable=:I,xVariable=:v;kwargs...)
    elseif n == :delay
        #d(ring::ring;kwargs...) = (typeof(ring.r) == Float64 && typeof(ring.ϕ) == Float64) ? tCloud(ring;kwargs...).*ring.I : tDisk(ring;kwargs...).*ring.I
        d(ring::ring;kwargs...) = t(ring;kwargs...).*ring.I 
        den(ring::ring;) = (typeof(ring.r) == Float64 && typeof(ring.ϕ) == Float64) ? ring.I*ring.η : ring.I.*ring.η 
        pNum = isnothing(dx) ? binModel(bins,m=m,yVariable=d,xVariable=:v;kwargs...) : binModel(bins,dx,m=m,yVariable=d,xVariable=:v;kwargs...)
        pDen = isnothing(dx) ? binModel(bins,m=m,yVariable=den,xVariable=:v;kwargs...) : binModel(bins,dx,m=m,yVariable=den,xVariable=:v;kwargs...)
        p = (pNum[1], pNum[2], pNum[3]./pDen[3])
    elseif n == :r
        r(ring::ring) = ring.r.*ring.I
        pNum = isnothing(dx) ? binModel(bins,m=m,yVariable=r,xVariable=:v;kwargs...) : binModel(bins,dx,m=m,yVariable=:r,xVariable=:v;kwargs...)
        pDen = isnothing(dx) ? binModel(bins,m=m,yVariable=:I,xVariable=:v;kwargs...) : binModel(bins,dx,m=m,yVariable=:I,xVariable=:v;kwargs...)
        p = (pNum[1], pNum[2], pNum[3]./pDen[3])
    elseif n == :ϕ
        ϕ(ring::ring) = ring.ϕ.*ring.I
        pNum = isnothing(dx) ? binModel(bins,m=m,yVariable=ϕ,xVariable=:v;kwargs...) : binModel(bins,dx,m=m,yVariable=:ϕ,xVariable=:v;kwargs...)
        pDen = isnothing(dx) ? binModel(bins,m=m,yVariable=:I,xVariable=:v;kwargs...) : binModel(bins,dx,m=m,yVariable=:I,xVariable=:v;kwargs...)
        p = (pNum[1], pNum[2], pNum[3]./pDen[3])
    elseif n == :phase
        edges,centers,avgPhase = phase(m,returnAvg=true;kwargs...)
        p = (edges, centers, avgPhase)
    elseif isa(name,Function)
        pNum = isnothing(dx) ? binModel(bins,m=m,yVariable=name,xVariable=:v;kwargs...) : binModel(bins,dx,m=m,yVariable=name,xVariable=:v;kwargs...)
        pDen = isnothing(dx) ? binModel(bins,m=m,yVariable=:I,xVariable=:v;kwargs...) : binModel(bins,dx,m=m,yVariable=:I,xVariable=:v;kwargs...)
        p = (pNum[1], pNum[2], pNum[3]./pDen[3])
    else
        throw(ArgumentError("profile $(name) not recognized -- choose from [:line, :delay, :r, :ϕ] or pass a function that can be applied to model.rings"))
    end
    return profile(name=n,binEdges=p[1],binCenters=p[2],binSums=p[3])
end

"""
    setProfile!(m::model, p::profile; overwrite::Bool=false)

Set a profile in the model's profiles dictionary. If the profile already exists and overwrite is false, a warning is issued.
"""
function setProfile!(m::model, p::profile; overwrite::Bool=false) #mutates the model.profiles dictionary based on passed on name passed (i.e. LP, delayProfile, etc.). can also take a function + kwargs to set a custom profile. stores bins + yvalues in each dict entry.
    if isdefined(m,:profiles) 
        if p.name ∈ keys(m.profiles) && !overwrite
            @warn("profile already exists in model.profiles -- to overwrite set overwrite=true")
        else
            m.profiles[p.name] = p
        end
    else
        m.profiles = Dict(p.name => p)
    end
    return m
end

