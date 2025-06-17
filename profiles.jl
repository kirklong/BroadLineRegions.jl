#!/usr/bin/env julia
function binnedSum(x::Array{Float64,}, y::Array{Float64, }; bins::Union{Int,Vector{Float64}}=100, overflow::Bool = false, centered::Bool=false, minX::Union{Float64,Nothing}=nothing, maxX::Union{Float64,Nothing}=nothing)
    """bin the x and y variables into a histogram, where each bin is the sum of the y values for the corresponding x values
    params:
        x: Array{Float64,}
            - x variable to bin over
        y: Array{Float64,}
            - y variable to bin
        bins: Union{Int,Vector{Float64}} = 100
            - number of bins or bin edges for binning
            - if {Int} then bins is the number of bins, and the edges of the bins will be equally spaced between the min and max of the x variable
            - if {Vector{Float64} then bins are the bin edges and the nubmer of bins will be length(bins)-1
                - left edge of each bin is inclusive, right edge is exclusive except for the last bin which is inclusive
        overflow: Bool = false
            - if true, will include values less than the minimum bin edge in the minimum bin and values greater than the maximum bin edge in the maximum bin
        centered: Bool = false
            - if true, shift the bin edges such that they are centered around middle value
        minX: Union{Float64,Nothing} = nothing
            - minimum value of x variable to use for binning
            - if nothing, will use the minimum value of x
        maxX: Union{Float64,Nothing} = nothing
            - maximum value of x variable to use for binning
            - if nothing, will use the maximum value of x
    returns:
        (binEdges, binCenters, result): Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}
            - bin edges and centers for the x variable of the histogram and the binned sums of the y variable
    """
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

"""bin the model into a histogram, where each bin is the integrated value of the yVariable as a function of the xVariable
params: 
    m: model
        - model object to bin
    yVariable: Union{String, Symbol, Function} -- if String then will be converted to symbols
        - y variable to bin
        - must be valid attribute of model.rings (e.g. :I, :v, :r, :e, :i, :ϕ) or a function that can be applied to model.rings
            - example: Keplerian disk time delays could be calculated like `t(ring) = ring.r*(1 .+ sin.(ring.ϕ).*ring.i))`
    bins: Union{Int, Vector{Float64}}
        - number of bins or bin edges for binning
        - if number of bins, then the bins will be equally spaced between the min and max of the xVariable
        - if bin edges then the number of bins will be length(bins)-1
            - leftmost edge is then the first bin and the rightmost edge being the last bin
        - left edge of each bin is inclusive, right edge is exclusive except for the last bin which is inclusive
    xVariable: Union{String, Symbol, Function}
        - variable to bin over
        - must be a valid attribute of model.rings (e.g. :I, :v, :r, :e, :i, :ϕ) or a function that can be applied to model.rings
    dx: each bin is an integral over the range of xVariable -- if dx is provided this is used as the associated integral element, otherwise defaults to ΔA in each ring struct
returns:
    (binEdges, binCenters): Tuple{Vector{Float64},Vector{Float64}}
        - bin edges and centers for the xVariable of the histogram 
    yBinned: Vector{Vector{Float64}}
        - binned values of the yVariables
        - length of yBinned is the same as the number of yVariables, with each element being a vector of length equal to the number of bins
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
end,
function binModel(bins::Vector{Float64}, dx::Array{Float64,}; m::model, yVariable::Union{String,Symbol,Function}, xVariable::Union{String,Symbol,Function}=:v, kwargs...)
    x = getVariable(m,xVariable)
    y = getVariable(m,yVariable) 
    if size(y) != size(dx)
        throw(ArgumentError("y and dx must be the same size, got $(size(y)) and $(size(dx))"))
    end

    return binnedSum(x,y.*dx,bins=bins;kwargs...)
end

"""
time delays for a point in a disk 
params:
    ring: ring
        - ring object to calculate time delays for
returns:
    t: Array{Float64,}
        - time delays for each ring in the model
"""
t(ring::ring) = ring.η.*ring.r.*(1 .+ cos.(ring.ϕ).*sin(ring.i)) # time delays for Keplerian disk [rₛ], or a cloud modelled as a point in a disk

function phase(m::model; returnAvg::Bool=false, offAxisInds::Union{Nothing,Vector{Int}}=nothing, 
    U::Vector{Float64}, V::Vector{Float64}, PA::Float64, BLRAng::Float64, kwargs...)    
    """
    Calculate diferential phase for a model ring
    params:
        m: {model}
            - model object to calculate phase for} 
        U: Vector{Float64}
            - U component of the complex visibility, units of [Mλ]
        V: Vector{Float64}
            - V component of the complex visibility, units of [Mλ]
        PA: Float64
            - on-sky position angle of the model in radians
        BLRAng: Float64
            - characteristic size of the BLR model in radians (i.e. conversion from rₛ to rad)
        returnAvg: Bool = false
            - if true, returns the average phase for the model, otherwise returns the differential phase on every baseline
        offAxisInds: Union{Nothing,Vector{Int}} = nothing
            - if provided, only calculates the phase for off axis baselines (u,v) at indexes in offAxisInds, otherwise calculates for all (u,v) pairs
    kwargs:
    returns:
        phase: either Tuple{Vector{Float64},Vector{Float64},Vector{Float64}} or Vector{Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}
            - if returnAvg is true, returns the edges, centers, and average phase for the model
            - if returnAvg is false, returns a list of tuples for each off-axis baseline (u,v) with edges, centers, and differential phase
            - edges: Vector{Float64}
                - bin edges for the velocity
            - centers: Vector{Float64}
                - bin centers for the velocity 
            - differential phase: Vector{Float64}
                - differential phase for the model, units of radians
                - calculated as the integral of the phase over the model at each velocity bin, weighted by the intensity and ΔA (i.e. area) of the ring
                - this is the phase that would be measured by an interferometer for a given baseline (u,v) and position angle PA
    """
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

#need to add raytrace flag for this -- if true use imgs instead of getVariable by rings
function getProfile(m::model, name::Union{String,Symbol,Function}; bins::Union{Int,Vector{Float64}}=100, dx::Union{Array{Float64,},Nothing}=nothing, kwargs...)
    """
    Return a profile for the model based on the name passed. The name can be one of the following:
        - :line: returns the line profile (i.e. integrated intensity as a function of velocity)
        - :delay: returns the delay profile (i.e. mean delays weighted by intensity and/or responsivity as a function of velocity)
        - :r: returns the mean radius (as weighted by intensity) as a function of velocity 
        - :ϕ: returns the mean azimuthal angle (as weighted by intensity) as a function of velocity
        - :phase: returns the phase profile (i.e. integrated phase as a function of velocity)
            -- must pass U [Mλ], V [Mλ], PA [rad], and BLRAng [rad] as minimum extra keyword arguments; see `phase` for more details
        - Function: if a function is passed, the intensity weighted mean of this function as a function of velocity is returned
    params:
        m: model
            - model object to get the profile from
        name: Union{String,Symbol,Function}
            - name of the profile to get
        bins: Union{Int,Vector{Float64}}
            - number of bins or bin edges for binning
            - if {Int} then bins is the number of bins, and the edges of the bins will be equally spaced between the min and max of the velocity in the model
            - if {Vector{Float64} then bins are the bin edges and the nubmer of bins will be length(bins)-1
            - left edge of each bin is inclusive, right edge is exclusive except for the last bin which is inclusive
        dx: Union{Array{Float64,},Nothing}
            - if provided, this is used as the associated integral element for each ring in the model, otherwise defaults to ΔA in each ring struct
        kwargs: optional keyword arguments passed to `binModel`
            - includes: `overflow=true` to include overflow bins, `centered=true` to center the bins around 0.0, `minX` and `maxX` to set the min and max of the bins, etc.
            - see `binnedSum` for more details on available keyword arguments
    returns:
        profile: profile
            - a profile object containing the bin edges, bin centers, and binned sums
            - assign to a model object using `setProfile!` to store the profile in the model
    """
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

function setProfile!(m::model, p::profile; overwrite::Bool=false) #mutates the model.profiles dictionary based on passed on name passed (i.e. LP, delayProfile, etc.). can also take a function + kwargs to set a custom profile. stores bins + yvalues in each dict entry.
    """
    Set a profile in the model's profiles dictionary. If the profile already exists and overwrite is false, a warning is issued.
    params:
        m: model
            - model object to set the profile in
        p: profile
            - profile object to set in the model
        overwrite: Bool = false
            - if true, will overwrite the existing profile with the same name, otherwise will issue a warning if the profile already exists
    returns:
        m: model
            - the model object with the profile set in its profiles dictionary
    """
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

