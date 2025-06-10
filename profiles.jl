#!/usr/bin/env julia
function binnedSum(x::Array{Float64,}, y::Array{Float64, }; bins::Union{Int,Vector{Float64}}=100, overflow::Bool = false, centered::Bool=false, minX::Union{Float64,Nothing}=nothing, maxX::Union{Float64,Nothing}=nothing)
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

t(ring::ring) = ring.η.*ring.r.*(1 .+ cos.(ring.ϕ).*sin(ring.i)) # time delays for Keplerian disk [rₛ], or a cloud modelled as a point in a disk

#tDisk(ring::ring) = ring.η.*ring.r.*(1 .+ cos.(ring.ϕ).*sin(ring.i)) # time delays for Keplerian disk [rₛ]
# tCloud(ring::ring) = begin -- redundant, this is equivalent to above
#     xyzSys = rotate3D(ring.r,ring.ϕ₀,ring.i,ring.rot,ring.θₒ,ring.reflect) 
#     ring.η*(ring.r - xyzSys[1]) # time delays for clouds [rₛ]
# end

#need to add raytrace flag for this -- if true use imgs instead of getVariable by rings
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
        #IMPLEMENT PHASE PROFILE FOR INTERFEROMETRY

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

