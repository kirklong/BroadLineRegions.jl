#!/usr/bin/env julia

function binModel(;m::model, yVariables::Union{Vector{String},Vector{Symbol}}, bins::Union{Int,Vector{Float64}}=100, xVariable::Union{String,Symbol}=:v)
    """bin the model into a histogram, where each bin is the integrated value of the yVariable as a function of the xVariable
    params: 
        m: model
            - model object to bin
        yVariables: Either Vector{String}, Vector{Symbol}, or Function -- if String then will be converted to symbols
            - list of variables to bin
            - must be valid attributes of model.rings (e.g. :I, :v, :r, :e, :i, :ϕ) or a function that can be applied to model.rings
                - example: Keplerian disk time delays could be calculated like `t(ring) = ring.r*(1 .+ sin.(ring.ϕ).*ring.i))`
        bins: Union{Int,Vector{Float64}}
            - number of bins or bin edges for binning
            - if number of bins, then the bins will be equally spaced between the min and max of the xVariable
            - if bin edges then the number of bins will be length(bins)-1
                - leftmost edge is then the first bin and the rightmost edge being the last bin
            - left edge of each bin is inclusive, right edge is exclusive except for the last bin which is inclusive
        xVariable: Union{String,Symbol}
            - variable to bin over
            - must be a valid attribute of model.rings (e.g. :I, :v, :r, :e, :i, :ϕ) or a function that can be applied to model.rings
    returns:
        (binEdges, binCenters): Tuple{Vector{Float64},Vector{Float64}}
            - bin edges and centers for the xVariable of the histogram 
        yBinned: Vector{Vector{Float64}}
            - binned values of the yVariables
            - length of yBinned is the same as the number of yVariables, with each element being a vector of length equal to the number of bins
    """
end