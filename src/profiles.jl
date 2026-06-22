#!/usr/bin/env julia
"""
    constructBinEdges(bins::Int, minX::Float64, maxX::Float64, centered::Bool) -> Vector{Float64}

Build the uniform velocity-bin edges for an integer bin count, given the (NaN-resolved) data range
`minX`/`maxX`. Shared by `binnedSum` and the GPU-resident profile path so the two cannot diverge on
the `centered` edge convention.
"""
function constructBinEdges(bins::Int, minX::Float64, maxX::Float64, centered::Bool)
    if centered
        middle = (maxX+minX)/2
        Δ = (maxX-minX)/bins
        nBins = 2*ceil(Int,(maxX-(middle+Δ/2))/Δ) + 1
        maxX = (middle + Δ/2) + nBins/2*Δ
        minX = (middle - Δ/2) - nBins/2*Δ
        return collect(range(minX,stop=maxX,length=nBins+1))
    else
        return collect(range(minX,stop=maxX,length=bins+1))
    end
end

"""
    binnedSum(x::Array{Float64,}, y::Array{Float64, }; bins=100,
            overflow=false, centered=true, minX=nothing, maxX=nothing)

Bin the x and y variables into a histogram, where each bin is the sum of the y values for the corresponding x values.

# Arguments
- `x::Array{Float64,}`: x variable to bin over
- `y::Array{Float64,}`: y variable to bin
- `bins::Union{Int,Vector{Float64}}=100`: Number of bins or bin edges for binning
  - If `Int`: Number of bins with edges equally spaced between min/max of x
  - If `Vector{Float64}`: Specific bin edges, with number of bins = `length(bins)-1`
  - Bins are left-edge exclusive, right-edge inclusive: a point exactly on an interior edge goes to the bin below it.
    Points at or below the first edge, or at or above the last edge, are only counted when `overflow=true` (into the first/last bin)
- `overflow::Bool=false`: If `true`, include values outside bin range in the first/last bins
- `centered::Bool=true`: If `true`, shift bin edges to center around middle value
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
    uniform = typeof(bins) == Int #edges we construct from a range are uniform -> direct index math; user-supplied edge vectors may not be
    if typeof(bins) == Int
        if isnothing(minX)
            minX = isnan(minimum(x)) ? minimum(i for i in x if !isnan(i)) : minimum(x)
        end
        if isnothing(maxX)
            maxX = isnan(maximum(x)) ? maximum(i for i in x if !isnan(i)) : maximum(x)
        end
        binEdges = constructBinEdges(bins, minX, maxX, centered)
    else
        binEdges = bins
    end
    binCenters = @. (binEdges[1:end-1] + binEdges[2:end])/2
    result = zeros(length(binEdges)-1)
    binAccumulate!(result, x, y, binEdges, overflow, uniform)
    return (binEdges, binCenters, result)
end

"""
    binAccumulate!(result, x, y, binEdges, overflow::Bool, uniform::Bool)

Typed accumulation loop for `binnedSum` (function barrier). When `uniform` is true the bin index is
computed directly from the spacing (`floor`) and then corrected against the actual edges, which keeps
the assignment *bit-identical* to `searchsortedfirst` — including the exact-edge convention (a point
exactly on an interior edge belongs to the bin to its left) and any floating-point rounding of the
direct computation. Non-uniform (user-supplied) edges use `searchsortedfirst` as before.
"""
function binAccumulate!(result::Vector{Float64}, x::Array{Float64,}, y::Array{Float64,}, binEdges::Vector{Float64}, overflow::Bool, uniform::Bool)
    nbins = length(binEdges)-1
    if uniform
        x0 = binEdges[1]
        invΔ = nbins/(binEdges[end]-binEdges[1])
        @inbounds for (xi,yi) in zip(x, y) #zip (not eachindex) preserves the old shape-leniency for callers
            if isfinite(xi) && isfinite(yi)
                if xi <= binEdges[1]
                    if overflow
                        result[1] += yi
                    end
                elseif xi >= binEdges[end]
                    if overflow
                        result[end] += yi
                    end
                else
                    b = clamp(floor(Int,(xi-x0)*invΔ)+1, 1, nbins)
                    while b > 1 && xi <= binEdges[b] #left-exclusive: a point exactly on an edge belongs to the bin to its left
                        b -= 1
                    end
                    while b < nbins && xi > binEdges[b+1] #guard against floor() rounding one bin low
                        b += 1
                    end
                    result[b] += yi
                end
            end
        end
    else
        @inbounds for (xi,yi) in zip(x, y) #zip (not eachindex) preserves the old shape-leniency for callers
            if isfinite(xi) && isfinite(yi)
                if xi <= binEdges[1]
                    if overflow
                        result[1] += yi
                    end
                elseif xi >= binEdges[end]
                    if overflow
                        result[end] += yi
                    end
                else
                    bin = searchsortedfirst(binEdges,xi)-1
                    result[bin] += yi
                end
            end
        end
    end
    return result
end

"""
    binnedVariance(x::Array{Float64,}, θ::Array{Float64,}, w::Array{Float64,};
            bins=100, overflow=false, kwargs...)

Per-bin weighted variance of θ: bin the x variable as in `binnedSum`, then in each bin compute the
w-weighted variance ``\\sigma^2 = \\sum_k w_k(\\theta_k-\\mu)^2 / \\sum_k w_k`` with μ the w-weighted mean of θ in that bin.

The variance is computed with a two-pass algorithm (pass 1: per-bin Σw and Σwθ → per-bin mean μ;
pass 2: per-bin Σw(θ-μ)²/Σw) rather than the raw-moment identity σ² = ⟨θ²⟩ - ⟨θ⟩², which suffers
catastrophic cancellation when the per-bin mean is large compared to the spread (e.g. a few clouds
clustered tightly far from the bin mean). The two-pass form is also guaranteed non-negative.
Bin membership in pass 2 follows the same rules as `binnedSum`, so the two can never disagree.

# Arguments
- `x::Array{Float64,}`: x variable to bin over
- `θ::Array{Float64,}`: variable to take the variance of
- `w::Array{Float64,}`: weights (e.g. intensity × area)
- `bins::Union{Int,Vector{Float64}}=100`: number of bins or bin edges, as in `binnedSum`
- `overflow::Bool=false`: if `true`, include values outside the bin range in the first/last bins
- Additional `kwargs` (`centered`, `minX`, `maxX`) are passed to `binnedSum`

# Returns
- `Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}`: A tuple containing:
  - `binEdges`: Bin edges for the x variable
  - `binCenters`: Bin centers for the x variable
  - `σ²`: per-bin weighted variance of θ (NaN for empty bins)
"""
function binnedVariance(x::Array{Float64,}, θ::Array{Float64,}, w::Array{Float64,};
    bins::Union{Int,Vector{Float64}}=100, overflow::Bool=false, kwargs...)
    binEdges, binCenters, Σw = binnedSum(x, w, bins=bins, overflow=overflow; kwargs...)
    _, _, Σwθ = binnedSum(x, w.*θ, bins=binEdges, overflow=overflow; kwargs...)
    μ = Σwθ./Σw #per-bin weighted mean -- 0/0 = NaN for empty bins
    Σwδ² = zeros(length(binCenters))
    for k in eachindex(x)
        xk = x[k]
        isfinite(xk) || continue
        bin = 0 #not in any bin -- same edge rules as binnedSum
        if xk <= binEdges[1]
            bin = overflow ? 1 : 0
        elseif xk >= binEdges[end]
            bin = overflow ? length(binEdges)-1 : 0
        else
            bin = searchsortedfirst(binEdges, xk) - 1
        end
        (bin == 0 || !isfinite(θ[k]) || !isfinite(w[k]) || !isfinite(μ[bin])) && continue
        Σwδ²[bin] += w[k]*(θ[k]-μ[bin])^2
    end
    return (binEdges, binCenters, Σwδ²./Σw) #0/0 = NaN for empty bins
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
  - Bins are left-edge exclusive, right-edge inclusive: a point exactly on an interior edge goes to the bin below it.
    Points at or below the first edge, or at or above the last edge, are only counted when `overflow=true` (into the first/last bin)
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

    Calculate time delays for a point in a disk as `` t = \\eta r \\left(1 - \\cos(\\phi) \\sin(i)\\right)``
"""
tDisk(ring::ring) = @. ring.η*ring.r*(1 - cos(ring.ϕ)*sin(ring.i)) #time delays for Keplerian disk [rₛ], or a cloud modelled as a point in a disk

"""
    tCloud(ring::ring)

    Calculate time delays from the general geometric formula ``t = \\eta (r - x)``, where ``x`` is the
    3D system coordinate of the point towards the camera (cached on the ring, see `getXYZ`).

    Despite the name this is exact for *any* geometry: for a flat ring (``\\theta_o = 0``) it reduces
    algebraically to the `tDisk` formula (``x = r\\cos\\phi\\sin i``), and unlike `tDisk` it also
    accounts for lifted (``\\theta_o \\neq 0``) and reflected points. Works on both point (scalar) and
    continuous (vector) rings.
"""
tCloud(ring::ring) = begin
    x = getXYZ(ring)[1] #cached system coordinates (computed once, reused every call)
    η = ring.η; r = ring.r
    return @. η*(r - x) #could also calculate new incliation angle based on θₒ, but this is simpler, +x
end
"""
    t(ring::ring)

    Calculate time delays for a point in a disk as `` t = \\eta r \\left(1 + \\cos(\\phi) \\sin(i)\\right)`` or a cloud with opening angle ``\\theta_o``
    as the x-coordinate of the point subtracted from the radial distance of the point ``t = r - x``.
"""
t(ring::ring) = begin 
    if ring.θₒ == 0.0 #if θₒ = 0.0, then this is a point in a disk,
        return tDisk(ring)
    else #otherwise it is a cloud with opening angle θₒ
        return tCloud(ring)  
    end
end
"""
#     t(ring::ring, subFxn::Function=tDisk)

#     Calculate time delays for a point in a disk or cloud using a custom function `subFxn` that takes a `ring` struct.
#     It is more peformant to pass the function directly rather than figure it out on the fly if known ahead of time.
# """
t(ring::ring,subFxn::Function) = subFxn(ring) #allow for custom time delay function, e.g. tDisk or tCloud

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
    X = m.camera.α.*BLRAng; Y = m.camera.β.*BLRAng #DEBUG: - to make consistent with old paper
    function getΔϕ(v::Array{Float64,},I::Array{Float64,},ΔA::Array{Float64,},x::Array{Float64,},y::Array{Float64},U::Float64,V::Float64,PA::Float64;kwargs...)
        U′ = cos(PA)*U+sin(PA)*V; V′ = -sin(PA)*U+cos(PA)*V
        edges,centers,LP = binnedSum(v, I.*ΔA;kwargs...)
        norm = maximum(LP)
        LP = LP./norm #normalize the line profile to 1.0
        I = I./norm #normalize the intensity to 1.0
        edges,centers,Δϕ = binnedSum(v, -2π*(x.*U′ .+ y.*V′).*I.*ΔA.*1e6;kwargs...) 
        return (edges,centers,Δϕ./(1.0 .+ LP))#./(1.0.+LP)) #divide out LP but then weight by f/(1+f) to get the differential phase
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
    secondMoment(m::model; U, V, PA, BLRAng, returnAvg=false, offAxisInds=nothing, kwargs...)

Calculate the second central moment (squared angular size) of the model image projected along interferometric
baseline directions, both per velocity channel and integrated over the whole line.

For each baseline the per-channel size is the intensity-weighted variance of the projected position û⋅x in each
velocity bin, ``\\sigma^2_\\theta(v, \\hat{u}) = \\langle (\\hat{u}\\cdot(x - \\bar{x}_v))^2 \\rangle_v``, i.e. the squared rms angular
extent of the channel image along the baseline direction. In the marginally-resolved limit this predicts the
differential visibility-amplitude dip via ``|V_l|(v) = 1 - 2\\pi^2 |u|^2 \\sigma^2_\\theta(v)`` with ``|u| = \\sqrt{U^2+V^2}\\times10^6``
(baseline length in units of the observed wavelength).

The line-integrated size is the same variance taken over the whole line image (all velocities in a single channel),
which by the law of total variance equals the flux-weighted mean of the per-channel sizes *plus* the flux-weighted
spread of the per-channel photocenters — for a rotation-dominated BLR the photocenter-spread term is comparable to
the mean per-channel size, so the line-integrated size is significantly larger than the flux-averaged channel size.

Note that these are moments of the line image only: unlike `phase`, no continuum dilution factor f/(1+f) is applied.

# Arguments
- `m::model`: Model object to calculate image second moments for
- `U::Vector{Float64}`: U component of complex visibility in [Mλ]
- `V::Vector{Float64}`: V component of complex visibility in [Mλ]
- `PA::Float64`: On-sky position angle of the model in radians
- `BLRAng::Float64`: Characteristic size of the BLR model in radians (conversion from ``r_s`` to radians)
- `returnAvg::Bool=false`: If `true`, returns the average across all baselines
- `offAxisInds::Union{Nothing,Vector{Int}}=nothing`: If provided, only calculates moments for baselines at specified indices
- Additional `kwargs` (`bins`, `overflow`, `centered`, `minX`, `maxX`) control the velocity binning as in `binnedSum`

# Returns
- If `returnAvg=true`: `Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Float64}` containing:
  - Bin edges for velocity
  - Bin centers for velocity
  - Average per-channel squared angular size σ²(v) (in rad², NaN for empty bins)
  - Average line-integrated squared angular size σ²tot (in rad²)
- If `returnAvg=false`: `Vector{Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Float64}}` containing:
  - For each baseline, a tuple of bin edges, bin centers, per-channel σ²(v), and line-integrated σ²tot
"""
function secondMoment(m::model; returnAvg::Bool=false, offAxisInds::Union{Nothing,Vector{Int}}=nothing,
    U::Vector{Float64}, V::Vector{Float64}, PA::Float64, BLRAng::Float64, kwargs...)
    if size(U) != size(V)
        throw(ArgumentError("U and V must be the same size, got $(size(U)) and $(size(V))"))
    end
    X = m.camera.α.*BLRAng; Y = m.camera.β.*BLRAng
    function getσ²(v::Array{Float64,},I::Array{Float64,},ΔA::Array{Float64,},x::Array{Float64,},y::Array{Float64},U::Float64,V::Float64,PA::Float64;kwargs...)
        U′ = cos(PA)*U+sin(PA)*V; V′ = -sin(PA)*U+cos(PA)*V
        s = (x.*U′ .+ y.*V′)./hypot(U′,V′) #projected position along baseline direction [rad] -- baseline length divides out
        edges,centers,σ² = binnedVariance(v, s, I.*ΔA; kwargs...)
        _,_,σ²tot = binnedVariance(v, s, I.*ΔA; bins=1, overflow=true, centered=false) #"one-channel" route: whole line in one bin
        return (edges,centers,σ²,σ²tot[1])
    end
    v = getVariable(m,:v,flatten=true); I = getVariable(m,:I,flatten=true); ΔA = getVariable(m,:ΔA,flatten=true)
    if isnothing(offAxisInds)
        offAxisInds = collect(1:length(U))
    end
    U = U[offAxisInds]; V = V[offAxisInds]
    σ²List = [getσ²(v,I,ΔA,X,Y,Ui,Vi,PA;kwargs...) for (Ui,Vi) in zip(U,V)]
    if returnAvg
        edges, centers, _, _ = σ²List[1] #get the edges and centers from the first baseline
        σ²Avg = mean([σ²[3] for σ² in σ²List],dims=1)[1]
        σ²totAvg = mean([σ²[4] for σ² in σ²List])
        return (edges,centers,σ²Avg,σ²totAvg)
    else
        return σ²List
    end
end

"""
    binWeightedMean(m::model, x, yNum, yDen, bins, dx; kwargs...)

Shared helper for the weighted-mean profile options of `getProfile` (`:delay`, `:r`, `:ϕ`):
bins `yNum.*dA` and `yDen.*dA` over `x` and returns `(edges, centers, num./den)`.
Replicates `binModel`'s `ΔA` shape fix-up (per-ring scalar `ΔA` broadcast across a matrix `y`).
Operates on precomputed arrays so the gathers can come from the model's memoized cache.
"""
function binWeightedMean(m::model, x, yNum, yDen, bins, dx; kwargs...)
    dA = isnothing(dx) ? getVariable(m,:ΔA) : dx
    if size(yNum) != size(dA)
        if size(dA)[1] == size(yNum)[1] && typeof(m.rings[1].ΔA) == Float64
            dA = stack([dA for _ in 1:size(yNum)[2]],dims=2)
        else
            throw(ArgumentError("y and dA must be the same size, got $(size(yNum)) and $(size(dA))\nconsider supplying custom dx = Array{Float64,}(size(y)) (use ones for equal weighting)"))
        end
    end
    pNum = binnedSum(x, yNum.*dA, bins=bins; kwargs...)
    pDen = binnedSum(x, yDen.*dA, bins=bins; kwargs...)
    return (pNum[1], pNum[2], pNum[3]./pDen[3])
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
  - `:moment2`: Returns the per-channel squared angular size of the image along the baseline direction(s) [rad²]
    - Requires `U` [Mλ], `V` [Mλ], `PA` [rad], and `BLRAng` [rad] as keyword arguments
    - Use `secondMoment` directly to also obtain the line-integrated size (a scalar per baseline, which does not fit in a profile struct)
  - `Function`: Returns the intensity weighted mean of this function vs. velocity
- `bins`: Number of bins or bin edges for binning
  - If `Int`: Number of bins with edges equally spaced between min/max velocity
  - If `Vector{Float64}`: Specific bin edges, with number of bins = `length(bins)-1`
  - Bins are left-edge exclusive, right-edge inclusive: a point exactly on an interior edge goes to the bin below it.
    Points at or below the first edge, or at or above the last edge, are only counted when `overflow=true` (into the first/last bin)
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
        #intensity-weighted mean delay per velocity bin, built from memoized arrays
        if length(m.subModelStartInds) == 1 #single model: getVariable(m,t) resolves tDisk/tCloud once and hits m.cache (shared with getΨ/getΨt)
            delays = getVariable(m,t); I = getVariable(m,:I); η = getVariable(m,:η); v = getVariable(m,:v)
            p = binWeightedMean(m, v, delays.*I, I.*η, bins, dx; kwargs...)
        else #combined model: getVariable(m,t) resolves to the general tCloud formula for every point (same delays as getΨ/getΨt -- exact consistency); expandPerPoint aligns per-ring scalar η with flattened I
            delays = getVariable(m,t); I = getVariable(m,:I,flatten=true); η = expandPerPoint(m,:η); v = getVariable(m,:v,flatten=true)
            p = binWeightedMean(m, v, delays.*I, I.*η, bins, dx; kwargs...)
        end
    elseif n == :r
        if isnothing(dx) #fast path from memoized arrays
            if length(m.subModelStartInds) == 1 #broadcasting aligns per-ring scalar r row-wise with the I matrix
                rr = getVariable(m,:r); I = getVariable(m,:I); v = getVariable(m,:v)
            else #combined: expandPerPoint aligns per-ring scalar r with the flattened I
                rr = expandPerPoint(m,:r); I = getVariable(m,:I,flatten=true); v = getVariable(m,:v,flatten=true)
            end
            p = binWeightedMean(m, v, rr.*I, I, bins, nothing; kwargs...)
        else #user-supplied dx: numerator is unweighted by I (pre-existing behavior, kept for compatibility)
            pNum = binModel(bins,dx,m=m,yVariable=:r,xVariable=:v;kwargs...)
            pDen = binModel(bins,dx,m=m,yVariable=:I,xVariable=:v;kwargs...)
            p = (pNum[1], pNum[2], pNum[3]./pDen[3])
        end
    elseif n == :ϕ
        if isnothing(dx) #fast path from memoized arrays (ϕ is always per-point: the ring constructor asserts length(ϕ) == length(I))
            ϕϕ = getVariable(m,:ϕ); I = getVariable(m,:I); v = getVariable(m,:v)
            p = binWeightedMean(m, v, ϕϕ.*I, I, bins, nothing; kwargs...)
        else #user-supplied dx: numerator is unweighted by I (pre-existing behavior, kept for compatibility)
            pNum = binModel(bins,dx,m=m,yVariable=:ϕ,xVariable=:v;kwargs...)
            pDen = binModel(bins,dx,m=m,yVariable=:I,xVariable=:v;kwargs...)
            p = (pNum[1], pNum[2], pNum[3]./pDen[3])
        end
    elseif n == :phase
        edges,centers,avgPhase = phase(m,returnAvg=true,bins=bins;kwargs...)
        p = (edges, centers, avgPhase)
    elseif n == :moment2
        edges,centers,σ²,σ²tot = secondMoment(m,returnAvg=true,bins=bins;kwargs...)
        p = (edges, centers, σ²)
    elseif isa(name,Function)
        pNum = isnothing(dx) ? binModel(bins,m=m,yVariable=name,xVariable=:v;kwargs...) : binModel(bins,dx,m=m,yVariable=name,xVariable=:v;kwargs...)
        pDen = isnothing(dx) ? binModel(bins,m=m,yVariable=:I,xVariable=:v;kwargs...) : binModel(bins,dx,m=m,yVariable=:I,xVariable=:v;kwargs...)
        p = (pNum[1], pNum[2], pNum[3]./pDen[3])
    else
        throw(ArgumentError("profile $(name) not recognized -- choose from [:line, :delay, :r, :ϕ, :phase, :moment2] or pass a function that can be applied to model.rings"))
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

