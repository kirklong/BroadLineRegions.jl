# GPU-resident observables: getProfile / getΨ / getΨt / phase / secondMoment evaluated on a
# device-resident model (see ResidentModel in gpu_arrays.jl). Each method mirrors its CPU
# counterpart's binning semantics exactly -- everything (profiles AND transfer functions) now shares
# binnedSum's left-exclusive / right-inclusive edge convention via `_rt_uniform_bin_index`. Only the
# small edge/result vectors move between host and device; the heavy per-point columns stay resident.

# NaN-aware extrema computed on the device (a reduction, not a full copy back).
_rt_device_nanmin(x) = mapreduce(xi -> isnan(xi) ? typemax(eltype(x)) : xi, min, x)
_rt_device_nanmax(x) = mapreduce(xi -> isnan(xi) ? typemin(eltype(x)) : xi, max, x)

# Resolve uniform bin edges (host Vector{Float64}) for a resident-model profile, matching binnedSum:
# Int bins -> constructBinEdges from the (NaN-resolved, device-reduced) data range; an edge vector is
# accepted only if (near-)uniform, since the GPU kernels assume uniform spacing.
function _rt_resident_edges(x, bins::Union{Int,Vector{Float64}}; centered::Bool=true,
        minX::Union{Float64,Nothing}=nothing, maxX::Union{Float64,Nothing}=nothing)
    if bins isa Int
        mn = isnothing(minX) ? Float64(_rt_device_nanmin(x)) : minX
        mx = isnothing(maxX) ? Float64(_rt_device_nanmax(x)) : maxX
        return constructBinEdges(bins, mn, mx, centered)
    else
        edges = bins
        n = length(edges) - 1
        n >= 1 || throw(ArgumentError("bins edge vector must have at least 2 entries"))
        Δ = (edges[end] - edges[1]) / n
        for k in 1:n
            isapprox(edges[k+1] - edges[k], Δ; rtol=1e-9, atol=1e-12) ||
                throw(ArgumentError("GPU-resident profiles require uniform bins; pass an Int or a uniformly spaced edge vector (got non-uniform edges)"))
        end
        return collect(Float64.(edges))
    end
end

_rt_bin_centers(edges) = [(edges[k] + edges[k+1]) / 2 for k in 1:length(edges)-1]

# Device weighted histogram Σ w over x-bins, returned as host (edges, centers, sums) like binnedSum.
function _rt_resident_binsum(rm::ResidentModel, x, w, bins; overflow::Bool=false,
        centered::Bool=true, minX=nothing, maxX=nothing)
    edges = _rt_resident_edges(x, bins; centered=centered, minX=minX, maxX=maxX)
    T = eltype(rm.ma.I)
    dEdges = _rt_backend_adapt(T.(edges), rm.ma)
    out = similar(rm.ma.I, length(edges) - 1)
    _rt_weighted_histogram!(out, x, w, dEdges; overflow=overflow, backend=rm.backend)
    return edges, _rt_bin_centers(edges), Float64.(Array(out))
end

# Device weighted mean num/den over x-bins (the :delay/:r/:ϕ/Function profile shape).
function _rt_resident_weighted_mean(rm::ResidentModel, x, wNum, wDen, bins; kwargs...)
    edges, centers, num = _rt_resident_binsum(rm, x, wNum, bins; kwargs...)
    _, _, den = _rt_resident_binsum(rm, x, wDen, bins; kwargs...)
    return edges, centers, num ./ den
end

# Device per-bin weighted variance of θ (the secondMoment shape), host (edges, centers, σ²).
function _rt_resident_binvar(rm::ResidentModel, x, θ, w, bins; overflow::Bool=false,
        centered::Bool=true, minX=nothing, maxX=nothing)
    edges = _rt_resident_edges(x, bins; centered=centered, minX=minX, maxX=maxX)
    nb = length(edges) - 1
    T = eltype(rm.ma.I)
    dEdges = _rt_backend_adapt(T.(edges), rm.ma)
    σ² = similar(rm.ma.I, nb); sumW = similar(σ²); sumWθ = similar(σ²); μ = similar(σ²); sumWδ² = similar(σ²)
    _rt_weighted_variance!(σ², sumW, sumWθ, μ, sumWδ², x, θ, w, dEdges; overflow=overflow, backend=rm.backend)
    return edges, _rt_bin_centers(edges), Float64.(Array(σ²))
end

# Reverberation delays as a resident device column (τ = η(r − x); see _rt_transfer_delays).
_rt_resident_delays(rm::ResidentModel) = _rt_transfer_delays(rm.ma; backend=rm.backend)

"""
    getProfile(rm::ResidentModel, name; bins=100, dx=nothing, kwargs...) -> profile

GPU-resident `getProfile`. Supports `:line`, `:delay`, `:r`, `:ϕ`, `:phase`, `:moment2`, and a custom
`Function`, matching the host `getProfile(::model, …)` semantics (binnedSum's left-exclusive bins).
`dx` is not supported on the resident path (the model's `ΔA` is always used); pass a host `model` for
custom integration elements.
"""
function getProfile(rm::ResidentModel, name::Union{String,Symbol,Function}; bins=100,
        dx::Union{Array{Float64,},Nothing}=nothing, kwargs...)
    isnothing(dx) || throw(ArgumentError("the GPU-resident getProfile does not accept a custom dx; use the host model path"))
    ma = rm.ma
    n = name isa Function ? :function : Symbol(name)
    w = ma.I .* ma.ΔA
    if n == :line
        edges, centers, sums = _rt_resident_binsum(rm, ma.v, w, bins; kwargs...)
        return profile(name=:line, binEdges=edges, binCenters=centers, binSums=sums)
    elseif n == :delay
        delays = _rt_resident_delays(rm)
        edges, centers, vals = _rt_resident_weighted_mean(rm, ma.v, delays .* ma.I .* ma.ΔA, ma.I .* ma.η .* ma.ΔA, bins; kwargs...)
        return profile(name=:delay, binEdges=edges, binCenters=centers, binSums=vals)
    elseif n == :r
        edges, centers, vals = _rt_resident_weighted_mean(rm, ma.v, ma.r .* w, w, bins; kwargs...)
        return profile(name=:r, binEdges=edges, binCenters=centers, binSums=vals)
    elseif n == :ϕ
        edges, centers, vals = _rt_resident_weighted_mean(rm, ma.v, ma.ϕ .* w, w, bins; kwargs...)
        return profile(name=:ϕ, binEdges=edges, binCenters=centers, binSums=vals)
    elseif n == :phase
        edges, centers, avg = phase(rm; returnAvg=true, bins=bins, kwargs...)
        return profile(name=:phase, binEdges=edges, binCenters=centers, binSums=avg)
    elseif n == :moment2
        edges, centers, σ², _ = secondMoment(rm; returnAvg=true, bins=bins, kwargs...)
        return profile(name=:moment2, binEdges=edges, binCenters=centers, binSums=σ²)
    elseif name isa Function
        throw(ArgumentError("custom-Function profiles are not supported on the GPU-resident path; use the host model"))
    else
        throw(ArgumentError("profile $(name) not recognized -- choose from [:line, :delay, :r, :ϕ, :phase, :moment2]"))
    end
end

"""
    getΨt(rm::ResidentModel, tEdges::AbstractVector, overflow::Bool=false) -> Ψt

GPU-resident 1D transfer function Ψ(t); uniform `tEdges` only. Matches the host `getΨt` (left-inclusive
delay bins, empty/NaN-poisoned bins floored to 1e-30, optional overflow into the edge bins).
"""
function getΨt(rm::ResidentModel, tEdges::AbstractVector, overflow::Bool=false)
    ma = rm.ma
    delays = _rt_resident_delays(rm)
    T = eltype(ma.I)
    dEdges = _rt_backend_adapt(T.(collect(tEdges)), ma)
    Ψt = similar(ma.I, length(tEdges) - 1); under = similar(ma.I, 1); over = similar(ma.I, 1)
    _rt_psit!(Ψt, under, over, delays, ma.I, ma.ΔA, dEdges; backend=rm.backend)
    Ψth = Float64.(Array(Ψt))
    @. Ψth = ifelse(Ψth > 0, Ψth, 1e-30)
    if overflow
        sU = Float64(Array(under)[1]); sO = Float64(Array(over)[1])
        Ψth[1] += sU > 0 ? sU : 1e-30
        Ψth[end] += sO > 0 ? sO : 1e-30
    end
    return Ψth
end

"""
    getΨ(rm::ResidentModel, vEdges::AbstractVector, tEdges::AbstractVector) -> Ψ

GPU-resident 2D transfer function Ψ(v, t); uniform edges only. Matches the host `getΨ`.
"""
function getΨ(rm::ResidentModel, vEdges::AbstractVector, tEdges::AbstractVector)
    ma = rm.ma
    delays = _rt_resident_delays(rm)
    T = eltype(ma.I)
    dV = _rt_backend_adapt(T.(collect(vEdges)), ma); dT = _rt_backend_adapt(T.(collect(tEdges)), ma)
    Ψ = similar(ma.I, length(vEdges) - 1, length(tEdges) - 1)
    _rt_psi2d!(Ψ, ma.v, delays, ma.I, ma.ΔA, dV, dT; backend=rm.backend)
    Ψh = Float64.(Array(Ψ))
    @. Ψh = ifelse(Ψh > 0, Ψh, 1e-30)
    return Ψh
end

"""
    phase(rm::ResidentModel; U, V, PA, BLRAng, returnAvg=false, offAxisInds=nothing, bins=100, kwargs...)

GPU-resident differential phase, matching `phase(::model, …)`.
"""
function phase(rm::ResidentModel; returnAvg::Bool=false, offAxisInds::Union{Nothing,Vector{Int}}=nothing,
        U::Vector{Float64}, V::Vector{Float64}, PA::Float64, BLRAng::Float64, bins=100, kwargs...)
    size(U) == size(V) || throw(ArgumentError("U and V must be the same size, got $(size(U)) and $(size(V))"))
    ma = rm.ma
    X = ma.α .* eltype(ma.I)(BLRAng); Y = ma.β .* eltype(ma.I)(BLRAng)
    w = ma.I .* ma.ΔA
    inds = isnothing(offAxisInds) ? collect(1:length(U)) : offAxisInds
    results = map(inds) do k
        Ui = U[k]; Vi = V[k]
        U′ = cos(PA)*Ui + sin(PA)*Vi; V′ = -sin(PA)*Ui + cos(PA)*Vi
        edges, centers, LP = _rt_resident_binsum(rm, ma.v, w, bins; kwargs...)
        norm = maximum(LP)
        proj = eltype(ma.I)(-2π) .* (X .* eltype(ma.I)(U′) .+ Y .* eltype(ma.I)(V′)) .* w .* eltype(ma.I)(1e6)
        _, _, Δϕ = _rt_resident_binsum(rm, ma.v, proj, edges; kwargs...)
        (edges, centers, (Δϕ ./ norm) ./ (1.0 .+ LP ./ norm))
    end
    if returnAvg
        edges, centers, _ = results[1]
        avg = sum(r -> r[3], results) ./ length(results)
        return (edges, centers, avg)
    else
        return results
    end
end

"""
    secondMoment(rm::ResidentModel; U, V, PA, BLRAng, returnAvg=false, offAxisInds=nothing, bins=100, kwargs...)

GPU-resident image second moment, matching `secondMoment(::model, …)`.
"""
function secondMoment(rm::ResidentModel; returnAvg::Bool=false, offAxisInds::Union{Nothing,Vector{Int}}=nothing,
        U::Vector{Float64}, V::Vector{Float64}, PA::Float64, BLRAng::Float64, bins=100, kwargs...)
    size(U) == size(V) || throw(ArgumentError("U and V must be the same size, got $(size(U)) and $(size(V))"))
    ma = rm.ma
    Tt = eltype(ma.I)
    X = ma.α .* Tt(BLRAng); Y = ma.β .* Tt(BLRAng)
    w = ma.I .* ma.ΔA
    inds = isnothing(offAxisInds) ? collect(1:length(U)) : offAxisInds
    results = map(inds) do k
        Ui = U[k]; Vi = V[k]
        U′ = cos(PA)*Ui + sin(PA)*Vi; V′ = -sin(PA)*Ui + cos(PA)*Vi
        s = (X .* Tt(U′) .+ Y .* Tt(V′)) ./ Tt(hypot(U′, V′))
        edges, centers, σ² = _rt_resident_binvar(rm, ma.v, s, w, bins; kwargs...)
        _, _, σ²tot = _rt_resident_binvar(rm, ma.v, s, w, 1; overflow=true, centered=false)
        (edges, centers, σ², σ²tot[1])
    end
    if returnAvg
        edges, centers, _, _ = results[1]
        σ²Avg = sum(r -> r[3], results) ./ length(results)
        σ²totAvg = sum(r -> r[4], results) / length(results)
        return (edges, centers, σ²Avg, σ²totAvg)
    else
        return results
    end
end
