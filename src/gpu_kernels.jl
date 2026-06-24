using KernelAbstractions
using KernelAbstractions: @atomic

function _rt_uniform_bin_index(xi, edges, nbins::Int, overflow::Bool)
    isfinite(xi) || return 0
    if xi <= edges[1]
        return overflow ? 1 : 0
    elseif xi >= edges[nbins+1]
        return overflow ? nbins : 0
    else
        x0 = edges[1]
        invΔ = nbins / (edges[nbins+1] - x0)
        bin = clamp(floor(Int, (xi - x0) * invΔ) + 1, 1, nbins)
        while bin > 1 && xi <= edges[bin]
            bin -= 1
        end
        while bin < nbins && xi > edges[bin+1]
            bin += 1
        end
        return bin
    end
end

@kernel function _rt_weighted_histogram_kernel!(out, x, y, edges, overflow, nbins)
    idx = @index(Global)
    xi = x[idx]
    yi = y[idx]
    if isfinite(xi) && isfinite(yi)
        bin = _rt_uniform_bin_index(xi, edges, nbins, overflow)
        bin != 0 && @atomic out[bin] += yi
    end
end

function _rt_weighted_histogram!(out::AbstractVector, x::AbstractArray, y::AbstractArray,
        edges::AbstractVector; overflow::Bool=false, backend=KernelAbstractions.CPU())
    length(out) == length(edges) - 1 || throw(DimensionMismatch("histogram output must have length length(edges)-1"))
    length(x) == length(y) || throw(DimensionMismatch("histogram inputs must have matching lengths"))
    fill!(out, zero(eltype(out)))
    kernel! = _rt_weighted_histogram_kernel!(backend)
    event = kernel!(out, x, y, edges, overflow, length(out); ndrange=length(x))
    event !== nothing && wait(event)
    return out
end

@kernel function _rt_weighted_product_histogram_kernel!(out, x, y1, y2, edges, overflow, nbins)
    idx = @index(Global)
    xi = x[idx]
    yi = y1[idx] * y2[idx]
    if isfinite(xi) && isfinite(yi)
        bin = _rt_uniform_bin_index(xi, edges, nbins, overflow)
        bin != 0 && @atomic out[bin] += yi
    end
end

function _rt_weighted_product_histogram!(out::AbstractVector, x::AbstractArray, y1::AbstractArray,
        y2::AbstractArray, edges::AbstractVector; overflow::Bool=false, backend=KernelAbstractions.CPU())
    length(out) == length(edges) - 1 || throw(DimensionMismatch("histogram output must have length length(edges)-1"))
    length(x) == length(y1) == length(y2) || throw(DimensionMismatch("histogram inputs must have matching lengths"))
    fill!(out, zero(eltype(out)))
    kernel! = _rt_weighted_product_histogram_kernel!(backend)
    event = kernel!(out, x, y1, y2, edges, overflow, length(out); ndrange=length(x))
    event !== nothing && wait(event)
    return out
end

_rt_line_profile!(out::AbstractVector, ma::ModelArrays, edges::AbstractVector;
    overflow::Bool=false, backend=KernelAbstractions.CPU()) =
    _rt_weighted_product_histogram!(out, ma.v, ma.I, ma.ΔA, edges; overflow=overflow, backend=backend)

@kernel function _rt_weighted_mean_pass1_kernel!(sumW, sumWθ, x, θ, w, edges, overflow, nbins)
    idx = @index(Global)
    xi = x[idx]
    θi = θ[idx]
    wi = w[idx]
    if isfinite(xi) && isfinite(wi)
        bin = _rt_uniform_bin_index(xi, edges, nbins, overflow)
        if bin != 0
            @atomic sumW[bin] += wi
            wθ = wi * θi
            isfinite(wθ) && @atomic sumWθ[bin] += wθ
        end
    end
end

@kernel function _rt_finalize_mean_kernel!(μ, sumW, sumWθ)
    idx = @index(Global)
    μ[idx] = sumWθ[idx] / sumW[idx]
end

@kernel function _rt_weighted_variance_pass2_kernel!(sumWδ², x, θ, w, μ, edges, overflow, nbins)
    idx = @index(Global)
    xi = x[idx]
    θi = θ[idx]
    wi = w[idx]
    if isfinite(xi) && isfinite(θi) && isfinite(wi)
        bin = _rt_uniform_bin_index(xi, edges, nbins, overflow)
        if bin != 0 && isfinite(μ[bin])
            δ = θi - μ[bin]
            @atomic sumWδ²[bin] += wi * δ * δ
        end
    end
end

@kernel function _rt_finalize_variance_kernel!(σ², sumWδ², sumW)
    idx = @index(Global)
    σ²[idx] = sumWδ²[idx] / sumW[idx]
end

function _rt_weighted_variance!(σ²::AbstractVector, sumW::AbstractVector, sumWθ::AbstractVector,
        μ::AbstractVector, sumWδ²::AbstractVector, x::AbstractArray, θ::AbstractArray,
        w::AbstractArray, edges::AbstractVector; overflow::Bool=false, backend=KernelAbstractions.CPU())
    nbins = length(edges) - 1
    length(σ²) == length(sumW) == length(sumWθ) == length(μ) == length(sumWδ²) == nbins ||
        throw(DimensionMismatch("variance outputs must all have length length(edges)-1"))
    length(x) == length(θ) == length(w) || throw(DimensionMismatch("variance inputs must have matching lengths"))
    fill!(sumW, zero(eltype(sumW)))
    fill!(sumWθ, zero(eltype(sumWθ)))
    fill!(μ, zero(eltype(μ)))
    fill!(sumWδ², zero(eltype(sumWδ²)))
    fill!(σ², zero(eltype(σ²)))

    pass1! = _rt_weighted_mean_pass1_kernel!(backend)
    event = pass1!(sumW, sumWθ, x, θ, w, edges, overflow, nbins; ndrange=length(x))
    event !== nothing && wait(event)

    mean! = _rt_finalize_mean_kernel!(backend)
    event = mean!(μ, sumW, sumWθ; ndrange=nbins)
    event !== nothing && wait(event)

    pass2! = _rt_weighted_variance_pass2_kernel!(backend)
    event = pass2!(sumWδ², x, θ, w, μ, edges, overflow, nbins; ndrange=length(x))
    event !== nothing && wait(event)

    final! = _rt_finalize_variance_kernel!(backend)
    event = final!(σ², sumWδ², sumW; ndrange=nbins)
    event !== nothing && wait(event)
    return σ²
end

# Reverberation delay τ = η·(r − x) (the general `tCloud` form; equals `tDisk` to rounding for a flat
# ring). One thread per point; reproduces `getVariable(m, t)` for the resident pipeline without a
# host round-trip. See the W1 unified-delay note (combined models already use this form).
@kernel function _rt_transfer_delays_kernel!(delays, η, r, x)
    idx = @index(Global)
    delays[idx] = η[idx] * (r[idx] - x[idx])
end

function _rt_transfer_delays!(delays::AbstractArray, η::AbstractArray, r::AbstractArray,
        x::AbstractArray; backend=KernelAbstractions.CPU())
    size(delays) == size(η) == size(r) == size(x) || throw(DimensionMismatch("delay arrays must have matching sizes"))
    kernel! = _rt_transfer_delays_kernel!(backend)
    event = kernel!(delays, η, r, x; ndrange=length(delays))
    event !== nothing && wait(event)
    return delays
end

_rt_transfer_delays(ma::ModelArrays; backend=KernelAbstractions.CPU()) =
    _rt_transfer_delays!(similar(ma.r), ma.η, ma.r, ma.x; backend=backend)

# 2D transfer function Ψ(v, t): single-pass atomic accumulation of I·ΔA into (velocity, delay) bins,
# uniform edges only, matching getΨ's `ΨAccumulate!` (binnedSum left-exclusive interior edges, no
# overflow bins so boundary points are dropped, NaN v/delay dropped, NaN weights poison their bin so
# the caller's `>0 ? Ψ : 1e-30` floor maps them to 1e-30).
@kernel function _rt_psi2d_kernel!(Ψ, v, delays, I, ΔA, vEdges, tEdges, nV, nT)
    idx = @index(Global)
    bv = _rt_uniform_bin_index(v[idx], vEdges, nV, false)
    bt = _rt_uniform_bin_index(delays[idx], tEdges, nT, false)
    if bv != 0 && bt != 0
        @atomic Ψ[bv + (bt - 1) * nV] += I[idx] * ΔA[idx]
    end
end

function _rt_psi2d!(Ψ::AbstractMatrix, v::AbstractArray, delays::AbstractArray, I::AbstractArray,
        ΔA::AbstractArray, vEdges::AbstractVector, tEdges::AbstractVector;
        backend=KernelAbstractions.CPU())
    nV = length(vEdges) - 1
    nT = length(tEdges) - 1
    size(Ψ) == (nV, nT) || throw(DimensionMismatch("Ψ must be ($nV, $nT)"))
    length(v) == length(delays) == length(I) == length(ΔA) || throw(DimensionMismatch("Ψ inputs must have matching lengths"))
    fill!(Ψ, zero(eltype(Ψ)))
    Ψflat = reshape(Ψ, nV * nT)
    kernel! = _rt_psi2d_kernel!(backend)
    event = kernel!(Ψflat, v, delays, I, ΔA, vEdges, tEdges, nV, nT; ndrange=length(v))
    event !== nothing && wait(event)
    return Ψ
end

# 1D transfer function Ψ(t): I·ΔA binned by delay (uniform edges, same binnedSum convention as getΨt),
# with optional overflow accumulation of under-/over-range weight into the first/last bin. `under`/
# `over` are length-1 device arrays so the single sweep also collects the overflow sums.
@kernel function _rt_psit_kernel!(Ψt, under, over, delays, I, ΔA, tEdges, nT)
    idx = @index(Global)
    di = delays[idx]
    if !isnan(di)
        w = I[idx] * ΔA[idx]
        if di <= tEdges[1]
            @atomic under[1] += w
        elseif di >= tEdges[nT+1]
            @atomic over[1] += w
        else
            @atomic Ψt[_rt_uniform_bin_index(di, tEdges, nT, false)] += w
        end
    end
end

function _rt_psit!(Ψt::AbstractVector, under::AbstractVector, over::AbstractVector,
        delays::AbstractArray, I::AbstractArray, ΔA::AbstractArray, tEdges::AbstractVector;
        backend=KernelAbstractions.CPU())
    nT = length(tEdges) - 1
    length(Ψt) == nT || throw(DimensionMismatch("Ψt must have length length(tEdges)-1"))
    length(delays) == length(I) == length(ΔA) || throw(DimensionMismatch("Ψt inputs must have matching lengths"))
    fill!(Ψt, zero(eltype(Ψt)))
    fill!(under, zero(eltype(under)))
    fill!(over, zero(eltype(over)))
    kernel! = _rt_psit_kernel!(backend)
    event = kernel!(Ψt, under, over, delays, I, ΔA, tEdges, nT; ndrange=length(delays))
    event !== nothing && wait(event)
    return Ψt, under, over
end

function _rt_disk_deproject_scalar(a, b, inc, rot, θₒ, m11, m12, m21, m22,
        r3d11, r3d12, r3d21, r3d22, r3d31, r3d32, rMin, rMax, ηₒ, η₁, αRM, rNorm)
    cosr = cos(rot)
    sinr = sin(rot)
    cosi = cos(inc)
    sini = sin(inc)
    cosθₒ = cos(θₒ)
    sinθₒ = sin(θₒ)
    xRing = -(b*cosr - a*cosi*sinr) / (cosi*cosθₒ + cosr*sini*sinθₒ)
    yRing = (a*(cosi*cosθₒ + sini/cosr*sinθₒ) + b*cosθₒ*sinr/cosr) /
        (cosi*cosθₒ/cosr + sini*sinθₒ)
    r = sqrt(xRing*xRing + yRing*yRing)
    if r < rMin || r > rMax
        nan = convert(typeof(r), NaN)
        return nan, nan, nan, nan, nan, nan, nan
    end
    ϕ₀ = atan(yRing, xRing)
    x = r3d11*xRing + r3d12*yRing
    y = r3d21*xRing + r3d22*yRing
    z = r3d31*xRing + r3d32*yRing
    xd = m11*xRing + m12*yRing
    yd = m21*xRing + m22*yRing
    ϕ = atan(yd, xd)
    η = ηₒ + η₁ * (r/rNorm)^αRM
    return r, ϕ, ϕ₀, η, x, y, z
end

@kernel function _rt_disk_deproject_kernel!(rSystem, ϕSystem, ϕ₀, η, xSystem, ySystem, zSystem,
        α, β, inc, rot, θₒ, m11, m12, m21, m22, r3d11, r3d12, r3d21, r3d22,
        r3d31, r3d32, rMin, rMax, ηₒ, η₁, αRM, rNorm)
    idx = @index(Global)
    rt, ϕt, ϕ₀t, ηt, xt, yt, zt = _rt_disk_deproject_scalar(α[idx], β[idx], inc, rot, θₒ,
        m11, m12, m21, m22, r3d11, r3d12, r3d21, r3d22, r3d31, r3d32, rMin, rMax,
        ηₒ, η₁, αRM, rNorm)
    rSystem[idx] = rt
    ϕSystem[idx] = ϕt
    ϕ₀[idx] = ϕ₀t
    η[idx] = ηt
    xSystem[idx] = xt
    ySystem[idx] = yt
    zSystem[idx] = zt
end

function _rt_disk_deproject!(rSystem::AbstractArray, ϕSystem::AbstractArray, ϕ₀::AbstractArray,
        η::AbstractArray, xSystem::AbstractArray, ySystem::AbstractArray, zSystem::AbstractArray,
        α::AbstractArray, β::AbstractArray, inc::Real, rot::Real, θₒ::Real, M::AbstractMatrix,
        r3D::AbstractMatrix, rMin::Real, rMax::Real, ηₒ::Real, η₁::Real, αRM::Real, rNorm::Real;
        backend=KernelAbstractions.CPU())
    size(rSystem) == size(α) == size(β) || throw(DimensionMismatch("disk deprojection arrays must have matching sizes"))
    kernel! = _rt_disk_deproject_kernel!(backend)
    event = kernel!(rSystem, ϕSystem, ϕ₀, η, xSystem, ySystem, zSystem, α, β, inc, rot, θₒ,
        M[1,1], M[1,2], M[2,1], M[2,2], r3D[1,1], r3D[1,2], r3D[2,1], r3D[2,2],
        r3D[3,1], r3D[3,2], round(rMin, sigdigits=9), round(rMax, sigdigits=9),
        ηₒ, η₁, αRM, rNorm; ndrange=length(α))
    event !== nothing && wait(event)
    return rSystem, ϕSystem, ϕ₀, η, xSystem, ySystem, zSystem
end

_rt_v_circular_disk_scalar(r, ϕ, inc, rₛ) = -sqrt(rₛ/(2*r)) * sin(inc) * sin(ϕ)

@kernel function _rt_v_circular_disk_kernel!(v, r, ϕ, inc, rₛ)
    idx = @index(Global)
    v[idx] = _rt_v_circular_disk_scalar(r[idx], ϕ[idx], inc, rₛ)
end

function _rt_v_circular_disk!(v::AbstractArray, r::AbstractArray, ϕ::AbstractArray, inc::Real;
        rₛ=1.0, backend=KernelAbstractions.CPU())
    size(v) == size(r) == size(ϕ) || throw(DimensionMismatch("velocity arrays must have matching sizes"))
    kernel! = _rt_v_circular_disk_kernel!(backend)
    event = kernel!(v, r, ϕ, inc, rₛ; ndrange=length(v))
    event !== nothing && wait(event)
    return v
end

function _rt_disk_wind_i_scalar(r, ϕ, inc, f1, f2, f3, f4, α)
    sinϕ = sin(ϕ)
    cosϕ = cos(ϕ)
    sini = sin(inc)
    cosi = cos(inc)
    pre = sqrt(1 / (2 * r^3))
    term12 = (3*sini^2) * cosϕ * (sqrt(2)*f1*cosϕ + f2/2*sinϕ)
    term3 = (-f3*3*sini*cosi) * cosϕ
    term4 = sqrt(2)*f4*cosi^2
    return r^(-α) * abs(pre * (term12 + term3 + term4))
end

@kernel function _rt_disk_wind_i_kernel!(I, r, ϕ, inc, f1, f2, f3, f4, α, rMin, rMax)
    idx = @index(Global)
    rv = r[idx]
    if rv >= rMin && rv <= rMax
        I[idx] = _rt_disk_wind_i_scalar(rv, ϕ[idx], inc, f1, f2, f3, f4, α)
    elseif isnan(rv)
        I[idx] = rv
    else
        I[idx] = zero(eltype(I))
    end
end

function _rt_disk_wind_i!(I::AbstractArray, r::AbstractArray, ϕ::AbstractArray, inc::Real,
        f1::Real, f2::Real, f3::Real, f4::Real, α::Real, rMin::Real, rMax::Real;
        backend=KernelAbstractions.CPU())
    size(I) == size(r) == size(ϕ) || throw(DimensionMismatch("intensity arrays must have matching sizes"))
    kernel! = _rt_disk_wind_i_kernel!(backend)
    event = kernel!(I, r, ϕ, inc, f1, f2, f3, f4, α, round(rMin, sigdigits=9),
        round(rMax, sigdigits=9); ndrange=length(I))
    event !== nothing && wait(event)
    return I
end

struct _RaytraceGridArrays{T<:Real,V<:AbstractVector{T},I<:AbstractVector{Int}}
    rMin::V
    rMax::V
    Δϕ::V
    nϕ::I
    pixelStarts::I
    pixelKeys::I
end

Adapt.@adapt_structure _RaytraceGridArrays

struct _RaytraceScanArrays{T<:Real,V<:AbstractVector{T},I<:AbstractVector{Int}}
    outputΔA::V
    IRatios::V
    submodel::I
    segmentStarts::I
    segmentStops::I
end

Adapt.@adapt_structure _RaytraceScanArrays

struct _RaytraceScanOutput{T<:Real,V<:AbstractVector{T},B<:AbstractVector{Bool}}
    I::V
    v::V
    r::V
    ϕ::V
    ϕ₀::V
    i::V
    rot::V
    θₒ::V
    τ::V
    η::V
    x::V
    y::V
    z::V
    reflect::B
    active::B
end

Adapt.@adapt_structure _RaytraceScanOutput

function _rt_grid_arrays(grids::Vector{_RaytraceGrid}; T=Float64)
    rMin = T[]
    rMax = T[]
    Δϕ = T[]
    nϕ = Int[]
    pixelStarts = Int[]
    pixelKeys = Int[]
    for grid in grids
        for ringPos in eachindex(grid.rMin)
            push!(rMin, T(grid.rMin[ringPos]))
            push!(rMax, T(grid.rMax[ringPos]))
            push!(Δϕ, T(grid.Δϕ))
            push!(nϕ, grid.nϕ)
            push!(pixelStarts, length(pixelKeys) + 1)
            for col in 1:grid.nϕ
                push!(pixelKeys, grid.localToPixel[(ringPos, col)])
            end
        end
    end
    return _RaytraceGridArrays{T,Vector{T},Vector{Int}}(rMin, rMax, Δϕ, nϕ, pixelStarts, pixelKeys)
end

# Test-only oracle: CPU reference for the GPU bin keys (used by test/gpu_*.jl, no src/ callers).
function _rt_reference_pixel_keys(points::Vector{_RaytracePoint}, grids::Vector{_RaytraceGrid})
    return [_rt_find_pixel(p, grids) for p in points]
end

@kernel function _rt_bin_assign_kernel!(keys, α, β, rMin, rMax, Δϕ, nϕ, pixelStarts, pixelKeys, nIntervals)
    idx = @index(Global)
    a = α[idx]
    b = β[idx]
    key = 0
    if isfinite(a) && isfinite(b)
        rCam = sqrt(a*a + b*b)
        ϕCam = atan(b, a)
        interval = 1
        while interval <= nIntervals
            if rCam >= rMin[interval] && rCam <= rMax[interval]
                dϕ = Δϕ[interval]
                shifted = mod(atan(sin(ϕCam), cos(ϕCam)) - (-π - dϕ/2), 2π)
                col = floor(Int, shifted/dϕ) + 1
                col = col > nϕ[interval] ? 1 : col
                key = pixelKeys[pixelStarts[interval] + col - 1]
                break
            end
            interval += 1
        end
    end
    keys[idx] = key
end

function _rt_bin_assign!(keys::AbstractVector{Int}, ma::ModelArrays, gridArrays::_RaytraceGridArrays;
        backend=KernelAbstractions.CPU())
    n = length(ma.I)
    length(keys) == n || throw(DimensionMismatch("keys has length $(length(keys)) but ModelArrays has length $n"))
    length(ma.α) == n && length(ma.β) == n || throw(DimensionMismatch("ModelArrays camera columns differ from I length $n"))
    kernel! = _rt_bin_assign_kernel!(backend)
    event = kernel!(keys, ma.α, ma.β, gridArrays.rMin, gridArrays.rMax, gridArrays.Δϕ,
        gridArrays.nϕ, gridArrays.pixelStarts, gridArrays.pixelKeys, length(gridArrays.rMin); ndrange=n)
    event !== nothing && wait(event)
    return keys
end

function _rt_bin_assign(ma::ModelArrays, gridArrays::_RaytraceGridArrays; backend=KernelAbstractions.CPU())
    keys = similar(ma.I, Int, length(ma.I))
    return _rt_bin_assign!(keys, ma, gridArrays; backend=backend)
end

_rt_sort_depth(x) = isfinite(x) ? x : -Inf

function _rt_sortperm_by_key_depth(keys::AbstractVector{Int}, x::AbstractVector)
    length(keys) == length(x) || throw(DimensionMismatch("keys has length $(length(keys)) but x has length $(length(x))"))
    depthOrder = sortperm(eachindex(keys); by=i -> _rt_sort_depth(x[i]), rev=true, alg=MergeSort)
    keyOrder = sortperm(keys[depthOrder]; alg=MergeSort)
    return depthOrder[keyOrder]
end

function _rt_sortperm_by_key_depth(ma::ModelArrays, keys::AbstractVector{Int})
    return _rt_sortperm_by_key_depth(keys, ma.x)
end

# Test-only oracle: CPU reference for the GPU sort order (used by test/gpu_*.jl, no src/ callers).
function _rt_sorted_key_depth_pairs(keys::AbstractVector{Int}, x::AbstractVector, perm::AbstractVector{Int})
    return [(keys[i], _rt_sort_depth(x[i])) for i in perm]
end

function _rt_submodel_indices(m::model)
    vals = Int[]
    sizehint!(vals, _flatten_length(m))
    for s in 1:length(m.subModelStartInds)
        rr = _rt_submodel_ring_range(m, s)
        nPer = _rt_point_length(m.rings[first(rr)])
        for _ in 1:nPer, _ in rr
            push!(vals, s)
        end
    end
    return vals
end

function _rt_tau_delta_v(m::model, ::Type{T}=Float64) where {T<:Real}
    vals = T[]
    sizehint!(vals, _flatten_length(m))
    for s in 1:length(m.subModelStartInds)
        rr = _rt_submodel_ring_range(m, s)
        nPer = _rt_point_length(m.rings[first(rr)])
        for col in 1:nPer
            for ringInd in rr
                ring = m.rings[ringInd]
                vVal = T(_rt_point_value(ring.v, col))
                τ_Δv = if ring.τ isa AbstractVector && nPer > 1
                    v2 = col == 1 ? T(_rt_point_value(ring.v, min(col+1, nPer))) : T(_rt_point_value(ring.v, col-1))
                    abs((vVal + v2)/2)
                else
                    T(Inf)
                end
                push!(vals, τ_Δv)
            end
        end
    end
    return vals
end

function _rt_scan_arrays(m::model, keys::AbstractVector{Int}, perm::AbstractVector{Int},
        pixels::Vector{_RaytracePixel}, IRatios::AbstractVector{<:Real}; T=Float64)
    nPixels = length(pixels)
    segmentStarts = zeros(Int, nPixels)
    segmentStops = zeros(Int, nPixels)
    for (pos, idx) in enumerate(perm)
        key = keys[idx]
        if 1 <= key <= nPixels
            segmentStarts[key] == 0 && (segmentStarts[key] = pos)
            segmentStops[key] = pos
        end
    end
    # Same guard as the CPU scan (`_rt_scan_bucket`): velocity-dependent optical depth is not yet
    # implemented, so error host-side before launching any kernel (kernels cannot `error`).
    _rt_velocity_dependent_τ(_rt_tau_delta_v(m, T)) && error(_RT_VELOCITY_τ_MSG)
    return _RaytraceScanArrays{T,Vector{T},Vector{Int}}(
        T.([p.ΔA for p in pixels]),
        T.(IRatios),
        _rt_submodel_indices(m),
        segmentStarts,
        segmentStops,
    )
end

function _rt_scan_output(n::Int; T=Float64)
    nan = T(NaN)
    return _RaytraceScanOutput{T,Vector{T},Vector{Bool}}(
        fill(nan, n), fill(nan, n), fill(nan, n), fill(nan, n), fill(nan, n),
        fill(nan, n), fill(nan, n), fill(nan, n), fill(nan, n), fill(nan, n),
        fill(nan, n), fill(nan, n), fill(nan, n), falses(n), falses(n),
    )
end

_rt_backend_model_arrays(m::model, ::KernelAbstractions.CPU; T=Float64) = flatten(m; T=T)
_rt_backend_model_arrays(m::model, backend; T=Float64) =
    error("no device ModelArrays transfer defined for backend $(typeof(backend)); load the matching extension (e.g. `using CUDA` for CUDABackend)")
_rt_backend_adapt(x, ma::ModelArrays) = x

@kernel function _rt_segmented_scan_kernel!(outI, outv, outr, outϕ, outϕ₀, outi, outrot, outθₒ,
        outτ, outη, outx, outy, outz, outreflect, outactive, perm, keys, outputΔA, IRatios,
        submodel, segmentStarts, segmentStops, r, ϕ, ϕ₀, i, rot, θₒ, v, I, ΔA, τ, η,
        x, α, β, reflect, τCutOff)
    pix = @index(Global)
    start = segmentStarts[pix]
    stop = segmentStops[pix]
    nan = convert(eltype(outI), NaN)
    outI[pix] = nan
    outv[pix] = nan
    outr[pix] = nan
    outϕ[pix] = nan
    outϕ₀[pix] = nan
    outi[pix] = nan
    outrot[pix] = nan
    outθₒ[pix] = nan
    outτ[pix] = nan
    outη[pix] = nan
    outx[pix] = nan
    outy[pix] = nan
    outz[pix] = nan
    outreflect[pix] = false
    outactive[pix] = false

    if start != 0
        firstIdx = 0
        pos = start
        while pos <= stop
            idx = perm[pos]
            if keys[idx] == pix && isfinite(I[idx]) && isfinite(x[idx])
                firstIdx = idx
                break
            end
            pos += 1
        end

        if firstIdx != 0
            outArea = outputΔA[pix]
            firstWeight = I[firstIdx] * IRatios[submodel[firstIdx]] * ΔA[firstIdx] / outArea
            newτ = τ[firstIdx]
            newI = firstWeight
            newv = v[firstIdx] * firstWeight
            newr = r[firstIdx] * firstWeight
            newϕ = ϕ[firstIdx] * firstWeight
            newϕ₀ = ϕ₀[firstIdx] * firstWeight
            newi = i[firstIdx] * firstWeight
            newrot = rot[firstIdx] * firstWeight
            newθₒ = θₒ[firstIdx] * firstWeight
            newη = η[firstIdx] * firstWeight
            newx = x[firstIdx] * firstWeight
            newy = α[firstIdx] * firstWeight
            newz = β[firstIdx] * firstWeight
            maxArea = ΔA[firstIdx]
            firstReflect = reflect[firstIdx]

            pos += 1
            while pos <= stop
                idx = perm[pos]
                if keys[idx] == pix && isfinite(I[idx]) && isfinite(x[idx])
                    area = ΔA[idx]
                    obscuredFrac = area / maxArea
                    (newτ < τCutOff || obscuredFrac > one(obscuredFrac)) || break
                    weight = I[idx] * IRatios[submodel[idx]] * area / outArea
                    tmpI = if obscuredFrac > one(obscuredFrac)
                        obscuredI = weight / obscuredFrac
                        unobscuredI = weight - obscuredI
                        exp(-newτ) * obscuredI + unobscuredI
                    else
                        exp(-newτ) * weight
                    end
                    newv += v[idx] * tmpI
                    newr += r[idx] * tmpI
                    newϕ += ϕ[idx] * tmpI
                    newϕ₀ += ϕ₀[idx] * tmpI
                    newi += i[idx] * tmpI
                    newrot += rot[idx] * tmpI
                    newθₒ += θₒ[idx] * tmpI
                    newη += η[idx] * tmpI
                    newx += x[idx] * tmpI
                    newy += α[idx] * tmpI
                    newz += β[idx] * tmpI
                    newI += tmpI
                    newτ += τ[idx]
                    maxArea = max(maxArea, area)
                end
                pos += 1
            end

            den = newI == zero(newI) ? one(newI) : newI
            outI[pix] = newI
            outv[pix] = newv / den
            outr[pix] = newr / den
            outϕ[pix] = newϕ / den
            outϕ₀[pix] = newϕ₀ / den
            outi[pix] = newi / den
            outrot[pix] = newrot / den
            outθₒ[pix] = newθₒ / den
            outτ[pix] = newτ
            outη[pix] = newη / den
            outx[pix] = newx / den
            outy[pix] = newy / den
            outz[pix] = newz / den
            outreflect[pix] = firstReflect
            outactive[pix] = true
        end
    end
end

function _rt_segmented_scan!(out::_RaytraceScanOutput, ma::ModelArrays, keys::AbstractVector{Int},
        perm::AbstractVector{Int}, scan::_RaytraceScanArrays; τCutOff=1.0,
        backend=KernelAbstractions.CPU())
    nPixels = length(scan.outputΔA)
    length(out.I) == nPixels || throw(DimensionMismatch("scan output has length $(length(out.I)) but there are $nPixels pixels"))
    nPixels == 0 && return out
    kernel! = _rt_segmented_scan_kernel!(backend)
    event = kernel!(out.I, out.v, out.r, out.ϕ, out.ϕ₀, out.i, out.rot, out.θₒ,
        out.τ, out.η, out.x, out.y, out.z, out.reflect, out.active, perm, keys,
        scan.outputΔA, scan.IRatios, scan.submodel, scan.segmentStarts,
        scan.segmentStops, ma.r, ma.ϕ, ma.ϕ₀, ma.i, ma.rot, ma.θₒ, ma.v, ma.I,
        ma.ΔA, ma.τ, ma.η, ma.x, ma.α, ma.β, ma.reflect, τCutOff; ndrange=nPixels)
    event !== nothing && wait(event)
    return out
end

function _rt_scan_output_result(out::_RaytraceScanOutput, pixInd::Int)
    return (I=out.I[pixInd], v=out.v[pixInd], r=out.r[pixInd], ϕ=out.ϕ[pixInd],
        ϕ₀=out.ϕ₀[pixInd], i=out.i[pixInd], rot=out.rot[pixInd], θₒ=out.θₒ[pixInd],
        τ=out.τ[pixInd], η=out.η[pixInd], x=out.x[pixInd], y=out.y[pixInd],
        z=out.z[pixInd], reflect=out.reflect[pixInd])
end

function _rt_apply_scan_output!(outRings::Vector{ring}, pixels::Vector{_RaytracePixel}, out::_RaytraceScanOutput)
    for pixInd in eachindex(pixels)
        out.active[pixInd] || continue
        pix = pixels[pixInd]
        _rt_set_point!(outRings[pix.ring], pix.col, _rt_scan_output_result(out, pixInd))
        outRings[pix.ring].y[pix.col] = pix.α
        outRings[pix.ring].z[pix.col] = pix.β
    end
    return outRings
end

function _rt_finalize_model(outRings::Vector{ring}, subStarts::Vector{Int}, freeRings)
    if !isempty(freeRings)
        push!(subStarts, length(outRings)+1)
        append!(outRings, freeRings)
    end
    outRings = _rt_compact_rings(outRings)
    isempty(outRings) && error("raytrace! produced an empty model")
    αout, βout, subStarts = _rt_rebuild_camera(outRings)
    out = model(outRings, Dict{Symbol,profile}(), camera(αout, βout, true), subStarts)
    out.cache = Dict{Any,Array}()
    return out
end

function _rt_backend_raytrace(m::model, IR::Vector{Float64}, τCutOff::Float64, raytraceFreeClouds::Bool,
        backend; T=Float64)
    camStartInds = getFlattenedCameraIndices(m)
    points = _rt_flatten_points(m, camStartInds)
    grids, pixels, outRings, _, _, subStarts = _rt_build_output(m, camStartInds)

    ma = _rt_backend_model_arrays(m, backend; T=T)
    gridArrays = _rt_backend_adapt(_rt_grid_arrays(grids; T=T), ma)
    keys = _rt_bin_assign(ma, gridArrays; backend=backend)
    perm = _rt_sortperm_by_key_depth(ma, keys)

    keysCPU = Adapt.adapt(Array, keys)
    permCPU = Adapt.adapt(Array, perm)
    scan = _rt_backend_adapt(_rt_scan_arrays(m, keysCPU, permCPU, pixels, IR; T=T), ma)
    scanOut = _rt_backend_adapt(_rt_scan_output(length(pixels); T=T), ma)
    _rt_segmented_scan!(scanOut, ma, keys, perm, scan; τCutOff=T(τCutOff), backend=backend)
    _rt_apply_scan_output!(outRings, pixels, Adapt.adapt(Array, scanOut))

    freeInds = Int[]
    for (idx, p) in enumerate(points)
        keysCPU[idx] == 0 && p.discrete && isfinite(p.I) && push!(freeInds, idx)
    end
    freeRings = if raytraceFreeClouds
        _rt_attenuate_free_clouds(points, freeInds, IR, τCutOff)
    else
        [_rt_copy_cloud_point(points[idx], points[idx].I * IR[points[idx].submodel]) for idx in freeInds]
    end
    return _rt_finalize_model(outRings, subStarts, freeRings)
end
