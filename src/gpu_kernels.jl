using KernelAbstractions

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
    τ_Δv::V
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
    return _RaytraceScanArrays{T,Vector{T},Vector{Int}}(
        T.([p.ΔA for p in pixels]),
        T.(IRatios),
        _rt_tau_delta_v(m, T),
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
_rt_backend_model_arrays(m::model, backend; T=Float64) = gpu(m; T=T)
_rt_backend_adapt(x, ma::ModelArrays) = x

@kernel function _rt_segmented_scan_kernel!(outI, outv, outr, outϕ, outϕ₀, outi, outrot, outθₒ,
        outτ, outη, outx, outy, outz, outreflect, outactive, perm, keys, outputΔA, IRatios,
        τ_Δv, submodel, segmentStarts, segmentStops, r, ϕ, ϕ₀, i, rot, θₒ, v, I, ΔA, τ, η,
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
        scan.outputΔA, scan.IRatios, scan.τ_Δv, scan.submodel, scan.segmentStarts,
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
