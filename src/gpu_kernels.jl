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
