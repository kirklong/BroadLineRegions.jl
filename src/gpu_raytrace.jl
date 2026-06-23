using KernelAbstractions
using Adapt

# Device-resident raytrace: take a ResidentModel and run the whole bin -> sort -> segmented-scan ->
# compact pipeline on its backend, with NO host round-trips during the call, producing a raytraced
# ResidentModel ready for the observables pipeline. This is the device dispatch of the mutating
# `raytrace!` generic (`raytrace!(rm::ResidentModel)`); like `raytrace!(m::model)` it returns a freshly
# combined handle (the point count changes: one combined point per active output pixel + surviving free
# clouds), so callers write `rm = raytrace!(rm)`.
#
# The metadata the bare flat columns lost (output grid edges, per-output-pixel ΔA/α/β, per-point
# submodel index and discrete/cloud flag) is carried on `rm.rt::RaytraceMeta`. It is built once when the
# handle is created -- on the host inside `gpu(m)` (cheap flat arrays, NOT during the raytrace call) or
# on-device by the construction path -- so the per-iteration raytrace stays entirely on the backend.

"""
    RaytraceMeta

Device-resident metadata that `raytrace!(::ResidentModel)` needs beyond the flat [`ModelArrays`](@ref):
the output-grid edges (`grid`, a `_RaytraceGridArrays`), the per-output-pixel area/camera coordinates
(`pixelΔA`/`pixelα`/`pixelβ`, length = number of output pixels), and the per-point submodel index and
discrete/cloud flag (`submodel`/`discrete`, length = number of input points). All fields live on the
same backend as the model columns.
"""
struct RaytraceMeta{GA,V<:AbstractVector,IV<:AbstractVector{Int},BV<:AbstractVector{Bool}}
    grid::GA
    pixelΔA::V
    pixelα::V
    pixelβ::V
    submodel::IV
    discrete::BV
    nPixels::Int
end

Adapt.@adapt_structure RaytraceMeta

# Build the (host-side) raytrace metadata for a combined model, reusing the audited host output-grid
# and submodel machinery. Flat arrays only (no per-point structs kept) -- `gpu(m)` adapts the result to
# the device alongside the columns, so the cost is paid once at transfer time, not per raytrace.
function _rt_build_meta(m::model; T=Float64)
    camStartInds = getFlattenedCameraIndices(m)
    grids, pixels, _, _, _, _ = _rt_build_output(m, camStartInds)
    gridArrays = _rt_grid_arrays(grids; T=T)
    pixelΔA = T[p.ΔA for p in pixels]
    pixelα = T[p.α for p in pixels]
    pixelβ = T[p.β for p in pixels]
    submodel = _rt_submodel_indices(m)
    discreteSub = [_rt_is_discrete(m, s) for s in 1:length(m.subModelStartInds)]
    discrete = Bool[discreteSub[s] for s in submodel]
    return RaytraceMeta(gridArrays, pixelΔA, pixelα, pixelβ, submodel, discrete, length(pixels))
end

# Segment boundaries from the key-sorted permutation. `perm` orders points by (key ascending, depth
# descending), so every output pixel's contributors are contiguous -- each pixel's first/last position
# is written exactly once, no atomics. keys==0 (NaN / free clouds) are skipped.
@kernel function _rt_segment_bounds_kernel!(segStarts, segStops, keys, perm, n, nPixels)
    pos = @index(Global)
    k = keys[perm[pos]]
    if 1 <= k <= nPixels
        if pos == 1 || keys[perm[pos-1]] != k
            segStarts[k] = pos
        end
        if pos == n || keys[perm[pos+1]] != k
            segStops[k] = pos
        end
    end
end

function _rt_device_segment_arrays(keys::AbstractVector{Int}, perm::AbstractVector{Int}, nPixels::Int;
        backend=KernelAbstractions.CPU())
    segStarts = KernelAbstractions.zeros(backend, Int, nPixels)
    segStops = KernelAbstractions.zeros(backend, Int, nPixels)
    n = length(keys)
    if n > 0 && nPixels > 0
        kernel! = _rt_segment_bounds_kernel!(backend)
        event = kernel!(segStarts, segStops, keys, perm, n, nPixels; ndrange=n)
        event !== nothing && wait(event)
    end
    return segStarts, segStops
end

# Allocate a NaN/false-initialized scan output of length n directly on `backend` (the segmented-scan
# kernel overwrites every slot, so the contents only need to be the right type/size).
function _rt_scan_output_device(backend, ::Type{T}, n::Int) where {T<:Real}
    f() = KernelAbstractions.allocate(backend, T, n)
    b() = KernelAbstractions.allocate(backend, Bool, n)
    return _RaytraceScanOutput{T,typeof(f()),typeof(b())}(
        f(), f(), f(), f(), f(), f(), f(), f(), f(), f(), f(), f(), f(), b(), b())
end

# Gather one column for the raytraced output: the active output pixels (from the segmented-scan result)
# followed by the surviving free clouds. `pixCol` is the per-pixel value, `cloudCol` the per-point value.
_rt_gather_out(pixCol, active, cloudCol, freeMask) = vcat(pixCol[active], cloudCol[freeMask])

function _rt_compact_resident(scanOut::_RaytraceScanOutput, meta::RaytraceMeta, ma::ModelArrays,
        keys::AbstractVector{Int}, IRpp, freeI, freeMask, ::Type{T}) where {T<:Real}
    a = scanOut.active
    g(pix, cloud) = _rt_gather_out(pix, a, cloud, freeMask)
    r = g(scanOut.r, ma.r)
    reflect = vcat(scanOut.reflect[a], ma.reflect[freeMask])
    out = ModelArrays{T,typeof(r),typeof(reflect)}(
        r,
        g(scanOut.ϕ, ma.ϕ),
        g(scanOut.ϕ₀, ma.ϕ₀),
        g(scanOut.i, ma.i),
        g(scanOut.rot, ma.rot),
        g(scanOut.θₒ, ma.θₒ),
        g(scanOut.v, ma.v),
        vcat(scanOut.I[a], freeI),                 # free-cloud I already scaled by IRatios
        vcat(meta.pixelΔA[a], ma.ΔA[freeMask]),    # output-pixel area for pixels, own ΔA for clouds
        g(scanOut.τ, ma.τ),
        g(scanOut.η, ma.η),
        g(scanOut.x, ma.x),
        g(scanOut.y, ma.y),
        g(scanOut.z, ma.z),
        vcat(meta.pixelα[a], ma.α[freeMask]),      # camera α/β: pixel coords for pixels, own for clouds
        vcat(meta.pixelβ[a], ma.β[freeMask]),
        reflect)
    return out
end

function _rt_resident_raytrace(rm::ResidentModel; IRatios::Union{Float64,Array{Float64,}}=1.0,
        τCutOff::Float64=1.0, raytraceFreeClouds::Bool=false)
    meta = rm.rt
    meta === nothing && error("raytrace!(::ResidentModel) needs raytrace metadata; build the handle " *
        "with gpu(m) (CUDA loaded) or the on-device constructors so it carries an output grid")
    ma = rm.ma
    backend = rm.backend
    T = eltype(ma.I)
    nPixels = meta.nPixels

    IRd = _rt_backend_adapt(Vector{T}(_rt_iratio_vector(IRatios, rm.nSubModels)), ma)
    keys = _rt_bin_assign(ma, meta.grid; backend=backend)
    perm = _rt_sortperm_by_key_depth(ma, keys)
    segStarts, segStops = _rt_device_segment_arrays(keys, perm, nPixels; backend=backend)
    scan = _RaytraceScanArrays{T,typeof(meta.pixelΔA),typeof(segStarts)}(
        meta.pixelΔA, IRd, meta.submodel, segStarts, segStops)
    scanOut = _rt_scan_output_device(backend, T, nPixels)
    _rt_segmented_scan!(scanOut, ma, keys, perm, scan; τCutOff=T(τCutOff), backend=backend)

    # surviving free clouds: out-of-grid (key==0), discrete, finite I. For the default path they pass
    # through with I scaled by their submodel's IRatio; with `raytraceFreeClouds` each is dimmed by the
    # τ of overlapping clouds in front of it (camera disks of radius √(ΔA/π)) and culled past τCutOff.
    IRpp = IRd[meta.submodel]
    freeMask = (keys .== 0) .& meta.discrete .& isfinite.(ma.I)
    freeKeep, freeI = if raytraceFreeClouds
        _rt_resident_free_clouds(ma, meta, freeMask, IRpp, T(τCutOff), backend)
    else
        (freeMask, ma.I[freeMask] .* IRpp[freeMask])
    end

    outMA = _rt_compact_resident(scanOut, meta, ma, keys, IRpp, freeI, freeKeep, T)
    nsub = length(freeI) > 0 ? 2 : 1
    return ResidentModel(outMA, backend, nsub)
end

# Cloud–cloud attenuation on the device (port of `_rt_attenuate_free_clouds`). Each free cloud is a
# camera disk of radius √(ΔA/π); it is dimmed by `exp(-Σ τ)` over overlapping clouds with larger depth
# `x` (ties broken by index, matching the host) and culled once that Σ τ exceeds `τCutOff`. Neighbour
# search is a uniform spatial hash on (α,β) sized to the max cloud radius, but with **no dense cell
# table** (the host uses a sparse Dict; clouds spread over hundreds of rₛ would need millions of cells):
# clouds are sorted by cell id and each cloud binary-searches the sorted ids for its 3×3 cell stencil —
# O(nF) memory, independent of the bounding box. Returns `(freeKeep, freeI)`: a full-length mask of the
# surviving free clouds (ascending global order) and their dimmed, IRatio-scaled intensities.

# First index `i` with `arr[i] >= val` (1-based; `n+1` if none) — device-safe binary search.
@inline function _rt_lower_bound(arr, val, n::Int)
    lo = 1; hi = n + 1
    while lo < hi
        mid = (lo + hi) >>> 1
        if arr[mid] < val
            lo = mid + 1
        else
            hi = mid
        end
    end
    return lo
end

@kernel function _rt_free_cloud_kernel!(τfront, αf, βf, xf, τf, radii, perm, sortedCellId,
        minx, miny, nx, ny, cellSize, nF)
    p = @index(Global)
    T = eltype(τfront)
    αp = αf[p]; βp = βf[p]; xp = xf[p]; rp = radii[p]
    cxp = unsafe_trunc(Int, floor(αp / cellSize))
    cyp = unsafe_trunc(Int, floor(βp / cellSize))
    acc = zero(T)
    for dj in -1:1, di in -1:1
        gx = cxp + di - minx
        gy = cyp + dj - miny
        if 0 <= gx < nx && 0 <= gy < ny
            cid = gx + gy * nx + 1
            lo = _rt_lower_bound(sortedCellId, cid, nF)
            hi = _rt_lower_bound(sortedCellId, cid + 1, nF) - 1
            pos = lo
            while pos <= hi
                q = perm[pos]
                if q != p
                    dα = αp - αf[q]; dβ = βp - βf[q]
                    dist = sqrt(dα * dα + dβ * dβ)
                    if dist < rp + radii[q] && (xf[q] > xp || (xf[q] == xp && q < p))
                        acc += τf[q]
                    end
                end
                pos += 1
            end
        end
    end
    τfront[p] = acc
end

function _rt_resident_free_clouds(ma::ModelArrays, meta::RaytraceMeta, freeMask, IRpp, τCutOff,
        backend)
    T = eltype(ma.I)
    idx = findall(freeMask)                         # ascending global indices of the free clouds
    nF = length(idx)
    nF == 0 && return (freeMask, similar(ma.I, 0))

    αf = ma.α[idx]; βf = ma.β[idx]; xf = ma.x[idx]; τf = ma.τ[idx]
    If = ma.I[idx] .* IRpp[idx]
    radii = sqrt.(ma.ΔA[idx] ./ T(π))
    cellSize = maximum(radii) * 2
    cellSize = cellSize == zero(T) ? one(T) : cellSize

    cx = unsafe_trunc.(Int, floor.(αf ./ cellSize))
    cy = unsafe_trunc.(Int, floor.(βf ./ cellSize))
    minx = minimum(cx); miny = minimum(cy)
    nx = maximum(cx) - minx + 1; ny = maximum(cy) - miny + 1
    cellId = (cx .- minx) .+ (cy .- miny) .* nx .+ 1
    perm = _rt_sortperm_by_key_depth(cellId, xf)     # group equal cellIds contiguously (total order)
    sortedCellId = cellId[perm]

    τfront = similar(If)
    kernel! = _rt_free_cloud_kernel!(backend)
    event = kernel!(τfront, αf, βf, xf, τf, radii, perm, sortedCellId, minx, miny, nx, ny,
        cellSize, nF; ndrange=nF)
    event !== nothing && wait(event)

    keep = τfront .<= τCutOff
    freeI = (If .* exp.(.-τfront))[keep]
    freeKeep = fill!(similar(freeMask), false)
    freeKeep[idx[keep]] .= true
    return (freeKeep, freeI)
end

"""
    raytrace!(rm::ResidentModel; IRatios=1.0, τCutOff=1.0, raytraceFreeClouds=false) -> ResidentModel

Device-resident dispatch of [`raytrace!`](@ref BLR.raytrace!): run the full bin→sort→segmented-scan→
compact combine on `rm`'s backend with no host round-trips, returning a freshly raytraced
[`ResidentModel`](@ref) (drop-in for the observable methods). `rm` must carry raytrace metadata — build
it with `gpu(m)` (CUDA loaded) or `resident(m; raytrace=true)`.

Like `raytrace!(::model)` this returns a new combined handle rather than mutating in place (the point
count changes: one combined point per active output pixel + surviving free clouds), so write
`rm = raytrace!(rm)`. Semantics match the host `raytrace!`: `IRatios` is a global per-submodel weight,
`τCutOff` the optical-depth cutoff, `raytraceFreeClouds` enables cloud–cloud attenuation.
"""
raytrace!(rm::ResidentModel; IRatios::Union{Float64,Array{Float64,}}=1.0, τCutOff::Float64=1.0,
        raytraceFreeClouds::Bool=false) =
    _rt_resident_raytrace(rm; IRatios=IRatios, τCutOff=τCutOff, raytraceFreeClouds=raytraceFreeClouds)

# ======================================================================================
# Phase 2: raytrace metadata built ON-DEVICE by the construction path, so a model assembled with
# residentDiskWindModel/residentCloudModel and combined with `+` can be raytraced without any host
# build. The output-pixel grid for a single DiskWind submodel is its own camera grid (pixel slot =
# disk point's flat index, ΔA/α/β = its columns), so the only new on-device array is `pixelKeys`.
# ======================================================================================

# Per-ring radial cell edges of the DiskWind camera grid (matches `_rt_ring_edges`): nr host values.
function _rt_diskwind_grid_edges(rMin::Real, rMax::Real, inc::Real, nr::Int, nϕ::Int, scale::Symbol,
        ::Type{T}) where {T<:Real}
    scaleLog = scale === :log
    rStart = scaleLog ? log(rMin * cos(inc)) : rMin * cos(inc)
    rEnd = scaleLog ? log(rMax) : rMax
    Δr = (rEnd - rStart) / (nr - 1)
    rMinRing = Vector{T}(undef, nr)
    rMaxRing = Vector{T}(undef, nr)
    for k in 1:nr
        rCam = scaleLog ? exp(rStart + (k - 1) * Δr) : (rStart + (k - 1) * Δr)
        ΔrUp = scaleLog ? rCam * (exp(Δr) - 1) : Δr
        ΔrDown = scaleLog ? min(rCam, rCam * (1 - 1 / exp(Δr))) : min(rCam, Δr)
        rMinRing[k] = rCam - ΔrDown / 2
        rMaxRing[k] = rCam + ΔrUp / 2
    end
    return rMinRing, rMaxRing, T(2π / nϕ)
end

# pixelKeys[(k-1)*nϕ + j] = k + (j-1)*nr  (cell (ring k, col j) -> its pixel slot = the disk flat index)
@kernel function _rt_disk_pixelkeys_kernel!(pixelKeys, nr, nϕ)
    idx = @index(Global)
    k = (idx - 1) ÷ nϕ + 1
    j = (idx - 1) % nϕ + 1
    pixelKeys[idx] = k + (j - 1) * nr
end

# RaytraceMeta for a freshly-constructed DiskWind submodel (one grid; no discrete points).
function _rt_diskwind_meta(ma::ModelArrays, rMin::Real, rMax::Real, inc::Real, nr::Int, nϕ::Int,
        scale::Symbol, backend, ::Type{T}) where {T<:Real}
    rMinRing, rMaxRing, Δϕv = _rt_diskwind_grid_edges(rMin, rMax, inc, nr, nϕ, scale, T)
    pixelStarts = Int[(k - 1) * nϕ + 1 for k in 1:nr]
    npix = nr * nϕ
    pixelKeys = KernelAbstractions.allocate(backend, Int, npix)
    kernel! = _rt_disk_pixelkeys_kernel!(backend)
    event = kernel!(pixelKeys, nr, nϕ; ndrange=npix)
    event !== nothing && wait(event)
    grid = _RaytraceGridArrays{T,typeof(ma.ΔA),typeof(pixelKeys)}(
        _rt_backend_adapt(rMinRing, ma), _rt_backend_adapt(rMaxRing, ma),
        _rt_backend_adapt(fill(Δϕv, nr), ma), _rt_backend_adapt(fill(nϕ, nr), ma),
        _rt_backend_adapt(pixelStarts, ma), pixelKeys)
    submodel = fill!(KernelAbstractions.allocate(backend, Int, npix), 1)
    discrete = fill!(KernelAbstractions.allocate(backend, Bool, npix), false)
    return RaytraceMeta(grid, ma.ΔA, ma.α, ma.β, submodel, discrete, npix)
end

# RaytraceMeta for a cloud submodel: no output grid (all points discrete). Empty grid arrays so the
# merge in `+` can adopt a real grid from a continuous submodel on the other side.
function _rt_cloud_meta(ma::ModelArrays, backend, ::Type{T}) where {T<:Real}
    n = length(ma.I)
    emptyT = _rt_backend_adapt(T[], ma)
    emptyI = _rt_backend_adapt(Int[], ma)
    grid = _RaytraceGridArrays{T,typeof(emptyT),typeof(emptyI)}(emptyT, emptyT, emptyT, emptyI, emptyI, emptyI)
    submodel = fill!(KernelAbstractions.allocate(backend, Int, n), 1)
    discrete = fill!(KernelAbstractions.allocate(backend, Bool, n), true)
    return RaytraceMeta(grid, emptyT, emptyT, emptyT, submodel, discrete, 0)
end

_rt_meta_has_grid(meta::RaytraceMeta) = meta.nPixels > 0
_rt_meta_has_grid(::Nothing) = false

# Merge two handles' raytrace metadata for `+`. The combined output grid is the (single) continuous
# submodel's grid (pixel slots are independent of point order); submodel indices shift by the left
# handle's submodel count and discrete flags concatenate. Two continuous grids (disk+disk) would need an
# on-device grid union -- not done here; combine those on the host (`gpu(m1 + m2)`).
function _rt_merge_meta(m1, nSub1::Int, m2, nSub2::Int)
    (m1 === nothing && m2 === nothing) && return nothing
    g1 = _rt_meta_has_grid(m1)
    g2 = _rt_meta_has_grid(m2)
    g1 && g2 && throw(ArgumentError(
        "combining two continuous (grid) submodels on-device is not supported; build on the host and " *
        "use `gpu(m1 + m2)` for the device-resident raytrace of multi-disk models"))
    sub1 = m1 === nothing ? nothing : m1.submodel
    sub2 = m2 === nothing ? nothing : m2.submodel
    sub1 === nothing && (sub1 = similar(sub2, 0))
    sub2 === nothing && (sub2 = similar(sub1, 0))
    disc1 = m1 === nothing ? similar(sub1, Bool, 0) : m1.discrete
    disc2 = m2 === nothing ? similar(sub2, Bool, 0) : m2.discrete
    submodel = vcat(sub1, sub2 .+ nSub1)
    discrete = vcat(disc1, disc2)
    src = g1 ? m1 : m2
    return RaytraceMeta(src.grid, src.pixelΔA, src.pixelα, src.pixelβ, submodel, discrete, src.nPixels)
end
