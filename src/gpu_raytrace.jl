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
a list of output-grid blocks (`grids`, one per continuous submodel — each a `_RtGridBlock` with grid
edges, per-output-pixel area/camera coords, and the block's inner radius for ordering) plus the
per-point submodel index and discrete/cloud flag (`submodel`/`discrete`, length = number of input
points). `raytrace!` assembles the blocks innermost-first into a single output grid at call time, so any
number/kind of submodels combine generically (mirroring the host `_rt_build_output`). All array fields
live on the same backend as the model columns.
"""
struct _RtGridBlock{GA,V<:AbstractVector}
    grid::GA          # _RaytraceGridArrays
    pixelΔA::V
    pixelα::V
    pixelβ::V
    nPixels::Int
    minRMin::Float64  # innermost cell radius (host scalar) -- sort key for overlap resolution
end

Adapt.@adapt_structure _RtGridBlock

struct RaytraceMeta{GB,IV<:AbstractVector{Int},BV<:AbstractVector{Bool}}
    grids::GB         # Vector{_RtGridBlock} (host vector; each block holds device arrays)
    submodel::IV
    discrete::BV
end

# Adapt each grid block's device arrays (and the per-point columns) but keep `grids` a host Vector.
Adapt.adapt_structure(to, m::RaytraceMeta) = RaytraceMeta(
    map(b -> Adapt.adapt(to, b), m.grids), Adapt.adapt(to, m.submodel), Adapt.adapt(to, m.discrete))

_rt_grid_block(gridArrays, pixelΔA, pixelα, pixelβ, npix) = _RtGridBlock(
    gridArrays, pixelΔA, pixelα, pixelβ, npix,
    isempty(gridArrays.rMin) ? Inf : Float64(minimum(gridArrays.rMin)))

# Build the (host-side) raytrace metadata for a combined model, reusing the audited host output-grid
# and submodel machinery. The host `_rt_build_output` already unions all continuous submodels into one
# grid, so this is a single block. Flat arrays only -- `gpu(m)` adapts the result to the device.
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
    block = _rt_grid_block(gridArrays, pixelΔA, pixelα, pixelβ, length(pixels))
    blocks = length(pixels) > 0 ? [block] : typeof(block)[]
    return RaytraceMeta(blocks, submodel, discrete)
end

# Assemble the grid blocks into one output grid for binning, ordered innermost-first so the bin-assign
# first-match resolves overlaps to the inner grid (matching the host union). Pixel slots are offset per
# block; covered outer pixels simply receive no points and are dropped in compaction (= the union).
function _rt_assemble_grid(grids, ma)
    T = eltype(ma.I)
    if isempty(grids)
        e = _rt_backend_adapt(T[], ma); ei = _rt_backend_adapt(Int[], ma)
        return (_RaytraceGridArrays{T,typeof(e),typeof(ei)}(e, e, e, ei, ei, ei), e, e, e, 0)
    end
    sorted = sort(grids; by=g -> g.minRMin)
    cat(f) = reduce(vcat, (getfield(g.grid, f) for g in sorted))
    keyOffset = 0; startOffset = 0
    pkParts = similar(sorted, Any); psParts = similar(sorted, Any)
    for (i, g) in enumerate(sorted)
        psParts[i] = g.grid.pixelStarts .+ startOffset
        pkParts[i] = g.grid.pixelKeys .+ keyOffset
        startOffset += length(g.grid.pixelKeys)
        keyOffset += g.nPixels
    end
    pixelKeys = reduce(vcat, pkParts)
    grid = _RaytraceGridArrays{T,typeof(cat(:rMin)),typeof(pixelKeys)}(
        cat(:rMin), cat(:rMax), cat(:Δϕ), cat(:nϕ), reduce(vcat, psParts), pixelKeys)
    pΔA = reduce(vcat, (g.pixelΔA for g in sorted))
    pα = reduce(vcat, (g.pixelα for g in sorted))
    pβ = reduce(vcat, (g.pixelβ for g in sorted))
    return (grid, pΔA, pα, pβ, keyOffset)
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

function _rt_compact_resident(scanOut::_RaytraceScanOutput, pixelΔA, pixelα, pixelβ, ma::ModelArrays,
        freeI, freeMask, ::Type{T}) where {T<:Real}
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
        vcat(pixelΔA[a], ma.ΔA[freeMask]),         # output-pixel area for pixels, own ΔA for clouds
        g(scanOut.τ, ma.τ),
        g(scanOut.η, ma.η),
        g(scanOut.x, ma.x),
        g(scanOut.y, ma.y),
        g(scanOut.z, ma.z),
        vcat(pixelα[a], ma.α[freeMask]),           # camera α/β: pixel coords for pixels, own for clouds
        vcat(pixelβ[a], ma.β[freeMask]),
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

    grid, pixelΔA, pixelα, pixelβ, nPixels = _rt_assemble_grid(meta.grids, ma)
    IRd = _rt_backend_adapt(Vector{T}(_rt_iratio_vector(IRatios, rm.nSubModels)), ma)
    keys = _rt_bin_assign(ma, grid; backend=backend)
    perm = _rt_sortperm_by_key_depth(ma, keys)
    segStarts, segStops = _rt_device_segment_arrays(keys, perm, nPixels; backend=backend)
    scan = _RaytraceScanArrays{T,typeof(pixelΔA),typeof(segStarts)}(
        pixelΔA, IRd, meta.submodel, segStarts, segStops)
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

    outMA = _rt_compact_resident(scanOut, pixelΔA, pixelα, pixelβ, ma, freeI, freeKeep, T)
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

# RaytraceMeta for a freshly-constructed DiskWind submodel: one grid block, no discrete points.
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
    block = _rt_grid_block(grid, ma.ΔA, ma.α, ma.β, npix)
    return RaytraceMeta([block], submodel, discrete)
end

# RaytraceMeta for a cloud submodel: no output grid (all points discrete). Empty grids list -- a `+`
# with a continuous submodel adopts that side's grid blocks.
function _rt_cloud_meta(ma::ModelArrays, backend, ::Type{T}) where {T<:Real}
    n = length(ma.I)
    Vt = typeof(_rt_backend_adapt(T[], ma))               # backend Float vector type
    It = typeof(_rt_backend_adapt(Int[], ma))             # backend Int vector type
    blocks = _RtGridBlock{_RaytraceGridArrays{T,Vt,It},Vt}[]
    submodel = fill!(KernelAbstractions.allocate(backend, Int, n), 1)
    discrete = fill!(KernelAbstractions.allocate(backend, Bool, n), true)
    return RaytraceMeta(blocks, submodel, discrete)
end

_rt_meta_has_grid(meta::RaytraceMeta) = !isempty(meta.grids)
_rt_meta_has_grid(::Nothing) = false

# Merge two handles' raytrace metadata for `+`: concatenate the grid-block lists (any number of
# continuous submodels -- `_rt_assemble_grid` orders them innermost-first at raytrace, so overlaps
# resolve to the inner grid exactly like the host union) and the per-point submodel/discrete columns
# (submodel indices shift by the left handle's submodel count). Generic over any N and kind of submodel.
# If either side lacks metadata the result is unraytraceable (`nothing`).
function _rt_merge_meta(m1, nSub1::Int, m2, nSub2::Int)
    (m1 === nothing || m2 === nothing) && return nothing
    grids = vcat(m1.grids, m2.grids)
    submodel = vcat(m1.submodel, m2.submodel .+ nSub1)
    discrete = vcat(m1.discrete, m2.discrete)
    return RaytraceMeta(grids, submodel, discrete)
end
