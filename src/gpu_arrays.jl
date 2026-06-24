using Adapt
using KernelAbstractions

"""
    ModelArrays{T,V,B}

Flat structure-of-arrays representation of a [`model`](@ref), with one entry per model point for
geometry, velocity, intensity, response, camera coordinates, and reflection state. This is the common
host/device payload used by [`resident`](@ref), [`gpu`](@ref), and the resident observable kernels.
"""
struct ModelArrays{T<:Real,V<:AbstractVector{T},B<:AbstractVector{Bool}}
    r::V
    ϕ::V
    ϕ₀::V
    i::V
    rot::V
    θₒ::V
    v::V
    I::V
    ΔA::V
    τ::V
    η::V
    x::V
    y::V
    z::V
    α::V
    β::V
    reflect::B
end

Adapt.@adapt_structure ModelArrays

# Number of flattened points (one per (ring, col)); used by the kernel helpers in gpu_kernels.jl.
_flatten_length(m::model) = sum(r -> length(r.I), m.rings)

# Populate every ring's cached 3D (x,y,z) so that `expandPerPoint(m, :x/:y/:z)` can read them.
# getXYZ is a no-op when the cache is already filled (W1 stores it at construction) and otherwise
# computes it once via rotate3D, matching the plan's "x,y,z from W1 storage or rotate3D per point".
function _ensure_xyz!(m::model)
    for r in m.rings
        getXYZ(r)
    end
    return m
end

function _validate_model_arrays(ma::ModelArrays)
    lengths = (; r=length(ma.r), ϕ=length(ma.ϕ), ϕ₀=length(ma.ϕ₀), i=length(ma.i),
        rot=length(ma.rot), θₒ=length(ma.θₒ), v=length(ma.v), I=length(ma.I),
        ΔA=length(ma.ΔA), τ=length(ma.τ), η=length(ma.η), x=length(ma.x),
        y=length(ma.y), z=length(ma.z), α=length(ma.α), β=length(ma.β),
        reflect=length(ma.reflect))
    expected = lengths.I
    all(==(expected), values(lengths)) || error("ModelArrays column lengths differ: $lengths")
    return ma
end

"""
    flatten(m::model; T=Float64) -> ModelArrays{T,Vector{T},Vector{Bool}}

Build a flat structure-of-arrays snapshot of `m`, one entry per (ring, point), as the host-side
bridge to the GPU (`gpu(m)` adapts the result to `CuArray`).

Each column is gathered through the package's existing, audited expansion machinery rather than a
GPU-specific copy: per-point fields and per-ring scalars are aligned elementwise by
[`expandPerPoint`](@ref) (so e.g. a scalar `η` is broadcast to every point of its ring), the 3D
coordinates come from the ring `getXYZ` cache, and the camera columns from `m.camera`. The point
ordering matches `getVariable(m, sym, flatten=true)`.
"""
function flatten(m::model; T=Float64)
    T <: Real || throw(ArgumentError("flatten element type T must be <: Real, got $T"))
    _ensure_xyz!(m)
    col(sym) = T.(expandPerPoint(m, sym))
    ma = ModelArrays{T,Vector{T},Vector{Bool}}(
        col(:r), col(:ϕ), col(:ϕ₀), col(:i), col(:rot), col(:θₒ),
        col(:v), col(:I), col(:ΔA), col(:τ), col(:η),
        col(:x), col(:y), col(:z),
        T.(vec(m.camera.α)), T.(vec(m.camera.β)),
        Bool.(expandPerPoint(m, :reflect)),
    )
    return _validate_model_arrays(ma)
end

cpu(ma::ModelArrays) = Adapt.adapt(Array, ma)
cpu(m::model; T=Float64) = flatten(m; T=T)

"""
    ResidentModel(ma::ModelArrays, backend, nSubModels::Int [, rt])

A model held resident on a compute device: a flat [`ModelArrays`](@ref) snapshot plus the
`KernelAbstractions` backend its columns live on. Built once (`gpu(m)` on the device, or
[`resident`](@ref) on the CPU) and reused across many observable calls — `getProfile`, `getΨ`,
`getΨt`, `phase`, `secondMoment` have methods on it that run the GPU kernels without re-flattening or
re-transferring the model each call.

The optional `rt` field carries device-resident raytrace metadata (output grid + per-point submodel /
discrete info; see `RaytraceMeta`) so `raytrace!(::ResidentModel)` can run the bin→sort→scan→compact
entirely on `backend`. It is `nothing` for handles that were not built for raytracing; the observable
methods ignore it.
"""
struct ResidentModel{MA<:ModelArrays,B,RT}
    ma::MA
    backend::B
    nSubModels::Int
    rt::RT
end

# Back-compat: a handle without raytrace metadata (every existing call site builds one this way).
ResidentModel(ma::ModelArrays, backend, nSubModels::Int) = ResidentModel(ma, backend, nSubModels, nothing)

"""
    resident(m::model; T=Float64, backend=KernelAbstractions.CPU(), raytrace=false) -> ResidentModel

Flatten `m` and wrap it as a [`ResidentModel`](@ref) on `backend`. The default `CPU()` backend keeps
everything on the host — useful for testing the resident pipeline without a GPU. Use `gpu(m)` (with
CUDA.jl loaded) to build a device-resident handle.

Pass `raytrace=true` to attach the metadata that `raytrace!(::ResidentModel)` needs. It is only built for
host models with more than one submodel (single-submodel handles keep `rt === nothing`); on-device
constructors carry their own metadata for later combination via `+`.
"""
function resident(m::model; T=Float64, backend=KernelAbstractions.CPU(), raytrace::Bool=false)
    ma = flatten(m; T=T)
    meta = (raytrace && length(m.subModelStartInds) > 1) ? _rt_build_meta(m; T=T) : nothing
    return ResidentModel(ma, backend, length(m.subModelStartInds), meta)
end

cpu(rm::ResidentModel) =
    ResidentModel(cpu(rm.ma), KernelAbstractions.CPU(), rm.nSubModels, Adapt.adapt(Array, rm.rt))

# Concatenate two flat snapshots column-by-column. `vcat` runs on whatever array type the columns
# already are, so for device columns (CuArray) it stays on the device — no host round-trip.
function _cat_model_arrays(a::ModelArrays, b::ModelArrays)
    eltype(a.I) == eltype(b.I) ||
        throw(ArgumentError("cannot combine ModelArrays with different element types: $(eltype(a.I)) and $(eltype(b.I))"))
    c(f) = vcat(getfield(a, f), getfield(b, f))
    r = c(:r)
    reflect = c(:reflect)
    return ModelArrays{eltype(a.I),typeof(r),typeof(reflect)}(
        r, c(:ϕ), c(:ϕ₀), c(:i), c(:rot), c(:θₒ), c(:v), c(:I), c(:ΔA), c(:τ), c(:η),
        c(:x), c(:y), c(:z), c(:α), c(:β), reflect)
end

"""
    +(rm1::ResidentModel, rm2::ResidentModel) -> ResidentModel

Combine two device-resident models by concatenating their columns **on their existing backend** — when
both live on the GPU the `vcat`s stay on the device, so this is the fast path for `gpu(m1) + gpu(m2)`
(or for merging submodels built on-device with [`residentDiskWindModel`](@ref)/[`residentCloudModel`](@ref))
without a host round-trip. Mirrors `+(::model, ::model)`: submodels are stacked and `nSubModels` adds.
Both handles must share the same backend type and element type.

Like the host `+`, this stacks submodels without raytracing, but it **does** preserve raytrace metadata
when both operands carry it (merged via `_rt_merge_meta`). Use `raytrace!(rm)` for the device-resident
τ-scan, or host `raytrace!(m1 + m2; backend=…)` for the one-shot host-orchestrated path.
"""
function Base.:+(rm1::ResidentModel, rm2::ResidentModel)
    if typeof(rm1.backend) != typeof(rm2.backend)
        throw(ArgumentError(string(
            "cannot combine ResidentModels that live on different backends ",
            "($(typeof(rm1.backend)) + $(typeof(rm2.backend))) — move them onto the same device first, then add.\n",
            "  • to combine on the CPU:  cpu(rm) brings a device-resident model back to the host (e.g. `cpu(gpuModel) + cpuModel`)\n",
            "  • to combine on the GPU:  build/move both with gpu(m) (e.g. `gpuModel + gpu(cpu(otherModel))`)")))
    end
    meta = _rt_merge_meta(rm1.rt, rm1.nSubModels, rm2.rt, rm2.nSubModels)
    return ResidentModel(_cat_model_arrays(rm1.ma, rm2.ma), rm1.backend,
        rm1.nSubModels + rm2.nSubModels, meta)
end

# Mixing a host `model` (lives on the CPU) with a device-resident handle is the other "different
# places" mistake — give the same actionable message instead of a bare MethodError.
const _RESIDENT_MIX_MSG = string(
    "cannot combine a host `model` with a `ResidentModel` — they live in different places. ",
    "Put both in the same representation first:\n",
    "  • flatten the host model with `resident(m)` (CPU) or `gpu(m)` (GPU) and add two ResidentModels, or\n",
    "  • bring the resident handle back with `cpu(rm)` and add two host models.")
Base.:+(::model, ::ResidentModel) = throw(ArgumentError(_RESIDENT_MIX_MSG))
Base.:+(::ResidentModel, ::model) = throw(ArgumentError(_RESIDENT_MIX_MSG))

"""
    gpu(m::model; T=Float32) -> ResidentModel
    gpu(ma::ModelArrays) -> ModelArrays

Move a host model or flat [`ModelArrays`](@ref) snapshot onto the GPU and return a device-resident
handle for repeated observable calls. Requires CUDA.jl to be loaded so the package extension can
provide the CUDA-backed methods.
"""
function gpu(::Any; kwargs...)
    error("GPU support requires CUDA.jl — run `using CUDA` (with a functional CUDA device) to activate the BroadLineRegions CUDA extension")
end
