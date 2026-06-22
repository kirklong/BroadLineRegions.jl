using Adapt
using KernelAbstractions

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
    ResidentModel(ma::ModelArrays, backend, nSubModels::Int)

A model held resident on a compute device: a flat [`ModelArrays`](@ref) snapshot plus the
`KernelAbstractions` backend its columns live on. Built once (`gpu(m)` on the device, or
[`resident`](@ref) on the CPU) and reused across many observable calls — `getProfile`, `getΨ`,
`getΨt`, `phase`, `secondMoment` have methods on it that run the GPU kernels without re-flattening or
re-transferring the model each call.
"""
struct ResidentModel{MA<:ModelArrays,B}
    ma::MA
    backend::B
    nSubModels::Int
end

"""
    resident(m::model; T=Float64, backend=KernelAbstractions.CPU()) -> ResidentModel

Flatten `m` and wrap it as a [`ResidentModel`](@ref) on `backend`. The default `CPU()` backend keeps
everything on the host — useful for testing the resident pipeline without a GPU. Use `gpu(m)` (with
CUDA.jl loaded) to build a device-resident handle.
"""
resident(m::model; T=Float64, backend=KernelAbstractions.CPU()) =
    ResidentModel(flatten(m; T=T), backend, length(m.subModelStartInds))

cpu(rm::ResidentModel) = ResidentModel(cpu(rm.ma), KernelAbstractions.CPU(), rm.nSubModels)

function gpu(::Any; kwargs...)
    error("GPU support requires loading CUDA.jl; use the CUDA extension on feature/gpu-framework Phase B")
end
