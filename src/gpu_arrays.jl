using Adapt

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

function _flatten_length(m::model)
    n = 0
    for r in m.rings
        n += length(r.I)
    end
    return n
end

_ma_point_value(x, col::Int) = x isa AbstractVector ? x[col] : x

function _flatten_chunks(m::model)
    chunks = UnitRange{Int}[]
    for s in 1:length(m.subModelStartInds)
        start = m.subModelStartInds[s]
        stop = s == length(m.subModelStartInds) ? length(m.rings) : m.subModelStartInds[s+1]-1
        push!(chunks, start:stop)
    end
    return chunks
end

function _flatten_per_point(m::model, variable::Symbol, ::Type{T}) where {T<:Real}
    vals = T[]
    sizehint!(vals, _flatten_length(m))
    for chunk in _flatten_chunks(m)
        nPer = length(m.rings[first(chunk)].I)
        for col in 1:nPer, ringInd in chunk
            r = m.rings[ringInd]
            push!(vals, T(_ma_point_value(getfield(r, variable), col)))
        end
    end
    return vals
end

function _flatten_reflect(m::model)
    vals = Bool[]
    sizehint!(vals, _flatten_length(m))
    for chunk in _flatten_chunks(m)
        nPer = length(m.rings[first(chunk)].I)
        for col in 1:nPer, ringInd in chunk
            r = m.rings[ringInd]
            push!(vals, Bool(_ma_point_value(getfield(r, :reflect), col)))
        end
    end
    return vals
end

function _flatten_xyz(m::model, ::Type{T}) where {T<:Real}
    x = T[]
    y = T[]
    z = T[]
    sizehint!(x, _flatten_length(m))
    sizehint!(y, _flatten_length(m))
    sizehint!(z, _flatten_length(m))
    for chunk in _flatten_chunks(m)
        nPer = length(m.rings[first(chunk)].I)
        xyz = [getXYZ(m.rings[ringInd]) for ringInd in chunk]
        for col in 1:nPer
            for (k, ringInd) in enumerate(chunk)
                push!(x, T(_ma_point_value(xyz[k][1], col)))
                push!(y, T(_ma_point_value(xyz[k][2], col)))
                push!(z, T(_ma_point_value(xyz[k][3], col)))
            end
        end
    end
    return x, y, z
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

function flatten(m::model; T=Float64)
    T <: Real || throw(ArgumentError("flatten element type T must be <: Real, got $T"))
    x, y, z = _flatten_xyz(m, T)
    ma = ModelArrays{T,Vector{T},Vector{Bool}}(
        _flatten_per_point(m, :r, T),
        _flatten_per_point(m, :ϕ, T),
        _flatten_per_point(m, :ϕ₀, T),
        _flatten_per_point(m, :i, T),
        _flatten_per_point(m, :rot, T),
        _flatten_per_point(m, :θₒ, T),
        _flatten_per_point(m, :v, T),
        _flatten_per_point(m, :I, T),
        _flatten_per_point(m, :ΔA, T),
        _flatten_per_point(m, :τ, T),
        _flatten_per_point(m, :η, T),
        x,
        y,
        z,
        T.(vec(m.camera.α)),
        T.(vec(m.camera.β)),
        _flatten_reflect(m),
    )
    return _validate_model_arrays(ma)
end

cpu(ma::ModelArrays) = Adapt.adapt(Array, ma)
cpu(m::model; T=Float64) = flatten(m; T=T)

function gpu(::Any; kwargs...)
    error("GPU support requires loading CUDA.jl; use the CUDA extension on feature/gpu-framework Phase B")
end
