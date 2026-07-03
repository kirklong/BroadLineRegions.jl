#!/usr/bin/env julia
"""
    CompositeModel

Holds an ordered collection of `line => model` pairs so that several broad emission
lines (e.g. Hα + Hβ + CIV) can be modeled together. The existing `model`/`ring`/
`raytrace!`/cache/GPU code paths are unmodified and know nothing about lines --
`CompositeModel` is purely a wrapper around independently valid `model` objects.

# Fields
- `lines::Vector{String}`: unique line names, in creation order. `lines[1]` is the
  *reference* line (the implicit `from` default for [`addLine!`](@ref) parameter reuse,
  and the line every `fluxRatio` is relative to).
- `models::Dict{String,model}`: the per-line `model` object.
- `lineCenters::Dict{String,Float64}`: central wavelength per line, in any units --
  consistent units across lines in one `CompositeModel` are the caller's
  responsibility (only Δλ/λ matters downstream).
- `fluxRatios::Dict{String,Float64}`: each line's velocity-integrated line flux
  relative to `lines[1]` (so `fluxRatios[lines[1]] == 1.0`).

There is deliberately no `lines(cm)` accessor -- `cm.lines` is the public API for the
ordered name list (a `lines(...)` function would also collide with Workstream 5's
`lines::Tuple` keyword).
"""
mutable struct CompositeModel
    lines::Vector{String}               # unique, order = creation order; lines[1] is the reference line
    models::Dict{String,model}
    lineCenters::Dict{String,Float64}   # central wavelength, any consistent units
    fluxRatios::Dict{String,Float64}    # integrated line flux relative to lines[1]; fluxRatios[lines[1]] == 1.0
end

#--- shared validation helpers (used by the constructor and both addLine! methods) ---

_assert_positive(name::Symbol, x::Real) = x > 0 || error("CompositeModel: `$name` must be positive (got $x)")

function _assert_new_line!(cm::CompositeModel, line::String)
    line in cm.lines && error("CompositeModel: line \"$line\" already exists (known lines: $(cm.lines)) -- " *
        "line names are unique and case-sensitive.")
end

"""
    CompositeModel(m::model; line::String, lineCenter::Float64)

Wrap an existing [`model`](@ref) as the first (reference) line of a new `CompositeModel`.
The reference line's `fluxRatio` is fixed at `1.0` -- every subsequent line's
`fluxRatio` (see [`addLine!`](@ref)) is defined relative to it.

# Parameters
- `m::model`: the model for this line (disk-wind, cloud, or a `+`-combined/raytraced model).
- `line::String`: unique name for this line (e.g. `"Hα"`).
- `lineCenter::Float64`: central wavelength (any units, consistent across lines you plan
  to add to this `CompositeModel`). Must be positive.

# Returns
- `CompositeModel` with one registered line.
"""
function CompositeModel(m::model; line::String, lineCenter::Float64)
    _assert_positive(:lineCenter, lineCenter)
    return CompositeModel([line], Dict(line => m), Dict(line => lineCenter), Dict(line => 1.0))
end

"""
    Base.getindex(cm::CompositeModel, line::String) -> model

Retrieve the [`model`](@ref) registered for `line`. Errors (listing the known line
names) if `line` is not present.
"""
function Base.getindex(cm::CompositeModel, line::String)
    haskey(cm.models, line) || error("CompositeModel: unknown line \"$line\" (known lines: $(cm.lines))")
    return cm.models[line]
end

"""
    Base.length(cm::CompositeModel) -> Int

Number of lines registered in `cm`.
"""
Base.length(cm::CompositeModel) = length(cm.lines)

"""
    Base.show(io::IO, cm::CompositeModel)

Print one row per line: name, `lineCenter`, `fluxRatio`, number of rings, and whether
profiles have been computed (`setProfile!`/`getProfile`) for that line's model.
"""
function Base.show(io::IO, cm::CompositeModel)
    println(io, "CompositeModel with $(length(cm)) line(s):")
    for line in cm.lines
        m = cm.models[line]
        hasProfiles = !isnothing(m.profiles) && length(m.profiles) > 0
        println(io, "\t- $line: lineCenter=$(cm.lineCenters[line]), fluxRatio=$(cm.fluxRatios[line]), " *
            "nrings=$(length(m.rings)), profiles set: $hasProfiles")
    end
end

#collect leaf-level params NamedTuples out of a (possibly nested :+/:raytrace!) params tree,
#mirroring rebuild's internal _rebuild_leaf_keys but returning the leaf records themselves --
#used below to find :cloudModel leaves recorded with rng=:legacy for the addLine! warning.
_composite_leaves(::Nothing) = NamedTuple[]
function _composite_leaves(params::NamedTuple)
    c = params.constructor
    if c === :+
        return vcat(_composite_leaves(params.left), _composite_leaves(params.right))
    elseif c === :raytrace!
        return _composite_leaves(params.parent)
    else
        return [params]
    end
end

"""
    addLine!(cm::CompositeModel, m::model; line::String, lineCenter::Float64, fluxRatio::Float64=1.0) -> cm

Register an independently built [`model`](@ref) (which may itself be a `+`-combined,
raytraced multi-component model) as a new line in `cm`. Returns `cm`.

# Parameters
- `line::String`: unique (case-sensitive) name for this line.
- `lineCenter::Float64`: central wavelength (same units convention as `cm`'s other lines). Must be positive.
- `fluxRatio::Float64=1.0`: this line's velocity-integrated flux relative to `cm.lines[1]`. Must be positive.
"""
function addLine!(cm::CompositeModel, m::model; line::String, lineCenter::Float64, fluxRatio::Float64=1.0)
    _assert_new_line!(cm, line)
    _assert_positive(:lineCenter, lineCenter)
    _assert_positive(:fluxRatio, fluxRatio)
    push!(cm.lines, line)
    cm.models[line] = m
    cm.lineCenters[line] = lineCenter
    cm.fluxRatios[line] = fluxRatio
    return cm
end

"""
    addLine!(cm::CompositeModel; line::String, lineCenter::Float64, fluxRatio::Float64=1.0,
             from::String=cm.lines[1], overrides...) -> cm

Add a new line by **parameter reuse**: [`rebuild`](@ref) the recorded construction
parameters (`model.params`, see the `model` docstring and Workstream 4-T1) of the
`from` line, merged with `overrides`, and register the result under `line`. Because
`rebuild` handles nested `:+`/`:raytrace!` records, this works for combined and
raytraced `from` lines too -- `overrides` broadcast to **every** leaf submodel that
has a matching parameter (so e.g. overriding `τ` on a raytraced disk+cloud line
reaches both submodels' intensity functions, if both accept `τ`).

Errors with a clear message if the `from` model has no recorded params
(`params === nothing`, i.e. it was not built by a public constructor) -- pass an
explicit `model` via the other `addLine!` method instead.

# Override semantics -- constructor semantics, no geometry invariance
`rebuild` does exactly what calling the original public constructor with the merged
arguments would do, nothing more -- **no geometry invariance is promised**. Whether an
override moves geometry depends on which parameterization was recorded:
- `:DiskWindModelMean` (i.e. `DiskWindModel(r̄, rFac, α, i; ...)`): `α` (and `r̄`, `rFac`)
  set `rMin`/`rMax` via `get_rMinMaxDiskWind`, so overriding any of them **moves the
  whole grid**.
- `:DiskWindModel` (i.e. explicit `DiskWindModel(rMin, rMax, i; ...)`): `α` reaches only
  the intensity function via `kwargs` -- geometry (`r`, `ϕ`, `v`) is untouched, only `I` changes.

# Cloud-model seed rule
The stored `rng`/`seed` live at the **leaves** of a (possibly nested) params tree, so
this rule is applied per leaf, recursively:
- `rng=:philox` and the caller does **not** override `seed` -> the stored seed is
  reused automatically (unmentioned override keys are left as recorded), producing
  **identical** clouds.
- `rng=:philox` and the caller overrides `seed` -> fresh gas from the new seed.
  Overriding with `seed=nothing` is rejected with a clear error (`cloudModel(...;
  rng=:philox)` requires an explicit integer seed) -- philox cannot draw without one.
- `rng=:legacy` (the model consumed `Random.GLOBAL_RNG` at construction time) -> a
  single `@warn` is emitted noting that gas geometry cannot be reproduced exactly and
  will be re-drawn (statistically equivalent, not identical) from `GLOBAL_RNG`.

# Not recorded
`params` captures construction history only. In-place mutators (`removeNaN!`,
`zeroDiskObscuredClouds!`) are not recorded and are not re-applied by `rebuild`/
`addLine!` -- see the [`rebuild`](@ref) docstring's "Not recorded" section. A caller
who ran such mutators on the `from` model must re-apply the same calls, in the same
order, to the newly added line's model.

# Parameters
- `line::String`: unique (case-sensitive) name for the new line.
- `lineCenter::Float64`: central wavelength. Must be positive.
- `fluxRatio::Float64=1.0`: this line's flux relative to `cm.lines[1]`. Must be positive.
- `from::String=cm.lines[1]`: which registered line's params to reuse.
- `overrides...`: forwarded to `rebuild` (constructor-argument overrides).
"""
function addLine!(cm::CompositeModel; line::String, lineCenter::Float64, fluxRatio::Float64=1.0,
        from::String=cm.lines[1], overrides...)
    fromModel = cm[from] #errors with the known-lines message on a miss
    fromModel.params === nothing && error("addLine!: the model for line \"$from\" was not built by a " *
        "public constructor (params === nothing) -- pass an explicit model instead.")
    leaves = _composite_leaves(fromModel.params)
    if any(l -> l.constructor === :cloudModel && l.rng === :legacy, leaves)
        @warn "addLine!: reusing parameters from a rng=:legacy cloud model (line \"$from\") -- exact gas " *
              "draws cannot be reproduced (GLOBAL_RNG state was consumed at construction); the new line's " *
              "gas geometry will be re-drawn as statistically-equivalent, different gas from GLOBAL_RNG."
    end
    m = rebuild(fromModel.params; overrides...)
    return addLine!(cm, m; line=line, lineCenter=lineCenter, fluxRatio=fluxRatio)
end
