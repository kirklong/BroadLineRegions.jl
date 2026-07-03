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

#=============================== W4-T4: per-line forwarding API ===============================#
# Thin wrappers: no logic changes inside the wrapped functions -- each forwards to the per-line
# `model` method after resolving `line`. The wrapped functions (getProfile, getVariable, setProfile!,
# reset!, removeNaN!, raytrace!) all live earlier in the load order (util.jl/profiles.jl) or are
# generic functions extended later (raytrace.jl); these bodies resolve them at call time.

"""
    getProfile(cm::CompositeModel, name; line=nothing, lines=nothing, kwargs...) -> profile

Per-line forwarding of [`getProfile`](@ref): compute `name`'s profile for the model registered under
`line` (`getProfile(cm[line], name; kwargs...)`).

Julia cannot dispatch on keywords, so this is a *single* method covering the whole positional signature
`getProfile(cm, name; ...)` -- branching happens at runtime:
- `name === :ratio` (velocity-resolved line ratio, which needs `lines::Tuple`) is **implemented in
  Workstream 5** and currently throws a clear error. This same method will be filled in there; two
  keyword-only "methods" would silently overwrite each other.
- every other `name` requires the `line` keyword and forwards to the per-line `model` method.
"""
function getProfile(cm::CompositeModel, name::Union{String,Symbol,Function};
        line::Union{Nothing,String}=nothing, lines=nothing, kwargs...)
    if name === :ratio
        error("getProfile(cm, :ratio; lines=...): the :ratio (velocity-resolved line-ratio) profile is " *
              "implemented in Workstream 5 (Balmer decrement) and is not available yet.")
    end
    line === nothing && error("getProfile(cm, $(repr(name)); line=...): the `line` keyword is required " *
        "(which line's profile do you want?). Known lines: $(cm.lines).")
    return getProfile(cm[line], name; kwargs...)
end

"""
    getVariable(cm::CompositeModel, variable; line::String, flatten=false)

Per-line forwarding of [`getVariable`](@ref): gather `variable` from the model registered under `line`.
"""
function getVariable(cm::CompositeModel, variable::Union{String,Symbol,Function}; line::String, flatten::Bool=false)
    return getVariable(cm[line], variable; flatten=flatten)
end

"""
    setProfile!(cm::CompositeModel, p::profile; line::String, overwrite=false) -> cm

Per-line forwarding of [`setProfile!`](@ref): store profile `p` on the model registered under `line`.
"""
function setProfile!(cm::CompositeModel, p::profile; line::String, overwrite::Bool=false)
    setProfile!(cm[line], p; overwrite=overwrite)
    return cm
end

"""
    reset!(cm::CompositeModel; line=nothing, profiles=true, img=false) -> cm

Per-line forwarding of [`reset!`](@ref). `line=nothing` (the default) resets every line; otherwise only
the named line. Extra keywords (`profiles`, `img`) forward to the per-line `model` method.
"""
function reset!(cm::CompositeModel; line::Union{Nothing,String}=nothing, kwargs...)
    for l in (line === nothing ? cm.lines : [line])
        reset!(cm[l]; kwargs...)
    end
    return cm
end

"""
    removeNaN!(cm::CompositeModel) -> cm

Per-line forwarding of [`removeNaN!`](@ref): drop NaN-sentinel points from *every* line's model.
"""
function removeNaN!(cm::CompositeModel)
    for l in cm.lines
        removeNaN!(cm[l])
    end
    return cm
end

"""
    raytrace!(cm::CompositeModel; line=nothing, kwargs...) -> cm

Per-line forwarding of [`raytrace!`](@ref). `line=nothing` (the default) raytraces every line; otherwise
only the named line.

`raytrace!` **returns a new model** (it does not mutate its argument), so this reassigns
`cm.models[line] = raytrace!(cm.models[line]; kwargs...)`. Lines that `raytrace!` would warn about and
hand back unaltered -- a single-submodel model (`subModelStartInds == [1]`) or an already-raytraced one
(`camera.raytraced`) -- are skipped silently so a multi-line loop does not warn-spam.

Raytracing happens *within* each line's own submodels only. Cross-line occlusion is deliberately **not**
modeled: different lines are at different wavelengths and each line's `τ` is that line's own optical
depth (see the Workstream-4 plan, "Out of scope").
"""
function raytrace!(cm::CompositeModel; line::Union{Nothing,String}=nothing, kwargs...)
    for l in (line === nothing ? cm.lines : [line])
        m = cm[l]
        (m.subModelStartInds == [1] || m.camera.raytraced) && continue #raytrace!'s warn-and-return cases
        cm.models[l] = raytrace!(m; kwargs...)
    end
    return cm
end

"""
    _fluxWeights(cm::CompositeModel) -> Dict{String,Float64}

Internal helper: per-line multiplicative weight `fluxRatios[line] / totalFlux(line)` that renormalizes
each line's binned profile to unit integral before scaling by its `fluxRatio` -- the integrated-flux
semantics of decision 3. `totalFlux(line) = Σ I[k]*ΔA[k]` over points where
`isfinite(v[k]) && isfinite(I[k]*ΔA[k])` -- **both** conditions, matching exactly what `binAccumulate!`
counts (finite-`I`/NaN-position points occur in practice, so finite-`I` alone is not enough).

Uses the memoized flattened gathers (`getVariable(m, …; flatten=true)` hits `model.cache`) and one
non-allocating loop. Deliberately **not** memoized itself -- it recomputes per call so nothing new joins
the `model.cache` invalidation contract; the O(N) sum is noise next to the binning it accompanies.
"""
function _fluxWeights(cm::CompositeModel)
    weights = Dict{String,Float64}()
    for line in cm.lines
        m = cm.models[line]
        v = getVariable(m, :v, flatten=true)
        I = getVariable(m, :I, flatten=true)
        ΔA = getVariable(m, :ΔA, flatten=true)
        total = 0.0
        @inbounds for k in eachindex(v, I, ΔA)
            iΔA = I[k]*ΔA[k]
            if isfinite(v[k]) && isfinite(iΔA)
                total += iΔA
            end
        end
        weights[line] = cm.fluxRatios[line] / total
    end
    return weights
end

#=========================== W4-T5: wavelength utilities + overlap check ===========================#

"""
    wavelength(v, lineCenter) -> λ

First-order velocity→wavelength map `λ = lineCenter * (1 + v)` (decision 4), with `v` the stored
line-of-sight velocity in units of c. Works for scalar or vector `v` (broadcasts). Stored `v` is
**redshift-positive** (positive v = receding), so `(1 + v)` is the correct sign -- do not re-derive it.
"""
wavelength(v, lineCenter) = lineCenter .* (1.0 .+ v)

"""
    _finiteVRange(m::model) -> (vmin, vmax)

Internal helper: **the** definition of a line's velocity range -- the min/max of `v` over points where
`isfinite(v) && isfinite(I*ΔA)` (the W4-T4 finiteness rule; finite-`I` alone would let a NaN `v` poison
the range). Single non-allocating pass over the memoized flattened `v`/`I`/`ΔA` gathers. Shared by
[`lineOverlap`](@ref), [`getSpectrum`](@ref) (and Workstream 5) so the three cannot disagree.

Returns `(NaN, NaN)` when the model has **no** finite points; callers treat that as "this line
contributes no range" (skip it).
"""
function _finiteVRange(m::model)
    v = getVariable(m, :v, flatten=true)
    I = getVariable(m, :I, flatten=true)
    ΔA = getVariable(m, :ΔA, flatten=true)
    vmin = Inf; vmax = -Inf
    @inbounds for k in eachindex(v, I, ΔA)
        if isfinite(v[k]) && isfinite(I[k]*ΔA[k])
            vk = v[k]
            vk < vmin && (vmin = vk)
            vk > vmax && (vmax = vk)
        end
    end
    return vmin > vmax ? (NaN, NaN) : (vmin, vmax) #vmin>vmax (Inf>-Inf) ⟺ no finite points
end

"""
    lineOverlap(cm::CompositeModel) -> Vector{NamedTuple}

Report wavelength-space overlaps between lines. Each line spans
`[lineCenter*(1+vmin), lineCenter*(1+vmax)]` from [`_finiteVRange`](@ref); every pair whose intervals
intersect yields an entry `(lineA, lineB, λlo, λhi)` (the overlap interval `[λlo, λhi]`). An **empty**
vector means no lines overlap. Lines with no finite points (`_finiteVRange` = `NaN`) are skipped.
"""
function lineOverlap(cm::CompositeModel)
    ranges = Dict{String,Tuple{Float64,Float64}}()
    for line in cm.lines
        vmin, vmax = _finiteVRange(cm.models[line])
        λc = cm.lineCenters[line]
        ranges[line] = (wavelength(vmin, λc), wavelength(vmax, λc)) #(NaN,NaN) if no finite points
    end
    result = NamedTuple[]
    n = length(cm.lines)
    for a in 1:n-1, b in a+1:n
        lineA = cm.lines[a]; lineB = cm.lines[b]
        loA, hiA = ranges[lineA]; loB, hiB = ranges[lineB]
        (isnan(loA) || isnan(hiA) || isnan(loB) || isnan(hiB)) && continue
        λlo = max(loA, loB); λhi = min(hiA, hiB)
        if λlo <= λhi #closed intervals intersect
            push!(result, (lineA=lineA, lineB=lineB, λlo=λlo, λhi=λhi))
        end
    end
    return result
end

#================================= W4-T6: combined spectrum =================================#

"""
    getSpectrum(cm::CompositeModel; bins=100, z=0.0, kwargs...)
        -> (edges, centers, flux::Dict{String,Vector{Float64}}, total::Vector{Float64})

Combined wavelength-space spectrum of all lines. Each line's flattened points are binned by
`wavelength(v, lineCenter)*(1+z)` weighted by `I*ΔA*_fluxWeights(cm)[line]`, over a **shared**
wavelength grid spanning every line's finite range (from [`_finiteVRange`](@ref)) times `(1+z)`.
`total` is the elementwise sum of the per-line vectors. With the unit-integral weighting, `sum(flux[line])`
equals `fluxRatios[line]` (decision 3).

# Arguments
- `bins::Union{Int,Vector{Float64}}=100`: bin count, or an explicit edge vector (passthrough as in
  [`binnedSum`](@ref)). For an integer count no edge vector is materialized -- the same `minX`/`maxX` go
  to every line's `binnedSum`, whose deterministic `constructBinEdges` yields bit-identical edges and
  preserves the uniform direct-index fast path.
- `z::Real=0.0`: redshift; `z=0` = rest frame, `z>0` shifts to the observed frame.
- `overflow` (via `kwargs`): honored **only** when `bins` is a user-supplied edge vector (default
  `false`, so out-of-range flux drops against a data grid). For internally-constructed integer-`bins`
  edges `overflow` is forced `true` -- non-centered edges built at exactly `[λmin, λmax]` place each
  line's extremal points on the boundary, which `binnedSum` otherwise drops (deliberate package-wide
  convention); `overflow=true` is the sanctioned mechanism to keep them and is what makes the integrated
  flux exact.
"""
function getSpectrum(cm::CompositeModel; bins::Union{Int,Vector{Float64}}=100, z::Real=0.0, kwargs...)
    #shared wavelength range across all lines (observed frame), from the per-line finite velocity ranges
    λmin = Inf; λmax = -Inf
    for line in cm.lines
        vmin, vmax = _finiteVRange(cm.models[line])
        (isnan(vmin) || isnan(vmax)) && continue
        λc = cm.lineCenters[line]
        lo = wavelength(vmin, λc)*(1.0+z); hi = wavelength(vmax, λc)*(1.0+z)
        lo < λmin && (λmin = lo)
        hi > λmax && (λmax = hi)
    end
    (isfinite(λmin) && isfinite(λmax)) || error("getSpectrum: no line has any finite points -- cannot " *
        "build a spectrum (every line's I*ΔA / v is NaN everywhere).")

    #overflow: forced true for internally-constructed integer-bins edges (keeps boundary points);
    #for a user edge vector, pass through untouched (default false -- out-of-range flux drops).
    overflow = (typeof(bins) == Int) ? true : get(kwargs, :overflow, false)

    weights = _fluxWeights(cm)
    flux = Dict{String,Vector{Float64}}()
    edges = nothing; centers = nothing
    for line in cm.lines
        m = cm.models[line]
        v = getVariable(m, :v, flatten=true)
        I = getVariable(m, :I, flatten=true)
        ΔA = getVariable(m, :ΔA, flatten=true)
        λc = cm.lineCenters[line]; w = weights[line]
        x = wavelength(v, λc) .* (1.0+z)
        y = (I .* ΔA) .* w
        #integer bins -> pass minX/maxX (no edge vector); constructBinEdges is deterministic so every
        #line gets bit-identical edges. A user edge vector ignores minX/maxX inside binnedSum.
        e, c, r = binnedSum(x, y; bins=bins, minX=λmin, maxX=λmax, centered=false, overflow=overflow)
        if edges === nothing
            edges = e; centers = c
        end
        flux[line] = r
    end
    total = zeros(length(centers))
    for line in cm.lines
        total .+= flux[line]
    end
    return (edges, centers, flux, total)
end
