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
- `name === :ratio` computes the **velocity-resolved line ratio** (below); it requires the
  `lines::Tuple{String,String}` keyword instead of `line`.
- every other `name` requires the `line` keyword and forwards to the per-line `model` method.

# `:ratio` -- velocity-resolved line ratio (Workstream 5)

    getProfile(cm, :ratio; lines=(a, b), bins=100, floor=0.0,
               minX=nothing, maxX=nothing, centered=false, overflow=nothing) -> profile

The per-velocity-bin flux ratio of line `a` (numerator) to line `b` (denominator) on **shared**
velocity bins -- e.g. the velocity-resolved Balmer decrement BD(v) = F_Hα(v)/F_Hβ(v) of
Chen et al. 2026 (arXiv:2606.04711 §5.3) with `lines=("Hα", "Hβ")`. Each line's flattened points are
binned by `v` weighted by `I*ΔA*_fluxWeights(cm)[line]` (the same integrated-flux weighting as
[`getSpectrum`](@ref), so the ratio of the two lines' *integrated* fluxes is pinned to
`fluxRatios[a]/fluxRatios[b]` = [`lineRatio`](@ref)`(cm, a, b)`; the per-line radial emissivity
profiles set only the velocity *shape*). Two lines with identical realized geometry and intensity
have a flat ratio equal to that constant.

The returned [`profile`](@ref) has `name = Symbol(a, "/", b)` (e.g. `Symbol("Hα/Hβ")`), so it can be
stored via `setProfile!(cm, p; line=...)` and plotted via the `profile`-struct recipe. Bins where the
denominator is empty (`0.0`), or below `floor` times the finite-only maximum of the denominator's
binned flux (an SNR-like guard for sparse cloud models, off by default), are set to `NaN` -- never
`Inf` -- per the package NaN-sentinel convention.

## `:ratio` keywords
- `lines::Tuple{String,String}`: `(numerator, denominator)` line names. Required.
- `bins::Union{Int,Vector{Float64}}=100`: bin count or explicit edge vector, as in [`binnedSum`](@ref).
  Both lines always share bit-identical edges.
- `floor::Real=0.0`: denominator threshold as a fraction of the denominator's largest binned flux;
  bins below it become `NaN`. `0.0` disables the guard (only empty bins are `NaN`).
- `minX`/`maxX::Union{Real,Nothing}=nothing`: pin one or both ends of the shared velocity grid instead
  of auto-computing them as the **union** of the two lines' finite ranges ([`_finiteVRange`](@ref)).
  Ignored when `bins` is an edge vector.
- `centered::Bool=false` and `overflow::Union{Bool,Nothing}=nothing`: forwarded to each line's
  [`binnedSum`](@ref) with exactly the [`getSpectrum`](@ref) semantics -- the `overflow` default
  resolves to `true` for integer `bins` (keeps boundary/out-of-window points, making the two binned
  totals exact) and `false` for a user edge vector (out-of-range model flux drops against a data grid).

## Comparing with data
- **Data-side velocity grids**: pass the measurement's bin edges as `bins` (a `Vector{Float64}`).
  Stored `v` is in units of c, so convert km/s edges first: `edges_c = edges_kms ./ 2.99792458e5`.
- **Velocity shift**: there is deliberately no `vShift` keyword -- shifting the model by `v_shift` is
  exactly evaluating it on shifted edges, so subtract it from the data grid instead:
  `bins = edges_kms ./ 2.99792458e5 .- v_shift_c`. (With auto-constructed edges a common shift is a
  no-op on the ratio: both lines and the edges move together.)
- **Instrumental broadening**: LSF convolution stays outside the package. `:ratio` divides the
  *unconvolved* binned profiles; to match LSF-matched data (the paper convolves the narrower line to
  the broader one's LSF before dividing), bin each line separately, convolve, and divide yourself.
"""
function getProfile(cm::CompositeModel, name::Union{String,Symbol,Function};
        line::Union{Nothing,String}=nothing, lines=nothing, kwargs...)
    if name === :ratio
        line === nothing || error("getProfile(cm, :ratio; ...): use `lines=(numerator, denominator)`, " *
            "not `line` -- the :ratio profile is computed from a PAIR of lines.")
        lines isa Tuple{String,String} || error("getProfile(cm, :ratio; lines=...): the `lines` keyword " *
            "is required and must be a Tuple{String,String} of (numerator, denominator) line names " *
            "(got $(repr(lines))). Known lines: $(cm.lines).")
        return _ratioProfile(cm, lines; kwargs...)
    end
    lines === nothing || error("getProfile(cm, $(repr(name)); ...): the `lines` keyword is only used by " *
        "the :ratio profile -- pass `line` to select a single line's profile.")
    line === nothing && error("getProfile(cm, $(repr(name)); line=...): the `line` keyword is required " *
        "(which line's profile do you want?). Known lines: $(cm.lines).")
    return getProfile(cm[line], name; kwargs...)
end

#W5-T1 implementation of the :ratio branch above (all user-facing documentation lives on getProfile).
#Mirrors getSpectrum's structure: per-side range fill from _finiteVRange, the same overflow default
#rule, integer `bins` + shared minX/maxX (bit-identical deterministic edges, direct-index fast path)
#or a user edge vector passed straight through.
function _ratioProfile(cm::CompositeModel, lines::Tuple{String,String};
        bins::Union{Int,Vector{Float64}}=100, floor::Real=0.0,
        minX::Union{Real,Nothing}=nothing, maxX::Union{Real,Nothing}=nothing,
        centered::Bool=false, overflow::Union{Bool,Nothing}=nothing, kwargs...)
    a, b = lines
    mA = cm[a]; mB = cm[b] #errors listing the known lines on a miss
    floor >= 0.0 || error("getProfile(cm, :ratio; ...): `floor` must be >= 0 (got $floor)")
    #shared velocity range: any side not pinned by the caller via minX/maxX is auto-computed as the
    #union of the two lines' finite velocity ranges (both pinned -> skip the pass)
    vmin = minX === nothing ? Inf : Float64(minX)
    vmax = maxX === nothing ? -Inf : Float64(maxX)
    if minX === nothing || maxX === nothing
        for m in (mA, mB)
            lo, hi = _finiteVRange(m)
            (isnan(lo) || isnan(hi)) && continue
            minX === nothing && lo < vmin && (vmin = lo)
            maxX === nothing && hi > vmax && (vmax = hi)
        end
    end
    (isfinite(vmin) && isfinite(vmax)) || error("getProfile(cm, :ratio; lines=$(repr(lines))): neither " *
        "line has any finite points -- cannot build shared velocity bins.")

    #overflow default: true for integer-bins edges (keeps boundary/out-of-window points -- the two
    #binned totals stay exact), false for a user edge vector (out-of-range flux drops against a data
    #grid). An explicit user overflow always wins. Same rule as getSpectrum.
    overflow = overflow === nothing ? (bins isa Int) : overflow

    weights = _fluxWeights(cm)
    edges = nothing; centers = nothing; num = nothing; den = nothing
    for line in (a, b)
        m = cm.models[line]
        v = getVariable(m, :v, flatten=true)
        I = getVariable(m, :I, flatten=true)
        ΔA = getVariable(m, :ΔA, flatten=true)
        y = (I .* ΔA) .* weights[line]
        #integer bins -> pass minX/maxX (no edge vector); constructBinEdges is deterministic so both
        #lines get bit-identical edges. A user edge vector ignores minX/maxX inside binnedSum.
        e, c, s = binnedSum(v, y; bins=bins, minX=vmin, maxX=vmax, centered=centered, overflow=overflow)
        if num === nothing
            edges = e; centers = c; num = s
        else
            den = s
        end
    end

    return profile(name=Symbol(a, "/", b), binCenters=centers, binEdges=edges, binSums=_ratioDivide(num, den, floor))
end

#shared num/den division for the :ratio profile -- used by both the host path above and the resident
#path (gpu_observables.jl, W5-G1; num/den are small host vectors on both), so the NaN/floor semantics
#cannot drift between them. Empty (0.0) or below-floor denominator bins become NaN, never Inf.
function _ratioDivide(num::Vector{Float64}, den::Vector{Float64}, floor::Real)
    #NaN-proof floor threshold: `floor` times a finite-only maximum of den. den comes straight from
    #the binned sums (zero-initialized, accumulate only finite values) so it is in fact always finite,
    #but the finite-only reduction costs nothing and keeps the guard correct if a stored/derived
    #vector is ever fed through this path.
    thresh = 0.0
    if floor > 0.0
        denMax = -Inf
        @inbounds for d in den
            isfinite(d) && d > denMax && (denMax = d)
        end
        thresh = isfinite(denMax) ? floor*denMax : 0.0
    end
    binSums = similar(den)
    @inbounds for k in eachindex(num, den, binSums)
        dk = den[k]
        binSums[k] = (dk == 0.0 || (floor > 0.0 && dk < thresh)) ? NaN : num[k]/dk
    end
    return binSums
end

"""
    lineRatio(cm::CompositeModel, a::String, b::String) -> Float64

**Integrated** flux ratio of line `a` to line `b`. By the integrated-flux semantics of `fluxRatio`
(each line's profile is normalized to unit integral before scaling, see [`_fluxWeights`](@ref)) this
is exactly `cm.fluxRatios[a] / cm.fluxRatios[b]` -- the constant that the velocity-resolved
`getProfile(cm, :ratio; lines=(a, b))` profile integrates to (and equals bin-by-bin in the degenerate
identical-geometry case). Errors (listing the known lines) if either name is not registered.

A [`ResidentCompositeModel`](@ref) is supported through the same shared implementation (the
`::ResidentCompositeModel` method lives with its resident siblings in `gpu_observables.jl`, W5-G1) --
`fluxRatios` is host metadata on both, so the two methods are identical.
"""
function lineRatio(cm::CompositeModel, a::String, b::String)
    return _lineRatio(cm, a, b)
end

#shared lineRatio implementation -- deliberately UNtyped, like `_lineOverlap`: the body only touches
#the host `fluxRatios`/`lines` metadata, which host and resident composites carry identically.
function _lineRatio(cm, a::String, b::String)
    for l in (a, b)
        haskey(cm.fluxRatios, l) || error("lineRatio: unknown line \"$l\" (known lines: $(cm.lines))")
    end
    return cm.fluxRatios[a] / cm.fluxRatios[b]
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

Errors (naming the offending line) if any line's total flux is not positive: a zero-emission line
cannot be normalized to unit integral, and the alternative -- an `Inf` weight whose non-finite
weighted fluxes the binning silently drops -- would make the line quietly vanish from the spectrum.
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
        #a zero (or negative) total cannot be normalized to unit integral -- error loudly by name
        #rather than produce an Inf weight, whose non-finite weighted fluxes binAccumulate! would
        #silently drop (the line would just vanish from the spectrum, breaking sum(flux) == fluxRatio)
        total > 0.0 || error("_fluxWeights: line \"$line\" has no positive integrated flux " *
            "(Σ I*ΔA = $total over its finite points) -- cannot normalize its profile to unit " *
            "integral. Check the line's intensity function/mask; a line with no emission cannot " *
            "carry a fluxRatio.")
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

A [`ResidentCompositeModel`](@ref) is supported through the same shared implementation (the
`::ResidentCompositeModel` method lives with its resident siblings in `gpu_observables.jl`, W4-G3);
the only difference is that each line's velocity range then comes from the device reductions of
`_finiteVRange(::ResidentModel)` instead of the host gathers.
"""
function lineOverlap(cm::CompositeModel)
    return _lineOverlap(cm)
end

#shared lineOverlap implementation -- deliberately UNtyped: the body is representation-generic (it
#only touches `lines`/`models`/`lineCenters` and calls `_finiteVRange`, which dispatches on the
#per-line model type, host `model` or `ResidentModel`), so the host and resident public methods
#cannot disagree.
function _lineOverlap(cm)
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
- `minX`/`maxX::Union{Real,Nothing}=nothing`: pin one or both ends of the shared wavelength grid
  (observed frame, i.e. after the `(1+z)` shift) instead of auto-computing them from the lines' finite
  ranges -- forwarded to every line's [`binnedSum`](@ref) like the other binning keywords the simple
  `model` methods forward. Ignored when `bins` is an edge vector (the edges pin the grid).
- `centered::Bool=false`: forwarded to `binnedSum`/`constructBinEdges`. The default `false` builds
  edges at exactly `[minX, maxX]`; `centered=true` pads them half a bin so the extremes fall in bin
  centers, as in the `getProfile` default.
- `overflow::Union{Bool,Nothing}=nothing`: an explicit `true`/`false` is forwarded to `binnedSum`
  untouched. The default (`nothing`) resolves to `true` for integer `bins` and `false` for a
  user-supplied edge vector: non-centered internal edges built at exactly `[λmin, λmax]` place each
  line's extremal points on the boundary, which `binnedSum` otherwise drops (deliberate package-wide
  convention) -- `overflow=true` is the sanctioned mechanism to keep them, is what makes the integrated
  flux exact, and (for a user-narrowed `minX`/`maxX` window) keeps out-of-window flux by accumulating
  it into the boundary bins. Pass `overflow=false` explicitly for drop-outside-the-window semantics;
  against a user edge vector (a data grid) out-of-range flux drops by default.
"""
function getSpectrum(cm::CompositeModel; bins::Union{Int,Vector{Float64}}=100, z::Real=0.0,
        minX::Union{Real,Nothing}=nothing, maxX::Union{Real,Nothing}=nothing,
        centered::Bool=false, overflow::Union{Bool,Nothing}=nothing, kwargs...)
    #shared wavelength range across all lines (observed frame): any side not pinned by the caller via
    #minX/maxX is auto-computed from the per-line finite velocity ranges (both pinned -> skip the pass)
    λmin = minX === nothing ? Inf : Float64(minX)
    λmax = maxX === nothing ? -Inf : Float64(maxX)
    if minX === nothing || maxX === nothing
        for line in cm.lines
            vmin, vmax = _finiteVRange(cm.models[line])
            (isnan(vmin) || isnan(vmax)) && continue
            λc = cm.lineCenters[line]
            lo = wavelength(vmin, λc)*(1.0+z); hi = wavelength(vmax, λc)*(1.0+z)
            minX === nothing && lo < λmin && (λmin = lo)
            maxX === nothing && hi > λmax && (λmax = hi)
        end
    end
    (isfinite(λmin) && isfinite(λmax)) || error("getSpectrum: no line has any finite points -- cannot " *
        "build a spectrum (every line's I*ΔA / v is NaN everywhere).")

    #overflow default: true for integer-bins edges (keeps boundary/out-of-window points -- preserves the
    #integrated-flux identity), false for a user edge vector (out-of-range flux drops against a data
    #grid). An explicit user overflow always wins.
    overflow = overflow === nothing ? (bins isa Int) : overflow

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
        e, c, r = binnedSum(x, y; bins=bins, minX=λmin, maxX=λmax, centered=centered, overflow=overflow)
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

#================================= W4-T7: visualization =================================#
# Include-order note: src/BroadLineRegions.jl includes util.jl (position 2) long before composite.jl
# (position after operators.jl) -- any top-level method *signature* annotated `::CompositeModel` would
# UndefVarError there, since the type doesn't exist yet at parse time. `@userplot Profile`'s CompositeModel
# support therefore lives as a runtime `isa CompositeModel` branch INSIDE that recipe's body (util.jl,
# fine because the name resolves when the recipe is actually applied, long after the module has finished
# loading) -- see the `profile`/`Profile` docstring there. Everything below needs a `::CompositeModel`
# signature (or calls composite-only helpers like `getSpectrum`/`lineOverlap`), so it lives here instead.
# `Plots`/`RecipesBase` (and their `@recipe`/`@userplot` macros) are already in scope module-wide from
# util.jl's `using Plots, RecipesBase` (position 2) by the time this file is included -- no include-order
# change needed for that part.

"""
    spectrum(cm::CompositeModel; z=0.0, bins=100, kwargs...)
    spectrum(rcm::ResidentCompositeModel; z=0.0, bins=100, kwargs...)

Combined wavelength-space spectrum plot: one series per line plus a black `"total"` series (both from
[`getSpectrum`](@ref)), with pairwise line overlaps ([`lineOverlap`](@ref)) shaded as gray bands behind
the curves.

A [`ResidentCompositeModel`](@ref) works through the same recipe body with no extra branch (W4-G3):
`getSpectrum` and `lineOverlap` resolve at call time and dispatch to the resident methods
(`gpu_observables.jl`), which run the per-line reductions/binning on the device. The resident
`getSpectrum` restrictions apply (a user-supplied `bins` edge vector must be uniform).

# Keywords
- `z::Real=0.0`: redshift, forwarded to `getSpectrum` (`z=0` = rest frame; matches `getSpectrum`'s
  observed-frame convention).
- `bins::Union{Int,Vector{Float64}}=100`: bin count or explicit edge vector, forwarded to `getSpectrum`.
- other `kwargs...`: ordinary `Plots` attributes.

`lineOverlap` reports rest-frame wavelength intervals; they are scaled by `(1+z)` here before shading so
the bands line up with the (possibly redshifted) spectrum.
"""
@userplot Spectrum
@recipe function f(s::Spectrum; z=0.0, bins=100)
    length(s.args) == 1 || error("spectrum(cm): expected a single CompositeModel argument, got $(s.args)")
    cm = s.args[1]
    edges, centers, flux, total = getSpectrum(cm; bins=bins, z=z)
    title --> "Composite spectrum"
    xlabel --> "λ [input units]"
    ylabel --> "flux [arb.]"
    legend --> :topright
    overlaps = lineOverlap(cm)
    if !isempty(overlaps)
        ymin = min(0.0, minimum(total))
        ymax = 1.05*maximum(total)
        for ov in overlaps
            lo = ov.λlo*(1.0+z); hi = ov.λhi*(1.0+z) #rest-frame -> observed frame (decision 4/W4-T5)
            @series begin
                subplot := 1
                seriestype := :shape
                fillcolor --> :gray
                fillalpha --> 0.25
                linealpha --> 0.0
                label --> ""
                x := [lo, hi, hi, lo]
                y := [ymin, ymin, ymax, ymax]
                ()
            end
        end
    end
    for (i,line) in enumerate(cm.lines)
        @series begin
            subplot := 1
            seriestype := :path
            color --> i
            label --> line
            x := centers
            y := flux[line]
            ()
        end
    end
    @series begin
        subplot := 1
        seriestype := :path
        color --> :black
        linewidth --> 2
        label --> "total"
        x := centers
        y := total
        ()
    end
end

"""
    image(cm::CompositeModel, [variable]; line::String, kwargs...)

Per-line forwarding of the [`image`](@ref) recipe: renders `cm[line]` exactly as
`image(cm[line], variable; kwargs...)` would. `line` is required (no default) -- unlike [`plot3d`](@ref),
an intensity image has no meaningful "overlay all lines" mode (each line's `I` is in its own arbitrary
units, so mixing them in one color scale would be misleading).
"""
function image(cm::CompositeModel, variable::Union{String,Symbol,Function}=:I; line::String, kwargs...)
    return image(cm[line], variable; kwargs...)
end

"""
    plot3d(cm::CompositeModel, [variable], [annotate]; line=nothing, kwargs...)

Per-line forwarding of the [`plot3d`](@ref) recipe for a [`CompositeModel`](@ref).

- `line::String`: identical to `plot3d(cm[line], variable, annotate; kwargs...)`.
- `line=nothing` (default): overlay **every** line's geometry on one 3D plot, one solid color per line.
  Reuses `_plot3d_submodels` -- the same per-submodel point gather the single-model `plot3d`
  recipe uses -- for each line's `(x,y,z)` (and its NaN mask on `variable`); unlike the single-model
  recipe this does NOT color points by `variable`'s value (a shared colormap across lines with
  independently-scaled intensities/variables would not be meaningful) -- points are colored by line
  instead, and `variable` only selects which column supplies the per-point finite/NaN mask (default
  `geometry`, i.e. no masking, matching the single-model recipe's default).
"""
function plot3d(cm::CompositeModel, args...; line::Union{Nothing,String}=nothing, kwargs...)
    line !== nothing && return plot3d(cm[line], args...; kwargs...)
    return _plot3d_composite(cm, args...; kwargs...)
end

#shared all-lines overlay implementation -- deliberately UNtyped, like `_lineOverlap`: the body only
#touches `lines`/getindex and calls `_plot3d_submodels`, which dispatches on the per-line model type
#(host `model` or `ResidentModel`), so the host and resident overlays cannot disagree. For a resident
#line the `_plot3d_submodels(::ResidentModel, …)` restrictions apply: a multi-submodel line needs
#raytrace metadata (build with `raytrace=true`) and custom-`Function` mask variables are not supported.
function _plot3d_composite(cm, args...; kwargs...)
    variable, annotate = geometry, true
    if length(args) == 1
        args[1] isa Bool ? (annotate = args[1]) : (variable = args[1])
    elseif length(args) == 2
        tmp1, tmp2 = args
        if tmp1 isa Bool
            annotate, variable = tmp1, tmp2
        else
            variable, annotate = tmp1, tmp2
        end
    elseif length(args) > 2
        throw(ArgumentError("plot3d(cm; line=nothing) overlay expects up to 2 positional arguments (variable, annotate), got $(length(args))"))
    end
    perLine = [_plot3d_submodels(cm[l], variable) for l in cm.lines] #(subs, camR) per line -- exactly what the single-model recipe gathers
    finiteMax(v) = maximum(t for t in v if !isnan(t))
    boxSize = maximum(camR for (_,camR) in perLine)
    for (subs,_) in perLine, sub in subs
        boxSize = max(boxSize, 1.1*max(finiteMax(sub.x), finiteMax(sub.y), finiteMax(sub.z)))
    end
    p = Plots.plot(; xlabel="x [rₛ]", ylabel="y [rₛ]", zlabel="z [rₛ]", aspect_ratio=:equal,
        title="System geometry (all lines)", foreground_color_legend=nothing,
        xlims=(-boxSize,boxSize), ylims=(-boxSize,boxSize), zlims=(-boxSize,boxSize), kwargs...)
    for (i,l) in enumerate(cm.lines)
        subs, _ = perLine[i]
        for (k,sub) in enumerate(subs)
            finite = .!isnan.(sub.mz)
            Plots.scatter3d!(p, sub.x[finite], sub.y[finite], sub.z[finite];
                markerstrokewidth=0.0, markersize=1.0, markeralpha=(sub.disk ? 0.9 : 0.1),
                color=i, label=(k == 1 ? l : ""))
        end
    end
    if annotate
        Plots.plot!(p, [0,boxSize], [0,0], [0,0]; seriestype=:path, linecolor=:dodgerblue,
            linewidth=1.0, linestyle=:dash, arrow=:arrow, label="camera")
    end
    return p
end

#===================== W4-G1: device-resident composite (gpu/resident forwarding) =====================#
# Include-order note: `ResidentCompositeModel` (and the per-model `resident`/`gpu` entry points) live in
# src/gpu_arrays.jl, which is included BEFORE this file -- but the `::CompositeModel` method signatures
# below cannot live there, because gpu_arrays.jl is parsed before the CompositeModel type exists (same
# rule as W4-T7's recipe placement). Both types are defined by the time this file is included, so the
# delegating methods go here. `gpu(::model)` itself lives in the CUDA package extension
# (ext/BroadLineRegionsCUDAExt.jl; CUDA is a weakdep) -- the delegating `gpu(cm::CompositeModel)` can
# live in src because it only calls `gpu` per line, so without CUDA loaded the existing `gpu(::Any)`
# fallback error still fires.

"""
    resident(cm::CompositeModel; T=Float64, backend=KernelAbstractions.CPU(), raytrace=false)
        -> ResidentCompositeModel

Flatten every line's model with the per-model [`resident`](@ref) (each line gets the same keyword
arguments) and wrap the per-line [`ResidentModel`](@ref) handles as a
[`ResidentCompositeModel`](@ref). The default `CPU()` backend keeps everything on the host — useful
for testing the resident composite pipeline without a GPU. Use `gpu(cm)` (with CUDA.jl loaded) to
build a device-resident composite.
"""
function resident(cm::CompositeModel; kwargs...)
    models = Dict{String,ResidentModel}(line => resident(cm.models[line]; kwargs...) for line in cm.lines)
    return ResidentCompositeModel(copy(cm.lines), models, copy(cm.lineCenters), copy(cm.fluxRatios))
end

"""
    gpu(cm::CompositeModel; T=Float32, kwargs...) -> ResidentCompositeModel

Move every line's model onto the GPU with the per-model [`gpu`](@ref) (each line gets the same
keyword arguments, e.g. `T=Float32`) and wrap the per-line [`ResidentModel`](@ref) handles as a
[`ResidentCompositeModel`](@ref). Requires CUDA.jl to be loaded so the package extension can provide
the CUDA-backed `gpu(::model)` method — without it the `gpu(::Any)` fallback error fires per line.
"""
function gpu(cm::CompositeModel; kwargs...)
    models = Dict{String,ResidentModel}(line => gpu(cm.models[line]; kwargs...) for line in cm.lines)
    return ResidentCompositeModel(copy(cm.lines), models, copy(cm.lineCenters), copy(cm.fluxRatios))
end

#===================== W4-G3: resident composite recipe forwarding (image/plot3d) =====================#
# These mirror the `::CompositeModel` visualization forwarding above (W4-T7) one-for-one, so they live
# next to those siblings; the per-line calls land on the existing ResidentModel recipe paths (host
# copies of the device columns, see gpu_arrays.jl). The `getProfile`/`lineOverlap` resident-composite
# forwarding lives with ITS resident siblings in gpu_observables.jl (W4-G2/G3), same split as the CPU
# phase: recipes/forwarding here, observables there. The `spectrum` and `profile` recipes need no new
# methods at all -- their bodies resolve `getSpectrum`/`lineOverlap`/`getProfile` at call time, which
# dispatch to the resident methods (see the runtime-branch notes on each recipe).

"""
    image(rcm::ResidentCompositeModel, [variable]; line::String, kwargs...)

Per-line forwarding of the [`image`](@ref) recipe for a [`ResidentCompositeModel`](@ref): renders
`rcm[line]` exactly as `image(rcm[line], variable; kwargs...)` would (the [`ResidentModel`](@ref)
recipe path -- host copies of that line's device columns). As on the host composite, `line` is
required (no default): an intensity image has no meaningful "overlay all lines" mode. `variable`
must be a `ModelArrays` column name (custom `Function` variables need the host `model` path).
"""
function image(rcm::ResidentCompositeModel, variable::Union{String,Symbol,Function}=:I; line::String, kwargs...)
    return image(rcm[line], variable; kwargs...)
end

"""
    plot3d(rcm::ResidentCompositeModel, [variable], [annotate]; line=nothing, kwargs...)

Per-line forwarding of the [`plot3d`](@ref) recipe for a [`ResidentCompositeModel`](@ref), matching
`plot3d(::CompositeModel, …)`:

- `line::String`: identical to `plot3d(rcm[line], variable, annotate; kwargs...)` (the
  [`ResidentModel`](@ref) recipe path).
- `line=nothing` (default): overlay **every** line's geometry on one 3D plot, one solid color per
  line, via the same shared implementation as the host composite overlay. The
  `_plot3d_submodels(::ResidentModel, …)` restrictions apply per line: a multi-submodel line must
  carry raytrace metadata (build the handle with `raytrace=true`), and custom-`Function` mask
  variables are not supported (pass a column `Symbol`, or use the host `model`).
"""
function plot3d(rcm::ResidentCompositeModel, args...; line::Union{Nothing,String}=nothing, kwargs...)
    line !== nothing && return plot3d(rcm[line], args...; kwargs...)
    return _plot3d_composite(rcm, args...; kwargs...)
end
