#!/usr/bin/env julia

"""
    rebuild(params::NamedTuple; overrides...) -> model
    rebuild(m::model; overrides...) -> model

Reconstruct a [`model`](@ref) from the construction provenance recorded in `params` (see the `model`
`params` field), optionally replacing individual constructor arguments via `overrides`.

`rebuild` dispatches on `params.constructor`:

- a leaf record (`:DiskWindModel`, `:DiskWindModelMean`, `:cloudModel`, `:cloudModelVectors`) calls the
  corresponding public constructor with its recorded arguments;
- a `:+` record rebuilds `left` and `right` and sums them;
- a `:raytrace!` record rebuilds `parent` and re-raytraces it with the recorded raytrace arguments.

# Overrides
`overrides...` broadcast to **every leaf** in which the key exists (so an override reaches all submodels
of a combined/raytraced model that share that parameter). An override key that matches **no** leaf is an
error (typo guard). **Override semantics are exactly constructor semantics** -- whether an override moves
geometry depends on the recorded parameterization. For example overriding `α` on a `:DiskWindModelMean`
record moves the whole grid (`α` sets `rMin`/`rMax` via `get_rMinMaxDiskWind`), while overriding `α` on a
`:DiskWindModel` (explicit `rMin`/`rMax`) record only reaches the intensity function -- geometry is
untouched. `rebuild` does exactly what calling the constructor with the merged arguments would do,
nothing more.

Node-level arguments (`IRatios`, `τCutOff`, `raytraceFreeClouds`) are **not** override targets in v1: to
change them, rebuild the parent and re-raytrace manually.

# Reproducibility
Disk-wind construction is deterministic, so a rebuilt disk-wind line has identical geometry. Cloud models
reproduce identical clouds only under `rng=:philox` with the recorded `seed`; a `rng=:legacy` record
consumed the `GLOBAL_RNG` state, so its rebuild re-draws statistically-equivalent (different) gas from
`GLOBAL_RNG`.

# Not recorded
`params` captures **construction** history only. In-place mutators (`removeNaN!`,
`zeroDiskObscuredClouds!`) are not recorded; a caller who ran them must re-apply the same calls in the
same order to the rebuilt model. `removeNaN!` is observable-neutral (binned observables skip non-finite
points), so the only unrecorded mutator with a physical effect is `zeroDiskObscuredClouds!`, which is
deterministic and cheap to re-run.

`rebuild(m::model; ...)` is a convenience method forwarding to `rebuild(m.params; ...)`; it errors when
`m.params === nothing` (the model was not built by a public constructor).
"""
function rebuild(params::NamedTuple; overrides...)
    ov = (; overrides...) #override kwargs as a NamedTuple
    leafKeys = _rebuild_leaf_keys(params)
    for k in keys(ov)
        k in leafKeys || error("rebuild: override key `$k` matches no leaf constructor parameter " *
            "(known leaf keys: $(sort!(collect(leafKeys)))). Node-level raytrace! arguments " *
            "(IRatios, τCutOff, raytraceFreeClouds) are not override targets.")
    end
    return _rebuild(params, ov)
end

function rebuild(m::model; overrides...)
    m.params === nothing && error("rebuild: this model has no recorded `params` (it was not built by a " *
        "public constructor). Pass an explicitly built model instead of relying on parameter reuse.")
    return rebuild(m.params; overrides...)
end

#collect the set of leaf-level override-able keys, recursing through :+ / :raytrace! nodes
_rebuild_leaf_keys(::Nothing) = Set{Symbol}()
function _rebuild_leaf_keys(params::NamedTuple)
    c = params.constructor
    if c === :+
        return union(_rebuild_leaf_keys(params.left), _rebuild_leaf_keys(params.right))
    elseif c === :raytrace!
        return _rebuild_leaf_keys(params.parent) #node-level args deliberately excluded
    else
        return Set(k for k in keys(params) if k !== :constructor)
    end
end

#a nothing branch means one operand of a + / raytrace! was built without a public constructor
_rebuild(::Nothing, ::NamedTuple) = error("rebuild: encountered a `+`/`raytrace!` branch whose operand " *
    "has no recorded params (built without a public constructor); that branch cannot be reconstructed.")

function _rebuild(params::NamedTuple, ov::NamedTuple)
    c = params.constructor
    if c === :+
        return _rebuild(params.left, ov) + _rebuild(params.right, ov)
    elseif c === :raytrace!
        parent = _rebuild(params.parent, ov)
        return raytrace!(parent; IRatios=params.IRatios, τCutOff=params.τCutOff,
            raytraceFreeClouds=params.raytraceFreeClouds)
    else
        return _rebuild_leaf(params, ov)
    end
end

#drop positional/marker keys, returning the remainder as a kwargs NamedTuple
_leaf_kwargs(nt::NamedTuple, drop::Tuple) = (; (k => nt[k] for k in keys(nt) if !(k in drop))...)

function _rebuild_leaf(params::NamedTuple, ov::NamedTuple)
    #apply only the overrides whose key exists in this leaf
    merged = merge(params, (; (k => ov[k] for k in keys(ov) if haskey(params, k))...))
    c = merged.constructor
    if c === :DiskWindModel
        return DiskWindModel(merged.rMin, merged.rMax, merged.i;
            _leaf_kwargs(merged, (:constructor, :rMin, :rMax, :i))...)
    elseif c === :DiskWindModelMean
        return DiskWindModel(merged.r̄, merged.rFac, merged.α, merged.i;
            _leaf_kwargs(merged, (:constructor, :r̄, :rFac, :α, :i))...)
    elseif c === :cloudModelVectors
        return cloudModel(merged.ϕ₀, merged.i, merged.rot, merged.θₒ, merged.θₒSystem, merged.ξ;
            _leaf_kwargs(merged, (:constructor, :ϕ₀, :i, :rot, :θₒ, :θₒSystem, :ξ))...)
    elseif c === :cloudModel
        return _rebuild_cloudModel(merged)
    else
        error("rebuild: unknown constructor `$c` in recorded params")
    end
end

function _rebuild_cloudModel(merged::NamedTuple)
    kw = _leaf_kwargs(merged, (:constructor, :nClouds, :rng, :seed))
    if merged.rng === :philox
        return cloudModel(merged.nClouds; rng=:philox, seed=merged.seed, kw...)
    else
        #legacy record: GLOBAL_RNG state was consumed, so exact draws cannot be reproduced --
        #re-draw statistically-equivalent gas from the default GLOBAL_RNG.
        return cloudModel(merged.nClouds; kw...)
    end
end
