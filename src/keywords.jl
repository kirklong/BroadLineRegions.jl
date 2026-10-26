# Keyword names are normalized before defaults or any model work is evaluated.
const _ASCII_KEYWORDS = Dict{Symbol,Symbol}(
    :mu=>:μ, :beta=>:β, :gamma=>:γ, :xi=>:ξ, :kappa=>:κ,
    :alpha=>:α, :tau=>:τ, :eta=>:η, :sigma=>:σ, :phi=>:ϕ,
    :theta_o=>:θₒ, :theta_e=>:θₑ, :theta_o_system=>:θₒSystem,
    :phi_0=>:ϕ₀, :eta_o=>:ηₒ, :eta_1=>:η₁, :r_s=>:rₛ,
    :r_bar=>:r̄, :n_phi=>:nϕ, :alpha_rm=>:αRM,
    :phi_min=>:ϕMin, :phi_max=>:ϕMax, :tau_cutoff=>:τCutOff,
    :v_r_frac=>:vᵣFrac, :sigma_rho_r=>:σρᵣ, :sigma_rho_c=>:σρc,
    :sigma_rho_relative=>:σρRelative,
    :sigma_theta_r=>:σΘᵣ, :sigma_theta_c=>:σΘc, :sigma_t=>:σₜ,
    :delta_a=>:ΔA, :delta_r=>:Δr, :delta_phi=>:Δϕ)

"""
    acceptedKeywords(f)

Optional declaration of additional canonical callback keywords. `nothing` (the default)
means unregistered; `()` registers a callback with no additional options. For example,
`BLR.acceptedKeywords(::typeof(myIntensityFx)) = (:level,)` enables strict checking.
Unicode and ASCII keyword names may be mixed, but passing both spellings of one parameter
raises `ArgumentError`. Use canonical Unicode symbols in declarations; shared geometry context
is supplied separately. Aliases are listed with the parameters in the relevant constructor/helper docstrings.
Aliases apply only to keyword arguments: struct fields and recorded parameters retain Unicode
names, e.g. `m.camera.α` and `m.params.μ` are allowed but `m.camera.alpha` is not.

Unknown keywords are errors when all selected callbacks are registered. If any callback is
unregistered, unchecked extras are forwarded to it with one grouped warning; recognized options
for inactive built-in physics are accepted with an unused-option warning. No warning occurs
merely because a callback is unregistered. Callback exceptions retain their original details.
Declarations are authoritative; BLR cannot determine whether callback code uses an option.
"""
acceptedKeywords(f) = nothing
_keyword_impl(f) = f
const _CALLBACK_CONTEXT = (:r, :rₛ, :i, :ϕ, :ϕ₀, :rot, :θₒ, :θₒSystem, :reflect,
    :rMin, :rMax, :rng, :ΔA, :Δr, :Δϕ, :scale, :x, :y, :z)
const _RING_KEYS = (_CALLBACK_CONTEXT..., :I, :v, :τ, :η)
const _BIN_KEYS = (:bins, :overflow, :centered, :minX, :maxX)

function _normalize_keywords(kwargs)
    ks = Symbol[]
    vs = Any[]
    for (k, v) in pairs(kwargs)
        canonical = get(_ASCII_KEYWORDS, k, k)
        canonical in ks && throw(ArgumentError("duplicate keyword spellings for `$canonical`"))
        push!(ks, canonical); push!(vs, v)
    end
    return NamedTuple{Tuple(ks)}(Tuple(vs))
end

function _check_keywords(kw, own, callbacks=(); reserved=(), inactive=false, warn=true, diagnostics=nothing)
    for k in keys(kw)
        k in reserved && throw(ArgumentError("keyword `$k` would overwrite a positional or generated parameter"))
    end
    declarations = map(acceptedKeywords, callbacks)
    for d in declarations
        d === nothing || (d isa Tuple && all(k -> k isa Symbol, d)) ||
            throw(ArgumentError("acceptedKeywords must return nothing or a tuple of symbols"))
    end
    active = Set{Symbol}(own)
    for d in declarations
        d === nothing || union!(active, d)
    end
    extra = [k for k in keys(kw) if !(k in active)]
    unknown = [k for k in extra if !(inactive && k in _physics_keywords())]
    unchecked = any(isnothing, declarations)
    if !isempty(unknown) && !unchecked
        hint = any(k -> k in (:θ, :θ₀), unknown) ? "; use θₒ/theta_o (letter o) for opening angle" : ""
        throw(ArgumentError("unknown keyword(s): $(join(unknown, ", "))$hint"))
    end
    if (warn || diagnostics !== nothing) && !isempty(extra)
        if unchecked
            f = callbacks[findfirst(isnothing, declarations)]
            message = "Keywords $(join(extra, ", ")) were not checked. To enable validation, define BLR.acceptedKeywords(::typeof($f)) = ($(join((repr(k) for k in extra), ", ")),)."
        else
            message = "Unused keyword(s) for the selected physics: $(join(extra, ", "))"
        end
        warn && @warn message
        diagnostics === nothing || push!(diagnostics, message)
    end
    return kw
end

# Built-in declarations are installed by @validated next to their implementations.
const _BUILTIN_CALLBACKS = Any[]
_physics_keywords() = Iterators.flatten(acceptedKeywords(f) for f in _BUILTIN_CALLBACKS)
function _callback(f, kw)
    d = acceptedKeywords(f)
    selected = d === nothing ? kw : (; (k=>v for (k,v) in pairs(kw)
        if k in d || (!(f in _BUILTIN_CALLBACKS) && k in _CALLBACK_CONTEXT))...)
    if f in _BUILTIN_CALLBACKS
        return _keyword_impl(f)(; selected...)
    end
    # A custom adapter may delegate its shared context to a built-in helper. The
    # construction boundary has already checked those keys; avoid per-point checks
    # and duplicate diagnostics in that delegation. Task-local state is thread safe.
    tls = task_local_storage()
    previous = get(tls, :BLR_callback_context, false)
    tls[:BLR_callback_context] = true
    try
        return f(; selected...)
    finally
        tls[:BLR_callback_context] = previous
    end
end

function _construction_keywords(kw, own, defaults; reserved=(), warn=true, diagnostics=nothing)
    callbacks = Function[]
    for (k, default) in pairs(defaults)
        f = get(kw, k, default)
        f isa Function && push!(callbacks, f)
    end
    for k in (:r, :η, :τ)
        if !haskey(defaults, k) && haskey(kw, k) && kw[k] isa Function
            push!(callbacks, kw[k])
        end
    end
    # Response is part of construction even when no response override is supplied.
    push!(callbacks, response)
    _check_keywords(kw, own, Tuple(callbacks); reserved=reserved, inactive=true, warn=warn, diagnostics=diagnostics)
end

# Generate a thin public normalizer and retain the original typed keyword implementation.
# mode is :plain, :callback, :cloud, :disk, or :bins. Internal model construction calls
# the typed implementations directly, so checks and warnings stay outside point loops.
macro validated(mode, def)
    def.head == :function || error("@validated requires a function definition")
    sig, body = def.args
    name = sig.args[1]
    pars = sig.args[2]
    pars isa Expr && pars.head == :parameters || error("keyword signature required")
    positional = sig.args[3:end]
    keyname(x) = x isa Symbol ? x : x.head in (:kw, :(::), :(...)) ? keyname(x.args[1]) : error("unsupported argument $x")
    keys = [keyname(k) for k in pars.args if !(k isa Expr && k.head == :(...))]
    poskeys = [keyname(k) for k in positional]
    impl = Symbol("_kw_", name)
    inner = Expr(:function, Expr(:call, impl, pars, positional...), body)
    kw = gensym(:kwargs); normalized = gensym(:normalized)
    own = Expr(:tuple, QuoteNode.(keys)...)
    reservedkeys = [poskeys..., :constructor]
    mvalue = mode isa QuoteNode ? mode.value : mode
    if mvalue in (:cloud, :disk)
        append!(reservedkeys, [k for k in (:r, :ϕ, :ϕ₀, :θₒ, :θₒSystem, :ΔA, :Δr, :Δϕ, :x, :y, :z, :nClouds) if !(k in keys)])
    end
    if mvalue == :cloud
        :rot in keys || push!(reservedkeys, :rot)
        push!(reservedkeys, :reflect)
    end
    reserved = Expr(:tuple, QuoteNode.(unique(reservedkeys))...)
    m = mode isa QuoteNode ? mode.value : mode
    check = if m in (:cloud, :disk)
        defaults = m == :cloud ? :( (; I=IsotropicIntensity, v=vCircularCloud) ) : :( (; I=DiskWindIntensity, v=vCircularDisk) )
        :(_construction_keywords($normalized, ($own..., :τ, :η, :rot, :reflect, :rng), $defaults; reserved=$reserved))
    elseif m == :resident
        target = name in (:gpuCloudModel, :residentCloudModel) ? :residentCloudModel : :residentDiskWindModel
        :(_check_resident_keywords($normalized, ($own..., _resident_keywords(Val($(QuoteNode(target))))...), Val($(QuoteNode(target))); reserved=$reserved))
    elseif m == :bins
        :(_check_keywords($normalized, ($own..., $(name == :getProfile ? :(_profile_keywords(name)) : :(_BIN_KEYS))...); reserved=$reserved))
    elseif m == :callback
        :(get(task_local_storage(), :BLR_callback_context, false) ||
            _check_keywords($normalized, $own; reserved=$reserved, inactive=true))
    else
        :(_check_keywords($normalized, $own; reserved=$reserved))
    end
    callargs = [k isa Expr && k.head == :kw ? keyname(k) : keyname(k) for k in positional]
    wrapper = Expr(:function, Expr(:call, name, Expr(:parameters, Expr(:(...), kw)), positional...), quote
        $normalized = _normalize_keywords($kw)
        $check
        $impl($(callargs...); $normalized...)
    end)
    spec = m in (:cloud, :disk) ? :(_constructor_spec(::Val{$(QuoteNode(name))}, ::Val{$(length(positional))}) = ($own, $reserved, $(QuoteNode(m)))) : nothing
    declaration = m == :callback ? quote
        acceptedKeywords(::typeof($name)) = $own
        _keyword_impl(::typeof($name)) = $impl
        push!(_BUILTIN_CALLBACKS, $name)
    end : nothing
    esc(quote
        $inner
        Base.@__doc__ $wrapper
        $declaration
        $spec
    end)
end

_resident_keywords(::Val{:residentDiskWindModel}) = (:nr, :nϕ, :scale, :rot, :θₒ,
    :f1, :f2, :f3, :f4, :α, :ηₒ, :η₁, :αRM, :rNorm, :rₛ, :vᵣFrac, :inflow,
    :τ, :intensity, :velocity, :backend, :T)
_resident_keywords(::Val{:residentCloudModel}) = (:μ, :β, :F, :rₛ, :θₒ, :γ, :ξ, :i,
    :rescale, :ηₒ, :η₁, :αRM, :rNorm, :intensity, :κ, :velocity, :σρᵣ, :σρc,
    :σΘᵣ, :σΘc, :θₑ, :fEllipse, :fFlow, :σₜ, :σρRelative, :τ, :backend, :T)

function _profile_keywords(name)
    n = name isa Function ? :function : Symbol(name)
    n === :ratio && return (_BIN_KEYS..., :floor)
    n in (:phase, :moment2) ? (_BIN_KEYS..., :U, :V, :PA, :BLRAng, :returnAvg, :offAxisInds) : _BIN_KEYS
end

function _check_resident_keywords(kw, own, target; reserved=())
    _check_keywords(kw, own; reserved=reserved)
    inactive = Symbol[]
    if target === Val(:residentCloudModel)
        if get(kw, :intensity, :isotropic) === :isotropic
            push!(inactive, :κ)
        elseif get(kw, :intensity, :isotropic) === :cloud
            push!(inactive, :rescale)
        end
        get(kw, :velocity, :circular) === :circular &&
            append!(inactive, (:σρᵣ, :σρc, :σΘᵣ, :σΘc, :θₑ, :fEllipse, :fFlow, :σₜ, :σρRelative))
    else
        get(kw, :intensity, nothing) === nothing || append!(inactive, (:f1, :f2, :f3, :f4, :α))
        get(kw, :velocity, nothing) === nothing || append!(inactive, (:rₛ, :vᵣFrac, :inflow))
    end
    unused = [k for k in keys(kw) if k in inactive]
    isempty(unused) || @warn "Unused keyword(s) for the selected physics: $(join(unused, ", "))"
    kw
end
