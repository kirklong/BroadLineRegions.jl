# Keyword aliases and callback validation

BLR accepts ASCII keyword names alongside its existing Unicode names. Function names,
fields, keyword values, units, positional arguments, and physical parameterizations are unchanged.
For example:

```julia
using BroadLineRegions
m = BLR.cloudModel(100_000;
    mu=500.0, beta=1.0, F=0.5, theta_o=π/6,
    gamma=1.0, xi=1.0, i=0.0,
    I=BLR.IsotropicIntensity, v=BLR.vCircularCloud,
    rescale=1e-5, tau=0.0)
```

Aliases apply wherever the corresponding canonical keyword is supported, including direct
physics helpers, host and resident constructors, GPU constructors, `DRW`, raytracing,
`rebuild` overrides, and parameter reuse in `addLine!`. Positional inputs remain positional;
`r_bar` is useful for reconstruction overrides of the existing `r̄` parameter.
Provenance in `model.params` uses canonical Unicode names.

| ASCII | Canonical | ASCII | Canonical |
|:--|:--|:--|:--|
| `mu` | `μ` | `beta` | `β` |
| `gamma` | `γ` | `xi` | `ξ` |
| `kappa` | `κ` | `alpha` | `α` |
| `tau` | `τ` | `eta` | `η` |
| `sigma` | `σ` | `phi` | `ϕ` |
| `theta_o` | `θₒ` | `theta_e` | `θₑ` |
| `theta_o_system` | `θₒSystem` | `phi_0` | `ϕ₀` |
| `eta_o` | `ηₒ` | `eta_1` | `η₁` |
| `r_s` | `rₛ` | `r_bar` | `r̄` |
| `n_phi` | `nϕ` | `alpha_rm` | `αRM` |
| `phi_min` | `ϕMin` | `phi_max` | `ϕMax` |
| `tau_cutoff` | `τCutOff` | `v_r_frac` | `vᵣFrac` |
| `sigma_rho_r` | `σρᵣ` | `sigma_rho_c` | `σρc` |
| `sigma_theta_r` | `σΘᵣ` | `sigma_theta_c` | `σΘc` |
| `sigma_t` | `σₜ` | `delta_a` | `ΔA` |
| `delta_r` | `Δr` | `delta_phi` | `Δϕ` |

The `o` in `theta_o` and `eta_o` is a letter. The `0` in `phi_0` is a digit.
Neither `θ` nor `θ₀` is an alias for opening angle: use `θₒ` or `theta_o`.
Supplying both spellings of one parameter raises `ArgumentError`, even for equal values.

Built-in calls reject unknown keywords before drawing random numbers, allocating model
storage, or mutating an input. Recognized physics options for an inactive built-in callback
remain accepted, with one grouped warning listing explicitly supplied unused options.
Internal context and defaults do not generate those warnings. Numerical binning options retain
their existing meanings; ordinary plotting attributes are still validated by Plots.

## Custom callbacks

Custom callbacks can continue accepting extra keywords:

```julia
myIntensity(; r, gain=1.0, kwargs...) = gain / r
m = BLR.cloudModel(10; I=myIntensity, gain=2.0)
```

Here `gain` is forwarded, with a grouped warning that it was not checked. Register additional
canonical keywords to enable strict validation:

```julia
BLR.acceptedKeywords(::typeof(myIntensity)) = (:gain,)
m = BLR.cloudModel(10; I=myIntensity, gain=2.0)  # no unchecked-key warning
```

`acceptedKeywords(f) === nothing` means unregistered, `()` registers no additional
keywords, and a tuple of symbols declares the additional accepted options. Shared geometry
context remains available; callbacks should accept `kwargs...` for that context. Declarations
are authoritative: BLR cannot prove that a callback uses a declared option.

If any selected callback is unregistered, undeclared extras are forwarded to those
unregistered callbacks with one grouped warning. An inactive built-in option is described
as unchecked when such a callback could consume it. Registered callbacks receive their
declared options and shared context; built-ins receive their own relevant arguments.
No warning occurs merely because an unregistered callback is selected. Once every selected
callback is registered, unknown extras raise `ArgumentError`. Original callback exceptions
retain their details.

Registration is also respected when reconstructing recorded parameters with `rebuild`
or `addLine!`. Overrides retain the existing rule: a key must already occur in at least one
leaf record, and it applies only to leaves containing that key.
