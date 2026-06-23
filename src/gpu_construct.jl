using KernelAbstractions

# On-device model construction (B2 pipeline). Builds the flat `ModelArrays` columns for a DiskWind
# model directly on a `KernelAbstractions` backend via ONE fused kernel — deprojection + velocity +
# intensity + reverberation response, one launch / one pass per pixel — and generates the camera
# (α,β,ΔA) grid in-kernel from scalars, so nothing touches a host `ring`. The result is wrapped as a
# `ResidentModel` ready for the observables pipeline. `gpu(m::model)` (host build → flatten →
# transfer) remains the path for custom `I`/`v` Functions; this path is for the built-in DiskWind
# physics (or a user-supplied GPU-safe scalar `(r, ϕ, inc) -> value`, see below).
#
# Point ordering matches `flatten(DiskWindModel(...))` exactly: flattened index `p = k + (j-1)*nr`
# with radial ring `k` the fast axis and azimuth `j` the slow axis (see `expandPerPoint`). Boundary
# masking follows the plan's B2 rule — bounds are host-rounded to 9 significant digits and compared
# with a plain `<`/`>` against the (unrounded) per-point radius, so a single boundary ring can differ
# from the host's per-point `round` by inclusion (documented, measure-zero for generic bounds).

"""
    _rt_disk_intensity_fn(::Type{T}, f1, f2, f3, f4, α) -> closure

GPU-safe scalar intensity closure for the built-in `DiskWindIntensity`, capturing the (T-typed)
shape factors so the fused construction kernel can call `Ifun(r, ϕ, inc)`. Mirrors
`DiskWind_I_(r, ϕ, i, f1, f2, f3, f4, α)`.
"""
_rt_disk_intensity_fn(::Type{T}, f1, f2, f3, f4, α) where {T<:Real} =
    let f1 = T(f1), f2 = T(f2), f3 = T(f3), f4 = T(f4), α = T(α)
        (r, ϕ, inc) -> _rt_disk_wind_i_scalar(r, ϕ, inc, f1, f2, f3, f4, α)
    end

"""
    _rt_disk_velocity_fn(::Type{T}, rₛ) -> closure

GPU-safe scalar velocity closure for the built-in `vCircularDisk`, capturing the (T-typed) Schwarzschild
radius so the fused construction kernel can call `vfun(r, ϕ, inc)`.
"""
_rt_disk_velocity_fn(::Type{T}, rₛ) where {T<:Real} =
    let rₛ = T(rₛ)
        (r, ϕ, inc) -> _rt_v_circular_disk_scalar(r, ϕ, inc, rₛ)
    end

# vCircularRadialDisk: a Keplerian + radial in/outflow split (`vᵣFrac` of v goes radial, `inflow`
# flips its sign). Reduces exactly to `vCircularDisk` at `vᵣFrac = 0`.
@inline function _rt_v_circular_radial_disk_scalar(r::T, ϕ::T, inc::T, vᵣFrac::T, inflowSign::T, rₛ::T) where {T<:Real}
    vsini = sqrt(rₛ / (2 * r)) * sin(inc)
    vᵣ = vsini * cos(ϕ) * vᵣFrac * inflowSign
    vϕ = vsini * sin(ϕ) * (one(T) - vᵣFrac)
    return -(vᵣ + vϕ)
end

"""
    _rt_disk_radial_velocity_fn(::Type{T}, vᵣFrac, inflow, rₛ) -> closure

GPU-safe scalar velocity closure for the built-in `vCircularRadialDisk` (Keplerian + radial
in/outflow). With `vᵣFrac = 0` it is identical to [`vCircularDisk`](@ref).
"""
_rt_disk_radial_velocity_fn(::Type{T}, vᵣFrac, inflow::Bool, rₛ) where {T<:Real} =
    let vf = T(vᵣFrac), s = inflow ? -one(T) : one(T), rₛ = T(rₛ)
        (r, ϕ, inc) -> _rt_v_circular_radial_disk_scalar(r, ϕ, inc, vf, s, rₛ)
    end

# One fused kernel: camera grid (from scalars) -> deproject -> velocity/intensity -> response.
# `Ifun`/`vfun` are GPU-safe scalar callables; out-of-range pixels deproject to NaN geometry and get
# NaN v/I (matching the host, where `DiskWind_I_`/`vCircularDisk` evaluate NaN radii to NaN). Per-ring
# scalars (i, rot, θₒ, τ=0, reflect=false) are written for every pixel, including out-of-range ones,
# exactly as `flatten` expands them.
@kernel function _rt_build_disk_kernel!(outr, outϕ, outϕ₀, outi, outrot, outθₒ, outv, outI, outΔA,
        outτ, outη, outx, outy, outz, outα, outβ, outreflect, Ifun, vfun,
        nr, scaleLog, rStart, Δr, Δϕ, inc, rot, θₒ, m11, m12, m21, m22,
        r3d11, r3d12, r3d21, r3d22, r3d31, r3d32, rMinR, rMaxR, ηₒ, η₁, αRM, rNorm, ΔAfac)
    p = @index(Global)
    T = eltype(outr)
    k = (p - 1) % nr + 1
    j = (p - 1) ÷ nr + 1
    rc = scaleLog ? exp(rStart + (k - 1) * Δr) : (rStart + (k - 1) * Δr)
    ϕc = -T(π) + (j - 1) * Δϕ
    sinϕc, cosϕc = sincos(ϕc)
    a = rc * cosϕc
    b = rc * sinϕc
    ΔA = (scaleLog ? rc * rc : rc) * ΔAfac
    r, ϕ, ϕ₀, η, x, y, z = _rt_disk_deproject_scalar(a, b, inc, rot, θₒ, m11, m12, m21, m22,
        r3d11, r3d12, r3d21, r3d22, r3d31, r3d32, rMinR, rMaxR, ηₒ, η₁, αRM, rNorm)
    nan = convert(T, NaN)
    if isnan(r)
        outv[p] = nan
        outI[p] = nan
    else
        outv[p] = vfun(r, ϕ, inc)
        outI[p] = Ifun(r, ϕ, inc)
    end
    outr[p] = r
    outϕ[p] = ϕ
    outϕ₀[p] = ϕ₀
    outη[p] = η
    outx[p] = x
    outy[p] = y
    outz[p] = z
    outα[p] = a
    outβ[p] = b
    outΔA[p] = ΔA
    outi[p] = inc
    outrot[p] = rot
    outθₒ[p] = θₒ
    outτ[p] = zero(T)
    outreflect[p] = false
end

"""
    _build_diskwind_modelarrays(rMin, rMax, inc, nr, nϕ, scale, rot, θₒ, Ifun, vfun;
        ηₒ, η₁, αRM, rNorm, backend, T) -> ModelArrays

Allocate the device columns on `backend` and run `_rt_build_disk_kernel!` to fill them, then
wrap as a [`ModelArrays`](@ref). `Ifun`/`vfun` are GPU-safe scalar callables `(r, ϕ, inc) -> value`.
"""
function _build_diskwind_modelarrays(rMin::Real, rMax::Real, inc::Real, nr::Int, nϕ::Int,
        scale::Symbol, rot::Real, θₒ::Real, Ifun, vfun; ηₒ::Real, η₁::Real, αRM::Real, rNorm::Real,
        backend=KernelAbstractions.CPU(), T=Float64)
    rMin < rMax || throw(ArgumentError("rMin must be less than rMax"))
    nr > 1 || throw(ArgumentError("nr must be greater than 1"))
    nϕ > 1 || throw(ArgumentError("nϕ must be greater than 1"))
    scale === :log || scale === :linear || throw(ArgumentError("scale must be :log or :linear"))

    Δϕ = T(2π / nϕ)                       # = ϕ[2]-ϕ[1] for range(-π,π,nϕ+1)[1:end-1]
    scaleLog = scale === :log
    rStart = scaleLog ? T(log(rMin * cos(inc))) : T(rMin * cos(inc))
    rEnd = scaleLog ? T(log(rMax)) : T(rMax)
    Δr = (rEnd - rStart) / (nr - 1)
    ΔAfac = Δr * Δϕ

    r3D = get_r3D(inc, rot, θₒ)
    undoTilt = [sin(inc) 0.0 -cos(inc); 0.0 1.0 0.0; cos(inc) 0.0 sin(inc)]
    M = undoTilt * r3D
    rMinR = T(round(rMin, sigdigits=9))
    rMaxR = T(round(rMax, sigdigits=9))

    n = nr * nϕ
    alloc() = KernelAbstractions.allocate(backend, T, n)
    outr, outϕ, outϕ₀, outi, outrot, outθₒ = alloc(), alloc(), alloc(), alloc(), alloc(), alloc()
    outv, outI, outΔA, outτ, outη = alloc(), alloc(), alloc(), alloc(), alloc()
    outx, outy, outz, outα, outβ = alloc(), alloc(), alloc(), alloc(), alloc()
    outreflect = KernelAbstractions.allocate(backend, Bool, n)

    kernel! = _rt_build_disk_kernel!(backend)
    event = kernel!(outr, outϕ, outϕ₀, outi, outrot, outθₒ, outv, outI, outΔA, outτ, outη,
        outx, outy, outz, outα, outβ, outreflect, Ifun, vfun, nr, scaleLog, rStart, Δr, Δϕ,
        T(inc), T(rot), T(θₒ), T(M[1, 1]), T(M[1, 2]), T(M[2, 1]), T(M[2, 2]),
        T(r3D[1, 1]), T(r3D[1, 2]), T(r3D[2, 1]), T(r3D[2, 2]), T(r3D[3, 1]), T(r3D[3, 2]),
        rMinR, rMaxR, T(ηₒ), T(η₁), T(αRM), T(rNorm), ΔAfac; ndrange=n)
    event !== nothing && wait(event)

    ma = ModelArrays{T,typeof(outr),typeof(outreflect)}(outr, outϕ, outϕ₀, outi, outrot, outθₒ,
        outv, outI, outΔA, outτ, outη, outx, outy, outz, outα, outβ, outreflect)
    return ma
end

"""
    residentDiskWindModel(rMin, rMax, i; nr=128, nϕ=256, scale=:log, rot=0.0, θₒ=0.0,
        f1, f2, f3, f4, α, ηₒ=0.5, η₁=0.5, αRM=0.0, rNorm=1.0, rₛ=1.0,
        intensity=nothing, velocity=nothing, backend=KernelAbstractions.CPU(), T=Float64) -> ResidentModel

Build a DiskWind model **directly on `backend`** as a device-resident [`ResidentModel`](@ref), without
constructing host `ring`s — the on-device counterpart of `DiskWindModel(...) |> resident`. One fused
kernel generates the camera grid and computes geometry, velocity, intensity, and reverberation
response per pixel. Reproduces `flatten(DiskWindModel(rMin, rMax, i; ...))` to floating-point rounding
(see the boundary-masking note in `src/gpu_construct.jl`).

The default physics is the built-in `DiskWindIntensity` (`f1`–`f4`, `α`) and `vCircularDisk` (`rₛ`). Pass
`vᵣFrac`/`inflow` to switch the velocity to the built-in `vCircularRadialDisk` (Keplerian + radial
in/outflow; `vᵣFrac = 0` is `vCircularDisk`). Advanced users may instead pass `intensity`/`velocity` as
**GPU-safe scalar callables** `(r, ϕ, inc) -> value` (e.g. a `let`-closure over isbits parameters) to run
custom physics in the same fused kernel; generic custom `Function`s that cannot run in a kernel should
use host construction + `gpu(m)` instead.

Use `gpu`-backed construction (CUDA loaded) via `gpuDiskWindModel`; the default `CPU()` backend builds
on the host and is mainly for testing the on-device pipeline without a GPU.
"""
function residentDiskWindModel(rMin::Real, rMax::Real, i::Real; nr::Int=128, nϕ::Int=256,
        scale::Symbol=:log, rot::Real=0.0, θₒ::Real=0.0,
        f1::Real=NaN, f2::Real=NaN, f3::Real=NaN, f4::Real=NaN, α::Real=NaN,
        ηₒ::Real=0.5, η₁::Real=0.5, αRM::Real=0.0, rNorm::Real=1.0, rₛ::Real=1.0,
        vᵣFrac::Real=0.0, inflow::Bool=true,
        intensity=nothing, velocity=nothing,
        backend=KernelAbstractions.CPU(), T=Float64)
    T <: Real || throw(ArgumentError("element type T must be <: Real, got $T"))
    Ifun = if intensity === nothing
        any(isnan, (f1, f2, f3, f4, α)) &&
            throw(ArgumentError("built-in DiskWind intensity requires f1, f2, f3, f4, α (or pass intensity=...)"))
        _rt_disk_intensity_fn(T, f1, f2, f3, f4, α)
    else
        intensity
    end
    vfun = velocity === nothing ? _rt_disk_radial_velocity_fn(T, vᵣFrac, inflow, rₛ) : velocity
    ma = _build_diskwind_modelarrays(rMin, rMax, i, nr, nϕ, scale, rot, θₒ, Ifun, vfun;
        ηₒ=ηₒ, η₁=η₁, αRM=αRM, rNorm=rNorm, backend=backend, T=T)
    meta = _rt_diskwind_meta(ma, rMin, rMax, i, nr, nϕ, scale, backend, T)
    return ResidentModel(ma, backend, 1, meta)
end

"""
    gpuDiskWindModel(rMin, rMax, i; T=Float32, kwargs...) -> ResidentModel
    gpuDiskWindModel(r̄, rFac, Sα, i; T=Float32, kwargs...) -> ResidentModel

Build a DiskWind model directly on the GPU (device-resident [`ResidentModel`](@ref)) — the CUDA
counterpart of [`residentDiskWindModel`](@ref), defaulting to `Float32` (GeForce FP64 is ~1/64 rate).
Requires CUDA.jl loaded (defined by the CUDA package extension); errors otherwise.
"""
function gpuDiskWindModel end

"""
    residentDiskWindModel(r̄, rFac, Sα, i; rot=0.0, nr=128, nϕ=256, scale=:log, kwargs...) -> ResidentModel

`r̄`/`rFac`/`Sα` parameterization (matching `DiskWindModel(r̄, rFac, α, i; ...)`): derives `rMin`/`rMax`
via [`get_rMinMaxDiskWind`](@ref) and sets the intensity power-law `α = Sα`.
"""
function residentDiskWindModel(r̄::Real, rFac::Real, Sα::Real, i::Real; rot::Real=0.0,
        nr::Int=128, nϕ::Int=256, scale::Symbol=:log, kwargs...)
    (Sα != 1 / 2 && Sα != 3 / 2) || throw(ArgumentError("Sα cannot be 1/2 or 3/2 (divides by zero)"))
    r̄ > 0 || throw(ArgumentError("r̄ must be greater than 0"))
    rMin, rMax = get_rMinMaxDiskWind(Float64(r̄), Float64(rFac), Float64(Sα))
    return residentDiskWindModel(rMin, rMax, i; nr=nr, nϕ=nϕ, scale=scale, rot=rot, α=Sα, kwargs...)
end

# ======================================================================================
# On-device cloud construction (B3 on the GPU). A `cloudModel` built directly on the device.
#
# Unlike the disk path this does NOT bit-match the host: the host `:philox` cloud path draws the
# radius from a `Gamma` via Distributions.jl's rejection sampler, which cannot be reproduced bit-for-
# bit inside a kernel. Instead each cloud gets an independent counter-based (Philox4x32) substream and
# samples the SAME target distributions on-device, so the cloud population is deterministic and
# seed-reproducible (and identical across CPU/CUDA backends) but statistically — not bitwise — equal to
# the host `:philox` realization. This matches the B3 design note ("you cannot bit-match MT and
# parallelize simultaneously"; the counter-based path is a different, fully deterministic stream).
# ======================================================================================

# --- Philox4x32-10 counter-based RNG (GPU-safe: only UInt32*UInt32 -> UInt64 widening) ---
const _PHILOX_M0 = 0xD2511F53 % UInt32
const _PHILOX_M1 = 0xCD9E8D57 % UInt32
const _PHILOX_W0 = 0x9E3779B9 % UInt32   # golden ratio
const _PHILOX_W1 = 0xBB67AE85 % UInt32   # sqrt(3)-1

@inline function _philox_mulhilo(a::UInt32, b::UInt32)
    prod = UInt64(a) * UInt64(b)
    return (UInt32(prod >> 32), UInt32(prod & 0xffffffff))   # (hi, lo)
end

@inline function _philox4x32_bijection(c0::UInt32, c1::UInt32, c2::UInt32, c3::UInt32,
        key0::UInt32, key1::UInt32)
    k0 = key0
    k1 = key1
    @inbounds for r in 1:10
        hi0, lo0 = _philox_mulhilo(_PHILOX_M0, c0)
        hi1, lo1 = _philox_mulhilo(_PHILOX_M1, c2)
        c0 = hi1 ⊻ c1 ⊻ k0
        c1 = lo1
        c2 = hi0 ⊻ c3 ⊻ k1
        c3 = lo0
        if r < 10
            k0 += _PHILOX_W0
            k1 += _PHILOX_W1
        end
    end
    return (c0, c1, c2, c3)
end

# The n-th uniform in (0,1) of cloud `cloudInd`'s substream: counter = (n÷4, cloudInd, 0, 0), word n%4.
@inline function _cloud_uniform(::Type{T}, key0::UInt32, key1::UInt32, cloudInd::UInt32, n::Integer) where {T<:Real}
    block = UInt32((n >> 2) & 0xffffffff)
    word = (n & 3) + 1
    out = _philox4x32_bijection(block, cloudInd, UInt32(0), UInt32(0), key0, key1)
    u = @inbounds out[word]
    return (T(u) + T(0.5)) * T(2.3283064365386963e-10)   # (u + 0.5) / 2^32  ∈ (0,1)
end

# --- Gamma(shape, 1) via Marsaglia–Tsang, consuming the substream; threads the draw counter `n`. ---
@inline function _cloud_gamma(::Type{T}, key0::UInt32, key1::UInt32, cloudInd::UInt32,
        shape::T, n0::Int) where {T<:Real}
    n = n0
    boost = one(T)
    a = shape
    if a < one(T)                          # shape < 1: draw Gamma(a+1) and rescale by U^(1/a)
        u = _cloud_uniform(T, key0, key1, cloudInd, n); n += 1
        boost = u^(one(T) / a)
        a += one(T)
    end
    d = a - T(1) / 3
    c = one(T) / sqrt(9 * d)
    g = d                                  # fallback (acceptance is ~0.95/iter; loop essentially always returns)
    @inbounds for _ in 1:128
        u1 = _cloud_uniform(T, key0, key1, cloudInd, n); n += 1
        u2 = _cloud_uniform(T, key0, key1, cloudInd, n); n += 1
        z = sqrt(-2 * log(u1)) * cos(T(2π) * u2)          # standard normal (Box–Muller)
        v = (one(T) + c * z)^3
        if v > zero(T)
            u3 = _cloud_uniform(T, key0, key1, cloudInd, n); n += 1
            if log(u3) < T(0.5) * z * z + d - d * v + d * log(v)
                g = d * v
                break
            end
        else
            n += 1                          # keep the per-iteration stream stride uniform (u3 slot)
        end
    end
    return boost * g, n
end

# --- GPU-safe cloud geometry (generic ports of rotate3D_scalar / reflect_scalar / vCircularCloud) ---
@inline function _rt_reflect3(x::T, y::T, z::T, sini::T, cosi::T) where {T<:Real}
    cti = cosi / sini
    den = one(T) + cti * cti
    xf = (x * (one(T) - cti * cti) - 2 * z * cti) / den
    zf = (z * (cti * cti - one(T)) - 2 * x * cti) / den
    return (xf, y, zf)
end

@inline function _rt_rotate3_vector(vx::T, vy::T, vz::T, sini::T, cosi::T, sinrot::T, cosrot::T,
        sinθₒ::T, cosθₒ::T) where {T<:Real}
    x = (cosθₒ * cosrot * sini - sinθₒ * cosi) * vx - sinrot * sini * vy + (sinθₒ * cosrot * sini + cosθₒ * cosi) * vz
    y = cosθₒ * sinrot * vx + cosrot * vy + sinθₒ * sinrot * vz
    z = -(cosθₒ * cosrot * cosi + sinθₒ * sini) * vx + sinrot * cosi * vy + (cosθₒ * sini - sinθₒ * cosrot * cosi) * vz
    return (x, y, z)
end

# Standard normal via Box–Muller; advances the draw counter by 2 (returns `(z, n+2)`).
@inline function _cloud_normal(::Type{T}, key0::UInt32, key1::UInt32, cloudInd::UInt32, n0::Int) where {T<:Real}
    u1 = _cloud_uniform(T, key0, key1, cloudInd, n0)
    u2 = _cloud_uniform(T, key0, key1, cloudInd, n0 + 1)
    return sqrt(-2 * log(u1)) * cos(T(2π) * u2), n0 + 2
end

# Per-cloud build: draw (ϕ₀, θ, rot, r, reflect) from the substream, then the deterministic geometry,
# velocity and intensity. Intensity is either isotropic (`rescale`) or `cloudIntensity` W(ϕw,κ); the
# velocity is either `vCircularCloud` (no draws) or `vCloudTurbulentEllipticalFlow` (consumes extra
# normals/uniform for turbulence + elliptical/flow). `useCloudI`/`useTurbulent` are scalar flags, so
# every thread takes the same branch (no warp divergence). Returns every ModelArrays column.
@inline function _rt_build_cloud_scalar(::Type{T}, key0::UInt32, key1::UInt32, cloudInd::UInt32,
        μ::T, β::T, F::T, rₛ::T, θₒ::T, γ::T, ξ::T, inc::T, rescale::T,
        ηₒ::T, η₁::T, αRM::T, rNorm::T, useCloudI::Bool, κ::T, useTurbulent::Bool,
        σρᵣ::T, σρc::T, σΘᵣ::T, σΘc::T, θₑ::T, fEllipse::T, fFlow::T, σₜ::T) where {T<:Real}
    sini, cosi = sincos(inc)
    ϕ₀ = _cloud_uniform(T, key0, key1, cloudInd, 0) * T(2π)
    uθ = _cloud_uniform(T, key0, key1, cloudInd, 1)
    cosθₒ_o = cos(θₒ)
    θ = acos(cosθₒ_o + (one(T) - cosθₒ_o) * uθ^γ)
    rot = _cloud_uniform(T, key0, key1, cloudInd, 2) * T(2π)
    g, n = _cloud_gamma(T, key0, key1, cloudInd, one(T) / (β * β), 3)
    r = rₛ + μ * F + μ * β * β * (one(T) - F) * g

    sinϕ₀, cosϕ₀ = sincos(ϕ₀)
    sinrot, cosrot = sincos(rot)
    sinθ, cosθ = sincos(θ)
    x, y, z = _rt_rotate3_vector(r * cosϕ₀, r * sinϕ₀, zero(T), sini, cosi, sinrot, cosrot, sinθ, cosθ)
    xPre = x                                   # pre-reflection x (cloudIntensity's ϕw uses the un-tilted frame)

    uref = _cloud_uniform(T, key0, key1, cloudInd, n); n += 1
    reflect = (z < -x * (cosi / sini)) && (uref > ξ)
    if reflect
        x, y, z = _rt_reflect3(x, y, z, sini, cosi)
    end

    xyPlaneX = sini * x - cosi * z
    ϕ = atan(y, xyPlaneX)
    η = ηₒ + η₁ * (r / rNorm)^αRM

    # intensity
    I = useCloudI ? (T(0.5) + κ * cos(atan(y, xPre))) : rescale   # W(ϕw,κ); ϕw in the pre-tilt frame

    # velocity
    vc = sqrt(rₛ / (2 * r))
    vLOS = if useTurbulent
        zt, n = _cloud_normal(T, key0, key1, cloudInd, n)
        vₜ = zt * σₜ * vc
        uE = _cloud_uniform(T, key0, key1, cloudInd, n); n += 1
        zρ, n = _cloud_normal(T, key0, key1, cloudInd, n)
        zΘ, n = _cloud_normal(T, key0, key1, cloudInd, n)
        ρ = uE < fEllipse ? vc + zρ * σρc : vc + zρ * σρᵣ
        Θ = uE < fEllipse ? T(π) / 2 + zΘ * σΘc :
            (fFlow < T(0.5) ? zΘ * σΘᵣ + (T(π) - θₑ) : zΘ * σΘᵣ + θₑ)
        vx = sqrt(T(2)) * ρ * cos(Θ)
        vy = ρ * sin(Θ)
        vX, _, _ = let
            a, b2, c = _rt_rotate3_vector(-(vx * cosϕ₀ + vy * sinϕ₀), -(vx * sinϕ₀ + vy * cosϕ₀),
                zero(T), sini, cosi, sinrot, cosrot, sinθ, cosθ)
            reflect ? _rt_reflect3(a, b2, c, sini, cosi) : (a, b2, c)
        end
        vX + vₜ
    else
        vX, _, _ = let
            a, b2, c = _rt_rotate3_vector(-vc * sinϕ₀, -vc * cosϕ₀, zero(T), sini, cosi, sinrot, cosrot, sinθ, cosθ)
            reflect ? _rt_reflect3(a, b2, c, sini, cosi) : (a, b2, c)
        end
        vX
    end

    return (r, ϕ, ϕ₀, inc, rot, θ, vLOS, I, one(T), zero(T), η, x, y, z, reflect)
end

@kernel function _rt_build_cloud_kernel!(outr, outϕ, outϕ₀, outi, outrot, outθₒ, outv, outI, outΔA,
        outτ, outη, outx, outy, outz, outα, outβ, outreflect, key0, key1,
        μ, β, F, rₛ, θₒ, γ, ξ, inc, rescale, ηₒ, η₁, αRM, rNorm, useCloudI, κ, useTurbulent,
        σρᵣ, σρc, σΘᵣ, σΘc, θₑ, fEllipse, fFlow, σₜ)
    p = @index(Global)
    T = eltype(outr)
    cloudInd = UInt32(p)
    r, ϕ, ϕ₀, i, rot, θ, v, I, ΔA, τ, η, x, y, z, reflect = _rt_build_cloud_scalar(
        T, key0, key1, cloudInd, μ, β, F, rₛ, θₒ, γ, ξ, inc, rescale, ηₒ, η₁, αRM, rNorm,
        useCloudI, κ, useTurbulent, σρᵣ, σρc, σΘᵣ, σΘc, θₑ, fEllipse, fFlow, σₜ)
    outr[p] = r
    outϕ[p] = ϕ
    outϕ₀[p] = ϕ₀
    outi[p] = i
    outrot[p] = rot
    outθₒ[p] = θ
    outv[p] = v
    outI[p] = I
    outΔA[p] = ΔA
    outτ[p] = τ
    outη[p] = η
    outx[p] = x
    outy[p] = y
    outz[p] = z
    outα[p] = y          # camera at +x: α = y, β = z (post-reflection), matching model(::Vector{ring})
    outβ[p] = z
    outreflect[p] = reflect
end

"""
    _build_cloud_modelarrays(nClouds, seed; μ, β, F, rₛ, θₒ, γ, ξ, i, rescale,
        ηₒ, η₁, αRM, rNorm, backend, T) -> ModelArrays

Allocate the device columns on `backend` and run `_rt_build_cloud_kernel!` to draw `nClouds`
clouds (one thread each) into a [`ModelArrays`](@ref).
"""
function _build_cloud_modelarrays(nClouds::Int, seed::Integer; μ::Real, β::Real, F::Real, rₛ::Real,
        θₒ::Real, γ::Real, ξ::Real, i::Real, rescale::Real, ηₒ::Real, η₁::Real, αRM::Real, rNorm::Real,
        useCloudI::Bool=false, κ::Real=0.0, useTurbulent::Bool=false, σρᵣ::Real=0.0, σρc::Real=0.0,
        σΘᵣ::Real=0.0, σΘc::Real=0.0, θₑ::Real=0.0, fEllipse::Real=0.0, fFlow::Real=0.0, σₜ::Real=0.0,
        backend=KernelAbstractions.CPU(), T=Float64)
    nClouds > 0 || throw(ArgumentError("nClouds must be positive"))
    s = seed % UInt64
    key0 = UInt32(s & 0xffffffff)
    key1 = UInt32((s >> 32) & 0xffffffff)

    n = nClouds
    alloc() = KernelAbstractions.allocate(backend, T, n)
    outr, outϕ, outϕ₀, outi, outrot, outθₒ = alloc(), alloc(), alloc(), alloc(), alloc(), alloc()
    outv, outI, outΔA, outτ, outη = alloc(), alloc(), alloc(), alloc(), alloc()
    outx, outy, outz, outα, outβ = alloc(), alloc(), alloc(), alloc(), alloc()
    outreflect = KernelAbstractions.allocate(backend, Bool, n)

    kernel! = _rt_build_cloud_kernel!(backend)
    event = kernel!(outr, outϕ, outϕ₀, outi, outrot, outθₒ, outv, outI, outΔA, outτ, outη,
        outx, outy, outz, outα, outβ, outreflect, key0, key1, T(μ), T(β), T(F), T(rₛ),
        T(θₒ), T(γ), T(ξ), T(i), T(rescale), T(ηₒ), T(η₁), T(αRM), T(rNorm),
        useCloudI, T(κ), useTurbulent, T(σρᵣ), T(σρc), T(σΘᵣ), T(σΘc), T(θₑ),
        T(fEllipse), T(fFlow), T(σₜ); ndrange=n)
    event !== nothing && wait(event)

    return ModelArrays{T,typeof(outr),typeof(outreflect)}(outr, outϕ, outϕ₀, outi, outrot, outθₒ,
        outv, outI, outΔA, outτ, outη, outx, outy, outz, outα, outβ, outreflect)
end

"""
    residentCloudModel(nClouds, seed; μ=500.0, β=1.0, F=0.5, rₛ=1.0, θₒ=π/2, γ=1.0, ξ=1.0, i=0.0,
        rescale=1.0, ηₒ=0.5, η₁=0.5, αRM=0.0, rNorm=1.0,
        backend=KernelAbstractions.CPU(), T=Float64) -> ResidentModel

Build a Pancoast-style cloud model **directly on `backend`** as a device-resident
[`ResidentModel`](@ref) — the on-device counterpart of `cloudModel(nClouds; rng=:philox, seed=…)`.
Each cloud is drawn from an independent counter-based (Philox4x32) substream keyed by `(seed, cloud
index)`, so the realization is deterministic, seed-reproducible, and identical across CPU/CUDA backends.

!!! note "Not bit-identical to the host"
    The host `Gamma` radius draw uses Distributions.jl's rejection sampler, which cannot be reproduced
    bitwise in a kernel. This path samples the same target distributions on-device, so it is
    *statistically* — not bitwise — equal to `cloudModel(...; rng=:philox, seed)`. Use the host path +
    `gpu(m)` if you need the exact legacy realization.

`seed` is required (the stream is counter-based).

Two built-in physics options are supported on-device (selected per call):
- **intensity** — isotropic (`intensity=:isotropic`, scaled by `rescale`) or the anisotropic
  `cloudIntensity` (`intensity=:cloud`, parameter `κ`).
- **velocity** — `vCircularCloud` (`velocity=:circular`) or the Pancoast
  `vCloudTurbulentEllipticalFlow` (`velocity=:turbulent`, parameters `σρᵣ, σρc, σΘᵣ, σΘc, θₑ, fEllipse,
  fFlow, σₜ`). The turbulent path draws extra normals/uniform per cloud from the same substream.

For physics beyond these built-ins, build on the host and call `gpu(m)`.
"""
function residentCloudModel(nClouds::Int, seed::Integer; μ::Real=500.0, β::Real=1.0, F::Real=0.5,
        rₛ::Real=1.0, θₒ::Real=π / 2, γ::Real=1.0, ξ::Real=1.0, i::Real=0.0, rescale::Real=1.0,
        ηₒ::Real=0.5, η₁::Real=0.5, αRM::Real=0.0, rNorm::Real=1.0,
        intensity::Symbol=:isotropic, κ::Real=0.0,
        velocity::Symbol=:circular, σρᵣ::Real=0.0, σρc::Real=0.0, σΘᵣ::Real=0.0, σΘc::Real=0.0,
        θₑ::Real=0.0, fEllipse::Real=0.0, fFlow::Real=0.0, σₜ::Real=0.0,
        backend=KernelAbstractions.CPU(), T=Float64)
    T <: Real || throw(ArgumentError("element type T must be <: Real, got $T"))
    intensity in (:isotropic, :cloud) ||
        throw(ArgumentError("intensity must be :isotropic or :cloud, got :$intensity"))
    velocity in (:circular, :turbulent) ||
        throw(ArgumentError("velocity must be :circular or :turbulent, got :$velocity"))
    ma = _build_cloud_modelarrays(nClouds, seed; μ=μ, β=β, F=F, rₛ=rₛ, θₒ=θₒ, γ=γ, ξ=ξ, i=i,
        rescale=rescale, ηₒ=ηₒ, η₁=η₁, αRM=αRM, rNorm=rNorm,
        useCloudI=(intensity === :cloud), κ=κ, useTurbulent=(velocity === :turbulent),
        σρᵣ=σρᵣ, σρc=σρc, σΘᵣ=σΘᵣ, σΘc=σΘc, θₑ=θₑ, fEllipse=fEllipse, fFlow=fFlow, σₜ=σₜ,
        backend=backend, T=T)
    return ResidentModel(ma, backend, 1, _rt_cloud_meta(ma, backend, T))
end

"""
    gpuCloudModel(nClouds, seed; T=Float32, kwargs...) -> ResidentModel

Build a cloud model directly on the GPU (device-resident [`ResidentModel`](@ref)) — the CUDA
counterpart of [`residentCloudModel`](@ref), defaulting to `Float32`. Requires CUDA.jl loaded.
"""
function gpuCloudModel end
