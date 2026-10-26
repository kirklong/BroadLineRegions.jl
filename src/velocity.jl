#!/usr/bin/env julia
using Random
"""    
    vCirc(r::Float64, rₛ::Float64=1.0)

Calculate circular velocity at radius `r` from central mass, with Schwarzschild radius `rₛ`.

Defaults to `rₛ=1.0` for unitless calculations.
"""
vCirc(r::Float64, rₛ::Float64=1.0) = √(rₛ/(2*r))

"""
    vCircularDisk(;r::Union{Float64,Vector{Float64}}, i::Float64, ϕ::Union{Vector{Float64},Float64}, rₛ=1.0, _...)

# Parameters
- `ϕ::Union{Vector{Float64},Float64}`: Azimuthal angle(s), in radians, ASCII alias=`phi`.
- `rₛ=1.0`: Scale radius (in terms of ``r_s``), ASCII alias=`r_s`.

Calculate line of sight velocity for circular orbit at radius `r` from central mass and inclined at angle `i` (rad) over grid of azimuthal angles `ϕ` (rad).
"""
@validated callback function vCircularDisk(;r::Union{Float64,Vector{Float64}}, i::Float64, ϕ::Union{Vector{Float64},Float64}, rₛ=1.0, _...)
    return @. -vCirc(r,$rₛ)*$sin(i)*sin(ϕ) #circular velocity where sides of disk are at ±π/2, flip so that right side is negative to match CM96
end

"""
    vCircularRadialDisk(;r::Union{Float64,Vector{Float64}}, i::Float64, ϕ::Union{Vector{Float64},Float64}, 
            vᵣFrac::Union{Vector{Float64},Float64}=0.0, inflow::Union{Vector{Bool},Bool}=true, rₛ=1.0, _...)

# Parameters
- `ϕ::Union{Vector{Float64},Float64}`: Azimuthal angle(s), in radians, ASCII alias=`phi`.
- `rₛ=1.0`: Scale radius (in terms of ``r_s``), ASCII alias=`r_s`.
- `vᵣFrac::Union{Vector{Float64},Float64}=0.0`: Radial velocity fraction, ASCII alias=`v_r_frac`.

Calculate line of sight velocity for circular orbit at radius `r` from central mass and inclined at angle `i` (rad) over grid of azimuthal angles `ϕ` (rad) with radial inflow/outflow.
"""
@validated callback function vCircularRadialDisk(;r::Union{Float64,Vector{Float64}}, i::Float64, ϕ::Union{Vector{Float64},Float64}, vᵣFrac::Union{Vector{Float64},Float64}=0.0, inflow::Union{Vector{Bool},Bool}=true, rₛ=1.0, _...)
    vsini = @. vCirc(r,$rₛ)*$sin(i)
    inflow = @. inflow ? -1.0 : 1.0 #at ϕ = 0 if inflow make v negative
    vᵣ = @. vsini*cos(ϕ)*vᵣFrac*inflow
    vϕ = @. vsini*sin(ϕ)*(1-vᵣFrac)
    return -1.0.*(vᵣ .+ vϕ) #flip so that right side is negative
end

"""
    vCircularCloud(;r::Float64, ϕ₀::Float64, i::Float64, rot::Float64, θₒ::Float64, rₛ::Float64=1.0, reflect::Bool=false, _...)

Calculate line of sight velocity for cloud in 3D space.

# Arguments
- `r::Float64`: radius from central mass (in terms of `rₛ`)
- `ϕ₀::Float64`: starting azimuthal angle in ring plane (rad), ASCII alias=`phi_0`.
- `i::Float64`: inclination angle of ring plane (rad)
- `rot::Float64`: rotation of system plane about z axis (rad)
- `θₒ::Float64`: opening angle of point, ASCII alias=`theta_o`.
- `rₛ::Float64=1.0`: Schwarzschild radius (optional, to convert to physical units), ASCII alias=`r_s`.
- `reflect::Bool=false`: whether the point is reflected across the midplane of the disk
- `_...`: extra kwargs, ignored

# Returns
- Line of sight velocity (`Float64`)
"""
@validated callback function vCircularCloud(;r::Float64, ϕ₀::Float64, i::Float64, rot::Float64, θₒ::Float64, rₛ::Float64=1.0, reflect::Bool=false, _...)
    v₀ = vCirc(r,rₛ)
    #proper circular-orbit tangent at ϕ₀: (-sin ϕ₀, +cos ϕ₀)*v₀ (perpendicular to the radius, r·v = 0),
    #then rotate into system coordinates and reflect if needed. The +cos ϕ₀ (vs -cos ϕ₀) is required so
    #the projected rotation sense survives the per-cloud `rot` and matches the disk (blue-shifted gas at +α).
    vXYZ = rotate3D_vector_scalar(-v₀*sin(ϕ₀), v₀*cos(ϕ₀), 0.0, i, rot, θₒ, reflect)
    return vXYZ[1] #line of sight velocity is x component after rotation (camera is at +x)
end

"""
    vCloudTurbulentEllipticalFlow(;σρᵣ::Float64, σρc::Float64, σΘᵣ::Float64, σΘc::Float64, 
            θₑ::Float64, fEllipse::Float64, fFlow::Float64, σₜ::Float64, r::Float64, 
            i::Float64, rot::Float64, θₒ::Float64, rₛ::Float64=1.0, ϕ₀::Float64=0.0, 
            reflect::Bool=false, rng::AbstractRNG=Random.GLOBAL_RNG, σρRelative::Bool=false, _...)

Calculate line of sight velocity for cloud in 3D space with potential for elliptical orbital velocities, in/outflow, and turbulence as in Pancoast+14.

# Arguments
- `σρᵣ::Float64`: Radial standard deviation around radial orbits, ASCII alias=`sigma_rho_r`.
- `σρc::Float64`: Radial standard deviation around circular orbits, ASCII alias=`sigma_rho_c`.
- `σΘᵣ::Float64`: Angular standard deviation around radial orbits, ASCII alias=`sigma_theta_r`.
- `σΘc::Float64`: Angular standard deviation around circular orbits, ASCII alias=`sigma_theta_c`.
- `θₑ::Float64`: Angle in v_ϕ-v_r plane, ASCII alias=`theta_e`.
- `fEllipse::Float64`: Fraction of elliptical orbits
- `fFlow::Float64`: If < 0.5, inflow, otherwise, outflow
- `σₜ::Float64`: Standard deviation of turbulent velocity, ASCII alias=`sigma_t`.
- `r::Float64`: Radius from central mass (in terms of `rₛ`)
- `i::Float64`: Inclination angle of ring plane (rad)
- `rot::Float64`: Rotation of system plane about z axis (rad)
- `θₒ::Float64`: Opening angle of point, ASCII alias=`theta_o`.
- `rₛ::Float64=1.0`: Schwarzschild radius (optional, to convert to physical units), ASCII alias=`r_s`.
- `ϕ₀::Float64=0.0`: Starting azimuthal angle in ring plane (rad), ASCII alias=`phi_0`.
- `reflect::Bool=false`: Whether the point is reflected across the midplane of the disk
- `rng::AbstractRNG=Random.GLOBAL_RNG`: Random number generator
- `σρRelative::Bool=false`: If true, `σρᵣ`/`σρc` are fractions of the local circular velocity (like `σₜ`) instead of absolute (units of c), ASCII alias=`sigma_rho_relative`.
- `_...`: Extra kwargs, ignored

# Returns
- Line of sight velocity (`Float64`)
"""
@validated callback function vCloudTurbulentEllipticalFlow(;σρᵣ::Float64,σρc::Float64, σΘᵣ::Float64, σΘc::Float64, θₑ::Float64, fEllipse::Float64, fFlow::Float64, σₜ::Float64,
    r::Float64, i::Float64, rot::Float64, θₒ::Float64, rₛ::Float64=1.0, ϕ₀::Float64=0.0, reflect::Bool=false, rng::AbstractRNG=Random.GLOBAL_RNG, σρRelative::Bool=false, _...)
    vc = vCirc(r,rₛ) #flip sign so that right side is negative
    vₜ = rand(rng,Normal(0.0,σₜ))*vc
    ρ = 0.0; Θ = 0.0
    if rand(rng) < fEllipse #elliptical orbit, distribution deviating from circular by σρc, σΘc, Pancoast 14 2.5.1
        ρ = rand(rng,Normal(vc,σρRelative ? σρc*vc : σρc))
        Θ = π/2 + rand(rng,Normal(0.0,σΘc))
    else #in/outflowing orbit, distribution deviating from circular by σρᵣ, σΘᵣ, Pancoast14 2.5.2
        ρ = rand(rng,Normal(vc,σρRelative ? σρᵣ*vc : σρᵣ))
        Θ = fFlow < 0.5 ? rand(rng,Normal(0.0,σΘᵣ)) + (π - θₑ) : rand(rng,Normal(0.0,σΘᵣ)) + θₑ
    end
    vx = √2*ρ*cos(Θ); vy = ρ*sin(Θ) #local (radial, tangential) components: radial along x, tangential along y at ϕ₀ = 0
    #rotate into the ring plane. The tangential (vy) term needs the CCW sign (+vy cos ϕ₀ in the 2nd
    #component); the radial (vx) terms are unchanged so the inflow/outflow convention is preserved
    #(near-side inflow ⟹ redshift, as in vCircularRadialDisk). Reduces to the corrected vCircularCloud
    #tangent in the circular limit (Θ = π/2 ⟹ (-sin ϕ₀, +cos ϕ₀)*vc).
    vXYZ = rotate3D_vector_scalar(-(vx*cos(ϕ₀)+vy*sin(ϕ₀)), vy*cos(ϕ₀)-vx*sin(ϕ₀), 0.0, i, rot, θₒ, reflect)
    return vXYZ[1]+vₜ #line of sight velocity is x component after rotation (camera is at +x), turbulence only along line of sight (see Pancoast14 2.5.3)
end

"""
    vElliptical(; a, i, ϕ, G=1.0, M=1.0, e=0.0)

# Parameters
- `ϕ::Union{Vector{Float64},Float64}`: Azimuthal angle(s), in radians, ASCII alias=`phi`.

Calculate line-of-sight velocities on an elliptical orbit with semimajor axis `a`,
eccentricity `e`, inclination `i` (radians), and azimuths `ϕ` (radians).
"""
@validated callback function vElliptical(;a::Float64, i::Float64, ϕ::Union{Vector{Float64},Float64}, G::Float64=1.0, M::Float64=1.0, e::Float64=0.0, _...)
    v² = @. G*M*((1 .+ 2*e.*cos.(ϕ) .+ e^2)./(a*(1 - e^2)))
    return @. sqrt.(v²).*sin(i).*cos.(ϕ)
end

"""
    rElliptical(; a, ϕ, e=0.0)

# Parameters
- `ϕ::Union{Vector{Float64},Float64}`: Azimuthal angle(s), in radians, ASCII alias=`phi`.

Calculate radii on an elliptical orbit with semimajor axis `a`, eccentricity `e`,
and azimuths `ϕ` (radians).
"""
@validated callback function rElliptical(;a::Float64, ϕ::Union{Vector{Float64},Float64}, e::Float64=0.0, _...)
    return @. (a*(1 - e^2))./(1 .+ e.*cos.(ϕ))
end