#!/usr/bin/env julia
"""
    DiskWind_I_(r::Vector{Float64}, ϕ::Vector{Float64}, i::Float64, f1::Float64, f2::Float64, f3::Float64, f4::Float64, α::Float64) 

Calculates the intensity of the disk wind at radius r from the central mass and inclined at angle i (rad) over grid of azimuthal angles ϕ (rad)
    Follows the prescription given in Long+ 2023 and 2025, similar to Chiang and Murray 1996 and 1997, assumes optically thick line emission limit

# Arguments:
- r::Vector{Float64} - radius from central mass (in terms of rₛ)
- ϕ::Vector{Float64} - list of azimuthal angles (rad)
- i::Float64 - inclination angle (rad)
- f1::Float64 - strength of radial velocity gradient dvᵣ/dr
- f2::Float64 - strength of Keplerian shear dvϕ/dr
- f3::Float64 - strength of velocity gradient dvθ/dr
- f4::Float64 - strength of velocity gradient dvθ/dθ (or isotropic emission)
- α::Float64 - power-law index of the source function S(r) ∝ r^(-α)
# Returns:
- intensity (arbitrary units) as a Vector{Float64}
"""
DiskWind_I_(r::Vector{Float64}, ϕ::Vector{Float64}, i::Float64, f1::Float64, f2::Float64, f3::Float64, f4::Float64, α::Float64) = begin #this is a performance hit
    sinϕ = sin.(ϕ); cosϕ = cos.(ϕ); sini = sin(i); cosi = cos(i)
    pre = @. sqrt(1 /(2 * r^3))
    term12 = @. (3*sini^2)*(cosϕ) * (√2*f1 * cosϕ + f2/2 *sinϕ)
    term3 = @. ((-f3*3*sini*cosi)) * cosϕ
    term4 = @. √2*f4*cosi^2 * $ones($length(ϕ))
    dvl_dl = @. pre * (term12 + term3 + term4)
    I = @. r^(-α) * abs(dvl_dl)
    I
end
"""
    DiskWind_I_(r::Float64, ϕ::Vector{Float64}, i::Float64, f1::Float64, f2::Float64, f3::Float64, f4::Float64, α::Float64) 

Same as `DiskWind_I_(r::Vector{Float64}, ϕ::Vector{Float64}, i::Float64, f1::Float64, f2::Float64, f3::Float64, f4::Float64, α::Float64)` but for a single radius `r` and a vector of azimuthal angles `ϕ`.
"""
DiskWind_I_(r::Float64, ϕ::Vector{Float64}, i::Float64, f1::Float64, f2::Float64, f3::Float64, f4::Float64, α::Float64) = begin #this is a performance hit
    sinϕ = sin.(ϕ); cosϕ = cos.(ϕ); sini = sin(i); cosi = cos(i)
    pre = sqrt(1 /(2 * r^3))
    term12 = @. (3*sini^2)*(cosϕ) * (√2*f1 * cosϕ + f2/2 *sinϕ)
    term3 = @. ((-f3*3*sini*cosi)) * cosϕ
    term4 = @. √2*f4*cosi^2 * $ones($length(ϕ))
    dvl_dl = @. pre * (term12 + term3 + term4)
    I = @. r^(-α) * abs(dvl_dl)
    I
end
"""
    DiskWind_I_(r::Vector{Float64}, ϕ::Float64, i::Float64, f1::Float64, f2::Float64, f3::Float64, f4::Float64, α::Float64)

Same as `DiskWind_I_(r::Vector{Float64}, ϕ::Vector{Float64}, i::Float64, f1::Float64, f2::Float64, f3::Float64, f4::Float64, α::Float64)` but for a vector of radii `r` and a single azimuthal angle `ϕ`.
"""
DiskWind_I_(r::Vector{Float64}, ϕ::Float64, i::Float64, f1::Float64, f2::Float64, f3::Float64, f4::Float64, α::Float64) = begin #this is a performance hit
    sinϕ = sin(ϕ); cosϕ = cos(ϕ); sini = sin(i); cosi = cos(i)
    pre = @. sqrt(1 /(2 * r^3))
    term12 = (3*sini^2)*(cosϕ) * (√2*f1 * cosϕ + f2/2 *sinϕ)
    term3 = ((-f3*3*sini*cosi)) * cosϕ
    term4 = √2*f4*cosi^2 
    dvl_dl = @. pre * (term12 + term3 + term4)
    I = @. r^(-α) * abs(dvl_dl)
    I
end
"""
    DiskWind_I_(r::Float64, ϕ::Float64, i::Float64, f1::Float64, f2::Float64, f3::Float64, f4::Float64, α::Float64)

Same as `DiskWind_I_(r::Vector{Float64}, ϕ::Vector{Float64}, i::Float64, f1::Float64, f2::Float64, f3::Float64, f4::Float64, α::Float64)` but for a single radius `r` and a single azimuthal angle `ϕ`.
"""
DiskWind_I_(r::Float64, ϕ::Float64, i::Float64, f1::Float64, f2::Float64, f3::Float64, f4::Float64, α::Float64) = begin #this is a performance hit
    sinϕ = sin(ϕ); cosϕ = cos(ϕ); sini = sin(i); cosi = cos(i)
    pre = sqrt(1 /(2 * r^3))
    term12 = (3*sini^2)*(cosϕ) * (√2*f1 * cosϕ + f2/2 *sinϕ)
    term3 = ((-f3*3*sini*cosi)) * cosϕ
    term4 = √2*f4*cosi^2 
    dvl_dl = pre * (term12 + term3 + term4)
    I = r^(-α) * abs(dvl_dl)
    I
end
"""
    DiskWindIntensity(;r::Union{Vector{Float64},Float64}, i::Float64, ϕ::Union{Vector{Float64},Float64}, 
            f1::Float64, f2::Float64, f3::Float64, f4::Float64, α::Float64, rMin::Float64 = 0.0, rMax::Float64 = Inf, _...)

Calculates the intensity from a disk-wind model of the BLR following the prescription given in Long+ 2023 and 2025, similar to Chiang and Murray 1996 and 1997, assumes optically thick line emission limit (Sobolev).

# Arguments:
- r::Union{Vector{Float64},Float64} - radius from central mass (in terms of rₛ)
- i::Float64 - inclination angle (rad)
- ϕ::Union{Vector{Float64},Float64} - list of azimuthal angles (rad)
- f1::Float64 - strength of radial velocity gradient ``\\mathrm{d}v_{r}/\\mathrm{d}r``
- f2::Float64 - strength of Keplerian shear velocity gradient ``\\mathrm{d}v_{\\phi}/\\mathrm{d}r``
- f3::Float64 - strength of radial lifting velocity gradient ``\\mathrm{d}v_{\\theta}/\\mathrm{d}r``
- f4::Float64 - strength of vertical lifting velocity gradient ``\\mathrm{d}v_{\\theta}/\\mathrm{d}\\theta`` (equivalent to isotropic emission)
- α::Float64 - power-law index of the source function ``S(r) \\propto r^{-\\alpha}``
- rMin::Float64 - minimum radius to consider (default: 0.0)
- rMax::Float64 - maximum radius to consider (default: Inf)

# Returns:
- intensity (arbitrary units) as a Vector{Float64}
"""
function DiskWindIntensity(;r::Union{Vector{Float64},Float64}, i::Float64, ϕ::Union{Vector{Float64},Float64}, f1::Float64, f2::Float64, f3::Float64, f4::Float64, α::Float64, rMin::Float64 = 0.0, rMax::Float64 = Inf, _...)
    """calculate the intensity of the disk wind at radius r from the central mass and inclined at angle i (rad) over grid of azimuthal angles ϕ (rad)
        Follows the prescription given in Long+ 2023 and 2024, similar to Chiang and Murray 1996, assumes optically thick line emission limit
        params: 
            r: radius from central mass (in terms of rₛ) Union{Vector{Float64},Float64}
            i: inclination angle (rad) {Float64}
            ϕ: list of azimuthal angles (rad) Union{Vector{Float64},Float64}
            f1: strength of radial velocity gradient dvᵣ/dr {Float64}
            f2: strength of Keplerian shear dvϕ/dr {Float64}
            f3: strength of velocity gradient dvθ/dr {Float64}
            f4: strength of velocity gradient dvθ/dθ (or isoptropic emission) {Float64}
            α: power-law index of the source function S(r) ∝ r^(-α) {Float64}
            _: extra kwargs, ignored
        returns:
            intensity (arbitrary units) {Vector{Float64}}
    """
    if typeof(r) == Float64 || typeof(ϕ) == Float64
        I = (r >= rMin && r <= rMax) ? DiskWind_I_(r, ϕ, i, f1, f2, f3, f4, α) : 0.0
    elseif typeof(r) == Vector{Float64} && typeof(ϕ) == Vector{Float64}
        I = vcat([DiskWind_I_(ri,ϕi,i,f1,f2,f3,f4,α) for (ri,ϕi) in zip(r,ϕ)]...) #one intensity value at each r,ϕ pair
        if round(minimum(r),sigdigits=9) <= round(rMin,sigdigits=9) || round(maximum(r),sigdigits=9) >= round(rMax,sigdigits=9)
            I = [(round(ri,sigdigits=9) >= round(rMin,sigdigits=9) && round(ri,sigdigits=9) <= round(rMax,sigdigits=9)) ? Ii : 0.0 for (ri,Ii) in zip(r,I)]
        end
    else
        error("got unsupported types for r and ϕ: got $(typeof(r)) and $(typeof(ϕ)), expected Float64 or Vector{Float64}")
    end
    return I
end

"""
    IsotropicIntensity(;r::Union{Vector{Float64},Float64}, ϕ::Union{Vector{Float64},Float64}, 
            rescale::Float64=1.0, rMin::Float64 = 0.0, rMax::Float64 = Inf, _...)

Returns a constant intensity value at all radii and azimuthal angles, rescaled by `rescale` factor.
"""
function IsotropicIntensity(;r::Union{Vector{Float64},Float64}, ϕ::Union{Vector{Float64},Float64}, rescale::Float64=1.0, rMin::Float64 = 0.0, rMax::Float64 = Inf, _...)
    """calculate the intensity of the isotropic emission over grid of azimuthal angles ϕ (rad) or at one azimuthal angle
    params: 
        r: radius from central mass (in terms of rₛ) Union{Vector{Float64},Float64}
        ϕ: list of azimuthal angles or single angle (rad) {Vector{Float64}}
        rescale: rescale factor for the intensity {Float64}
        _: extra kwargs, ignored
    returns:
        intensity (arbitrary units) Union{Float64,Vector{Float64},Matrix{Float64}}
    """
    if typeof(ϕ) == Float64 && typeof(r) == Float64
        I = (round(r,sigdigits=9) >= round(rMin,sigdigits=9) && round(r,sigdigits=9) <= round(rMax,sigdigits=9)) ? rescale : 0.0
    else
        if typeof(r) == Float64
            I = (round(r,sigdigits=9) >= round(rMin,sigdigits=9) && round(r,sigdigits=9) <= round(rMax,sigdigits=9)) ? rescale.*ones(length(ϕ)) : zeros(length(ϕ))
        elseif typeof(ϕ) == Vector{Float64}
            I = vcat([round(ri,sigdigits=9) >= round(rMin,sigdigits=9) && round(ri,sigdigits=9) <= round(rMax,sigdigits=9) ? rescale : 0.0 for (ri,ϕi) in zip(r,ϕ)]...) #one intensity value at each r,ϕ pair
        else
            error("got unsupported types for r and ϕ: got $(typeof(r)) and $(typeof(ϕ)), expected Float64 or Vector{Float64}")
        end
    end
    return I
end    

W(ϕ::Float64,κ::Float64) = 1/2+κ*cos(ϕ)

"""
    cloudIntensity(;r::Float64, ϕ::Float64, θₒ::Float64, ϕ₀::Float64, rot::Float64, i::Float64, κ::Float64=0.0, _...)

Calculate the intensity of the cloud at radius `r` from the central mass and inclined at angle `i` (rad) over grid of azimuthal angles `ϕ` (rad) following Pancoast+ 2011 and 2014.

# Arguments
- `r::Float64`: radius from central mass (in terms of r₍ₛ₎)
- `ϕ::Float64`: azimuthal angle (rad)
- `θₒ::Float64`: opening angle of cloud
- `ϕ₀::Float64`: initial azimuthal angle
- `rot::Float64`: rotation angle
- `i::Float64`: inclination angle
- `κ::Float64`: anisotropy parameter

# Returns
- Intensity (arbitrary units) as a `Float64`
"""
function cloudIntensity(;r::Float64,ϕ::Float64,θₒ::Float64,ϕ₀::Float64,rot::Float64,i::Float64,κ::Float64=0.0, _...)
    xyzSys = rotate3D(r,ϕ₀,i,rot,θₒ) #system coordinates xyz
    ϕw = atan(xyzSys[2],xyzSys[1]) #angle between cloud and BH line of sight
    I = W(ϕw,κ) #clouds don't have scaling -- the scaling is set by geometry itself, just whether they are on or off 
    return I
end

"""
    IϕCloudMask(;r::Float64, ϕ::Float64, θₒ::Float64, ϕ₀::Float64, rot::Float64, i::Float64, κ::Float64=0.0, 
            ϕMin::Float64, ϕMax::Float64, overdense::Bool=false, _...)

Calculates the intensity using `cloudIntensity` but with an extra layer of masking based on azimuthal angle `ϕ` and the specified range `[ϕMin, ϕMax]`. If `ϕ` is outside this range, the intensity is set to 0.0.
If overdense is `true`, the intensity is multiplied by 2.0 to account for the "lost" cloud and thus preserve total intensity in the model.
"""
function IϕCloudMask(;r::Float64,ϕ::Float64,θₒ::Float64,ϕ₀::Float64,rot::Float64,i::Float64,κ::Float64=0.0,ϕMin::Float64,ϕMax::Float64,overdense::Bool=false, _...)
    if ϕ > π || ϕ < -π
        ϕ = atan(sin(ϕ),cos(ϕ))
    end
    if ϕMin > π || ϕMin < -π
        ϕMin = atan(sin(ϕMin),cos(ϕMin))
    end
    if ϕMax > π || ϕMax < -π
        ϕMax = atan(sin(ϕMax),cos(ϕMax))
    end
    I = 0.0
    if ϕ >= ϕMin && ϕ <= ϕMax
        I = cloudIntensity(r=r,ϕ=ϕ,θₒ=θₒ,ϕ₀=ϕ₀,rot=rot,i=i,κ=κ)
        if overdense
            I *= 2.0 #x2 to replace cloud "lost" 
        end
    end
    return I 
end

"""
    IϕDiskWindMask(;r::Union{Vector{Float64},Float64}, ϕ::Union{Vector{Float64},Float64}, i::Float64, 
            f1::Float64, f2::Float64, f3::Float64, f4::Float64, α::Float64, ϕMin::Float64, ϕMax::Float64, _...)

Calculates the intensity using `DiskWindIntensity` but with an extra layer of masking based on azimuthal angle `ϕ` and the specified range `[ϕMin, ϕMax]`. If `ϕ` is outside this range, the intensity is set to 0.0.
"""
function IϕDiskWindMask(;r::Union{Vector{Float64},Float64},ϕ::Union{Vector{Float64},Float64},i::Float64,f1::Float64,f2::Float64,f3::Float64,f4::Float64,α::Float64,ϕMin::Float64,ϕMax::Float64, _...)
    I = nothing
    if typeof(ϕ) == Float64
        if ϕ > π || ϕ < -π
            ϕ = atan(sin(ϕ),cos(ϕ))
        end
        if ϕMin > π || ϕMin < -π
            ϕMin = atan(sin(ϕMin),cos(ϕMin))
        end
        if ϕMax > π || ϕMax < -π
            ϕMax = atan(sin(ϕMax),cos(ϕMax))
        end
        I = (ϕ >= ϕMin && ϕ <= ϕMax) ? DiskWindIntensity(r=r,i=i,ϕ=ϕ,f1=f1,f2=f2,f3=f3,f4=f4,α=α) : 0.0 #left side of disk goes 0 -> -π, right side goes 0 -> π, where 0 (top) is tilted towards the observer (exclude back of disk)
    else
        I = vcat([IϕDiskWindMask(r=ri,ϕ=ϕi,i=i,f1=f1,f2=f2,f3=f3,f4=f4,α=α,ϕMin=ϕMin,ϕMax=ϕMax) for (ri,ϕi) in zip(r,ϕ)]...)
    end
    return I
end