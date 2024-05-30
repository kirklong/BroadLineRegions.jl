#!/usr/bin/env julia

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
    ϕ′ = ϕ .+ π/2 # shift ϕ by 90° to match convention in Long+ 2023 and Waters+ 2016
    I = nothing
    I_(r, ϕ, f1, f2, f3, f4, α) = begin 
        sinϕ = sin.(ϕ); cosϕ = cos.(ϕ); sini = sin(i); cosi = cos(i)
        pre = @. sqrt(1 /(2 * r^3))
        term12 = @. (3*sini^2)*(cosϕ) * (√2*f1 * cosϕ + f2/2 *sinϕ)
        term3 = @. ((-f3*3*sini*cosi)) * cosϕ
        term4 = @. √2*f4*cosi^2 * $ones($length(ϕ))
        dvl_dl = @. pre * (term12 + term3 + term4)
        I = @. r^(-α) * abs(dvl_dl)
        I = (typeof(r) == Float64 && typeof(ϕ) == Float64) ? I[1] : I
    end
    if typeof(r) == Float64 || typeof(ϕ) == Float64
        I = (r >= rMin && r <= rMax) ? I_(r, ϕ′, f1, f2, f3, f4, α) : 0.0
    elseif typeof(r) == Vector{Float64} && typeof(ϕ) == Vector{Float64}
        I = vcat([I_(ri,ϕi,f1,f2,f3,f4,α) for (ri,ϕi) in zip(r,ϕ′)]...) #one intensity value at each r,ϕ pair
        if minimum(r) < rMin || maximum(r) > rMax
            I = [(ri >= rMin && ri <= rMax) ? Ii : 0.0 for (ri,Ii) in zip(r,I)]
        end
    else
        error("got unsupported types for r and ϕ: got $(typeof(r)) and $(typeof(ϕ)), expected Float64 or Vector{Float64}")
    end
    return I
end

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
        I = (r >= rMin && r <= rMax) ? rescale : 0.0
    else
        if typeof(r) == Float64
            I = (r >= rMin && r <= rMax) ? rescale.*ones(length(ϕ)) : zeros(length(ϕ))
        elseif typeof(ϕ) == Vector{Float64}
            I = vcat([(ri >= rMin && ri <= rMax) ? rescale : 0.0 for (ri,ϕi) in zip(r,ϕ)]...) #one intensity value at each r,ϕ pair
        else
            error("got unsupported types for r and ϕ: got $(typeof(r)) and $(typeof(ϕ)), expected Float64 or Vector{Float64}")
        end
    end
    return I
end    