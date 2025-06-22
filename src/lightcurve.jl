#!/usr/bin/env julia
"""
    function DRW(;t::Array{Float64}=collect(range(0,stop=100,length=1001)), μ::Float64=0.0, τ::Float64=1.0, σ::Float64=0.1)
    
Generate a damped random walk continuum lightcurve following [MacLeod+2010](https://iopscience.iop.org/article/10.1088/0004-637X/721/2/1014/pdf)
# Arguments:
-`t::Vector{Float64}``: time array that returned lightcurve will be sampled at
-`μ::Float64`: mean value of random walk
-`τ::Float64`: characteristic time scale of random walk
-`σ::Float64`: standard deviation of random walk
returns:
-`C::Vector{Float64}`: damped random walk array corresponding to times t
"""
function DRW(;t::Array{Float64}=collect(range(0,stop=100,length=1001)), μ::Float64=0.0, τ::Float64=1.0, σ::Float64=0.1)
    if nt == 0
        nt = length(t)
    end
    SF∞ = √2*σ #eq 4 of https://iopscience.iop.org/article/10.1088/0004-637X/721/2/1014/pdf
    C = zeros(nt)
    C[1] = μ + rand(Normal(0,σ)) #first point
    for i in 2:length(t)
        Δt = t[i] - t[i-1]
        #below from eq 5 of https://iopscience.iop.org/article/10.1088/0004-637X/721/2/1014/pdf
        E = exp(-Δt/τ)*C[i-1]+μ*(1-exp(-Δt/τ))
        V = 0.5*SF∞^2*(1-exp(-2*Δt/τ))
        C[i] = E + rand(Normal(0,sqrt(V)))
    end
    return C
end

"""
    function syntheticLC(Ct,CLC,Ψτ;tStart=0.0,normalize=false)
Generate a synthetic lightcurve from a continuum lightcurve and a 1D transfer function.
# Arguments:
- `Ct::Vector{Float64}`: time array for the continuum lightcurve
- `CLC::Vector{Float64}`: continuum lightcurve values at times `Ct`
- `Ψτ::Vector{Float64}`: transfer function values at times `Ct`
- `tStart::Float64`: time after which the lightcurve should be normalized to zero
- `normalize::Bool`: whether to normalize the lightcurve to the extent of variations in the continuum
# Returns:
- `Δlc::Vector{Float64}`: synthetic lightcurve sampled at times `Ct` with zero point ~at `tStart` and optionally normalized
"""
function syntheticLC(Ct,CLC,Ψτ;tStart=0.0,normalize=false) #assumes Ψτ is already interpolated to continuum timestamps as in functions.jl/getProfiles
    t = Ct
    if length(Ψτ) < length(t)
        Ψτ = vcat(Ψτ, zeros(length(t) - length(Ψτ))) #extend Ψτ to match t length
    elseif length(Ψτ) > length(t)
        Ψτ = Ψτ[1:length(t)] #truncate Ψτ to match t length
    end
    C = C./maximum(C)
    span = maximum(C) - minimum(C) 
    ΔC = C .- C[1]
    ΔC = ΔC ./ span #make span 1, normalize to first point so that + = brighter and - = dimmer
    function LC(t,Ψτ,ΔC) #ΔLC variations at times t following formula for convolution of transfer function with continuum variations
        N = length(t)
        LC = zeros(N)
        for ti=1:N
            integrand = 0
            for τ = 1:ti
                dτ = ti < N ? t[τ+1]-t[τ] : t[end] - t[end-1]
                integrand += Ψτ[τ]*ΔC[ti-τ+1]*dτ
            end
            LC[ti] = integrand
        end
        return LC
    end
    Δlc = LC(t,Ψτ,ΔC) #these are variations, not properly normalized
    zeroPoint = findfirst(t.>=tStart) #find first index after tStart, this is value that should be "zero" 
    Δlc = Δlc .- Δlc[zeroPoint] #normalize to zero point
    norm = normalize ? maximum(Δlc) - minimum(Δlc) : 1.0 #when normalized the extent of variations in continuum should be the same as in emission line
    return Δlc./norm
end