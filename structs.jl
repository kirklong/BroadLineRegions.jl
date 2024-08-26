#!/usr/bin/env julia
using Random

struct ring{V,F}
    """structure to hold parameters of model ring
    attributes:
        r: distance from central mass (in terms of rₛ) Union{Vector{Float64},Float64
            - can be Float64 for constant r across ϕ or single point at (r,ϕ), otherwise it is a Vector{Float64} of distances corresponding to azimuthal angles ϕ
            - user can (optionally) supply a function to calculate r from other parameters (see constructor) that returns a Vector{Float64} or Float64
        i: inclination angle (rad) {Float64}
            - must be between 0 and π/2, with 0 being face-on and π/2 being edge-on
        rot: rotation of system plane about z axis (rad) {Float64}
        θₒ: opening angle (rad) of ring {Float64}
            - should be between 0 and π/2
            - optional, defaults to 0
        v: line of sight velocity {Union{Vector{Float64},Float64}}
            - Float64 for a single angle, otherwise a Vector{Float64} of velocities corresponding to azimuthal angles ϕ
            - user can (optionally) supply a function to calculate v from other parameters (see constructor) that returns a Vector{Float64}
        I: intensity {Union{Vector{Float64},Float64,Matrix{Float64}}
            - Float64 for a single angle, otherwise a Vector{Float64} of intensities corresponding to azimuthal angles ϕ
            - user can (optionally) supply a function to calculate I from other parameters (see constructor) that returns a Vector{Float64}
        ϕ: azimuthal angle (rad) {Union{Vector{Float64},Float64}}
            - Float64 for a single point, otherwise a Vector{Float64} of azimuthal angles
        ϕₒ (optional): azimuthal angle of system before rotation (rad) {Union{Vector{Float64},Float64}}
            - Float64 for a single point, otherwise a Vector{Float64} of azimuthal angles
            - defaults to 0.0 if not provided or if rot is 0.0 (if no rotation ϕ = ϕₒ)
        ΔA: projected area of ring in image (used in calculating profiles) {Union{Vector{Float64},Float64}}
            - Float64 for a single point, otherwise a Vector{Float64} of projected areas corresponding to azimuthal angles ϕ
            - defaults to 1.0 if not provided
    
    constructor:
        ring(;r, i, rot, θₒ, v, I, ϕ, ϕₒ, kwargs...)
        r: Float64, Vector{Float64} or Function
            - if function, must return a Vector{Float64} or Float64 corresponding to ϕ
        i: Float64
            - must be between 0 and π/2
        rot: Float64
            - optional (defaults to 0.0)
            - must be between 0 and 2π (or -π to π)
            - describes 3D rotation about z axis of ring plane
        θₒ: Float64
            - optional (defaults to 0.0)
            - must be between 0 and π/2
            - opening angle of ring (i.e. thin disk has opening angle near 0.0 and spherical shell has opening angles distributed between 0 and π/2)
        v: Float64, Vector{Float64} or Function
            - if function, must return a Vector{Float64} corresponding to ϕ
        I: Float64, Vector{Float64}, Matrix{Float64}, or Function
            - if function, must return a Vector{Float64} corresponding to ϕ for each r (if r is a Vector{Float64} must return a Matrix{Float64} of intensities for each r)
        ϕ: Float64 or Vector{Float64}
            - must be between 0 and 2π
            - if Float64, r, v, and I must be Float64 (or functions that return Float64)
        ϕₒ: Float64 or Vector{Float64}
            - optional (defaults to 0.0)
        ΔA: Float64 or Vector{Float64}
            - projected area of ring in image (used in calculating profiles)
            - i.e. for log-spherical ring ΔA = r^2*Δr*Δϕ, for linear-spherical ring ΔA = r*Δr*Δϕ, for cloud ΔA could be size of cloud (note here r would be the image r not the physical r)
            - defaults to 1.0 if not provided
        kwargs: contain extra keyword arguments for v, I, and/or r if they are functions (see examples)
    """

    r::Union{V,F,Function}
    i::F
    rot::F
    θₒ::F
    v::Union{V,F,Function}
    I::Union{V,F,Matrix{F},Function}
    ϕ::Union{V,F}
    ϕₒ::Union{V,F}
    ΔA::Union{V,F}
    reflect::Bool

    function ring(;kwargs...) #could re-write this to use multiple dispatch? i.e. ring(;r::Float64, i::Float64, e::Float64, v::Float64, I::Float64, ϕ::Float64) etc.
        """
        constructor for ring struct -- takes in kwargs (detailed above) and returns a ring object (detailed above) while checking for errors
        """
        r = nothing; i = nothing; v = nothing; I = nothing; ϕ = nothing; ΔA = 1.0; rot = 0.0; θₒ = 0.0; ϕₒ = 0.0; reflect = false
        try; r = kwargs[:r]; catch; error("r must be provided as kwarg"); end
        try; i = kwargs[:i]; catch; error("i must be provided as kwarg"); end
        try; v = kwargs[:v]; catch; error("v must be provided as kwarg"); end
        try; I = kwargs[:I]; catch; error("I must be provided as kwarg"); end
        try; ϕ = kwargs[:ϕ]; catch; error("ϕ must be provided as kwarg"); end
        try; ΔA = kwargs[:ΔA]; catch; error("ΔA not provided, defaulting to 1.0"); end
        try; rot = kwargs[:rot]; catch; println("rot not provided: defaulting to 0.0"); end
        try; θₒ = kwargs[:θₒ]; catch; println("θₒ not provided: defaulting to 0.0"); end
        try; ϕₒ = kwargs[:ϕₒ]; catch; println("ϕₒ not provided: defaulting to 0.0"); end
        try; reflect = kwargs[:reflect]; catch; println("reflect not provided: defaulting to false"); end
        kwargs = values(kwargs)
        
        @assert typeof(reflect) == Bool "reflect must be Bool, got $(typeof(reflect))"
        @assert typeof(i) == Float64 "i must be Float64, got $(typeof(i))"
        @assert (i >= -π/2) && (i <= π/2) "i must be between -π/2 and π/2, got $i"
        i = -i #flip i to match convention of bottom of disk being tilted "towards" observer, relic of how rotation matrix was implemented.
        @assert (typeof(rot) == Float64) "rot must be Float64, got $(typeof(rot))"
        @assert (typeof(θₒ) == Float64) "θₒ must be Float64, got $(typeof(θₒ))"
        @assert (θₒ >= 0) && (θₒ <= π/2) "θₒ must be between 0 and π/2, got $θₒ"
        @assert (typeof(ΔA) == Float64) || (typeof(ΔA) == Vector{Float64}) "ΔA must be Float64 or Vector{Float64}, got $(typeof(ΔA))"
        if typeof(ΔA) == Vector{Float64}
            @assert length(ΔA) == length(ϕ) "ΔA must be the same length as ϕ"
        end
        @assert (typeof(r) == Vector{Float64}) || (typeof(r) == Float64) || (isa(r,Function)) "r must be Float64, Vector{Float64} or Function, got $(typeof(v))"        

        @assert (typeof(ϕ) == Vector{Float64}) || (typeof(ϕ) == Float64) "ϕ must be Float64 or Vector{Float64}, got $(typeof(ϕ))"
        @assert (typeof(ϕₒ) == Vector{Float64}) || (typeof(ϕₒ) == Float64) "ϕₒ must be Float64 or Vector{Float64}, got $(typeof(ϕₒ))"

        @assert (typeof(v) == Vector{Float64}) || (typeof(v) == Float64) || (isa(v,Function)) "v must be Float64, Vector{Float64} or Function, got $(typeof(v))"
        if isa(v,Function)
            try; v = v(;kwargs...); catch; error("error in function v -- check kwargs"); end
            @assert length(v) == length(ϕ) "v must return a vector of the same length as ϕ"
        elseif typeof(v) == Vector{Float64}
            @assert length(v) == length(ϕ) "v must be the same length as ϕ"
        end
        @assert typeof(v) == typeof(ϕ) "v and ϕ must be the same type"

        @assert (typeof(I) == Vector{Float64}) || (typeof(I) == Float64) || (typeof(I) == Matrix{Float64}) || (isa(I,Function)) "I must be Float64, Vector{Float64}, Matrix{Float64}, or Function, got $(typeof(I))"
        if isa(I,Function)
            try; I = I(;kwargs...); catch; error("error in function I -- check kwargs"); end
            if length(r) == 1
                @assert length(I) == length(ϕ) "I must return a vector of the same length as ϕ"
                @assert typeof(I) == typeof(ϕ) "I and ϕ must be the same type"
            elseif typeof(I) == Matrix{Float64}
                @assert size(I) == (length(r),length(ϕ)) "I must return a matrix of size (length(r),length(ϕ))"
            elseif typeof(I) == Vector{Float64}
                @assert length(I) == length(ϕ) && length(I) == length(r) "I must be the same length as ϕ and r"
            else
                error("Invalid return type for function I (got $(typeof(I))): must return Vector{Float64} with length(ϕ) = length(r) or Matrix{Float64} with size (length(r),length(ϕ))")
            end
        elseif typeof(I) == Vector{Float64}
            @assert length(I) == length(ϕ) "I must be the same length as ϕ"
        elseif typeof(I) == Matrix{Float64}
            @assert size(I) == (length(r),length(ϕ)) "I must be a matrix of size (length(r),length(ϕ))"
        end

        new{Vector{Float64},Float64}(r,i,rot,θₒ,v,I,ϕ,ϕₒ,ΔA,reflect)
    end
end

Base.show(io::IO, r::ring) = begin
    println(io, "ring struct with inclination $(round(r.i,sigdigits=3)) rad, rotation $(round(r.rot,sigdigits=3)) rad, and opening angle $(round(r.θₒ,sigdigits=3)) rad")
    xMin = nothing; xMax = nothing
    if typeof(r.ϕ) == Float64
        if r.ϕ != r.ϕₒ
            println(io, "cloud with final azimuthal angle: $(round(r.ϕ,sigdigits=3)) rad and initial azimuthal angle: $(round(r.ϕₒ,sigdigits=3)) rad")
        else
            println(io, "cloud with final azimuthal angle: $(round(r.ϕ,sigdigits=3)) rad")
        end
        if r.reflect
            println("\t--cloud originally on back side of disk, reflected across disk mid-plane to front")
        end
    else
        try 
            xMin = minimum(i for i in r.ϕ if !isnan(i))
            xMax = maximum(i for i in r.ϕ if !isnan(i)) 
        catch 
            xMin = NaN
            xMax = NaN 
        end
        println(io, "intensity distributed over $(length(r.ϕ)) azimuthal angles ($(round(xMin,sigdigits=3)) < ϕ < $(round(xMax,sigdigits=3)) rad)")
    end
    if typeof(r.r) == Float64
        println(io, "cloud at radius: $(r.r) rₛ")
    else
        try 
            xMin = minimum(i for i in r.r if !isnan(i))
            xMax = maximum(i for i in r.r if !isnan(i)) 
        catch 
            xMin = NaN
            xMax = NaN 
        end
        println(io, "intensity distributed over $(length(r.r)) radii ($(round(xMin,sigdigits=3)) < r < $(round(xMax,sigdigits=3)) rₛ)")
    end
    if typeof(r.v) == Float64
        println(io, "cloud line of sight velocity: $(round(r.v,sigdigits=3)) c")
    else
        try 
            xMin = minimum(i for i in r.v if !isnan(i))
            xMax = maximum(i for i in r.v if !isnan(i)) 
        catch 
            xMin = NaN
            xMax = NaN 
        end
        println(io, "line of sight velocity: $(round(xMin,sigdigits=3)) < v < $(round(xMax,sigdigits=3)) c")
    end
    if typeof(r.I) == Float64
        println(io, "cloud intensity: $(round(r.I,sigdigits=3)) arbitrary units")
    else
        try
            xMin = minimum(i for i in r.I if !isnan(i))
            xMax = maximum(i for i in r.I if !isnan(i))
        catch
            xMin = NaN
            xMax = NaN
        end
        println(io, "intensity: $(round(xMin,sigdigits=3)) < I < $(round(xMax,sigdigits=3)) arbitrary units")
    end
    if typeof(r.ΔA) == Float64
        println(io, "projected area: $(round(r.ΔA,sigdigits=3)) rₛ²")
    else
        try
            xMin = minimum(i for i in r.ΔA if !isnan(i))
            xMax = maximum(i for i in r.ΔA if !isnan(i))
        catch
            xMin = NaN
            xMax = NaN
        end
        println(io, "projected area: $(round(xMin,sigdigits=3)) < ΔA < $(round(xMax,sigdigits=3)) rₛ²")
    end
end

@kwdef struct profile
    """profile struct to hold binned data, usually set to model struct with profiles.jl#setProfile
    attributes:
        name: name of profile {Symbol}
        centers: bin centers {Vector{Float64}}
        edges: bin edges {Vector{Float64}}
        sums: sum of values in each bin {Vector{Float64}}
    """
    name::Symbol
    binCenters::Vector{Float64}
    binEdges::Vector{Float64}
    binSums::Vector{Float64}
end

Base.show(io::IO, p::profile) = begin
    println(io, "$(p.name) profile struct with $(length(p.binCenters)) bins")
end

struct camera #need to modify to include "imgs" of each quantity -- most importantly v and I after raytracing
    """camera coordinates struct
    attributes:
        name: name of image {Symbol}
        x: x values Union{Vector{Float64}, Matrix{Float64}} 
        y: y values Union{Vector{Float64}, Matrix{Float64}
    """
    α::Union{Vector{Float64},Matrix{Float64}}
    β::Union{Vector{Float64},Matrix{Float64}}
end

Base.show(io::IO, c::camera) = begin
    pix = length(c.α)
    println(io, "camera with $pix pixels and range: $(round(minimum(c.α),sigdigits=3)) < α < $(round(maximum(c.α),sigdigits=3)) and $(round(minimum(c.β),sigdigits=3)) < β < $(round(maximum(c.β),sigdigits=3))")
end 

meshgrid(x,y) = (reshape(repeat(x,outer=length(y)),length(x),length(y)), reshape(repeat(y,inner=length(x)),length(x),length(y)))

mutable struct model #restructure to make mutable struct -- Δr doesn't need to be kept track of it only matters for binning profiles? profiles should be a dict i.e. LP => array, delayProfile => array, custom => array
    #make immutable to see if improves performance cost
    #add image as keyword arg to constructors -- only initialize if true to save on performance
    
    """structure to hold many rings and their parameters
    attributes:
        rings: Vector{ring}
            - list of ring objects
    """
    rings::Vector{ring}
    profiles::Union{Nothing,Dict{Symbol,profile}}
    camera::Union{Nothing,camera}

    function model(rings::Vector{ring{Vector{Float64},Float64}})
        """
        constructor for model struct -- takes in rings and returns a model object (detailed above) while checking for errors
        """
        r = [ring.r for ring in rings]; ϕₒ = [ring.ϕₒ for ring in rings]; i = [ring.i for ring in rings]; rot = [ring.rot for ring in rings]; θₒ = [ring.θₒ for ring in rings]; reflect = [ring.reflect for ring in rings]
        α = zeros(length(r)); β = zeros(length(r))
        for (i,(ri,ϕi,ii,roti,θₒi,reflecti)) in enumerate(zip(r,ϕₒ,i,rot,θₒ,reflect))
            α[i], β[i] = photograph(ri,ϕi,ii,roti,θₒi,reflecti)  #get camera coordinates from physical 
        end
        new(rings,Dict{Symbol,profile}(),camera(stack(α,dims=1),stack(β,dims=1)))
    end

    function model(rMin::Float64, rMax::Float64, i::Float64, nr::Int, nϕ::Int, I::Function, v::Function, scale::Symbol; kwargs...)
        """constructor for model struct -- takes in rMin, rMax, i, nr, nϕ, I, v, and scale and returns a model object (detailed above) while checking for errors
        params:
            rMin: minimum radius of model (in terms of rₛ) {Float64}
            rMax: maximum radius of model (in terms of rₛ) {Float64}
            i: inclination angle (rad) {Float64} (all rings have the same inclination)
            rot: rotation of system plane about z axis (rad) {Float64}
            nr: number of radial bins {Int}
            nϕ: number of azimuthal bins {Int}
            I: intensity function {Function} (defaults to DiskWindIntensity)
            v: velocity function {Function} (defaults to vCircularDisk)
            scale: radial binning scale (:log or :linear)
            kwargs: extra keyword arguments for I and v if they are functions (see examples)
        returns:
            model object {model}
        """
        r = nothing; ΔA = nothing
        @assert rMin < rMax "rMin must be less than rMax"
        @assert nr > 0 "nr must be greater than 0"
        @assert nϕ > 0 "nϕ must be greater than 0"

        ϕ = collect(range(-π,stop=π,length=nϕ)) #non-rotated frame
        Δϕ = ϕ[2] - ϕ[1]
        if scale == :log        
            logr = collect(range(log(rMin*cos(i)),stop=log(rMax),length=nr))
            Δr = logr[2] - logr[1]
            r = exp.(logr)
        elseif scale == :linear
            r = collect(range(rMin*cos(i),stop=rMax,length=nr))
            Δr = r[2] - r[1]
        else
            error("scale must be :log or :linear")
        end

        rMesh, ϕMesh = meshgrid(r,ϕ) #camera r,ϕ
        α = rMesh .* cos.(ϕMesh); β = rMesh .* sin.(ϕMesh) #camera coordinates
        ΔA = scale == :log ? rMesh.^2 .* (Δr * Δϕ) : rMesh .* (Δr * Δϕ) #projected disk area, normalization doesn't matter
        rSystem = zeros(nr,nϕ); ϕSystem = zeros(nr,nϕ); ϕₒ = zeros(nr,nϕ)
        θₒ = 0.0; rot = 0.0
        r3D = get_r3D(i,rot,θₒ) 
        xyz = [0.0;0.0;0.0]
        matBuff = zeros(3,3)
        colBuff = zeros(3)
        rt = 0.0; ϕt = 0.0; ϕₒt = 0.0 #preallocate raytracing variables
        for ri in 1:nr
            for ϕi in 1:nϕ
                rt, ϕt, ϕₒt = raytrace(α[ri,ϕi], β[ri,ϕi], i, rot, θₒ, r3D, xyz, matBuff, colBuff) #flip i to match convention of +z being up, relic
                # println("RAYTRACE: rt = $rt, ϕt = $ϕt, ϕₒt = $ϕₒt")
                # x = β[ri,ϕi]/cos(i); y = α[ri,ϕi]; z = 0.0 #system coordinates from camera coordinates, raytraced back to disk plane
                # rt = sqrt(x^2 + y^2 + z^2); ϕt = atan(y,x); ϕₒt = atan(y,x) #convert to polar coordinates
                # println("OLD WAY: rt = $rt, ϕt = $ϕt, ϕₒt = $ϕₒt")
                # exit()
                if rt < rMin || rt > rMax #exclude portions outside of (rMin, rMax)
                    rSystem[ri,ϕi], ϕSystem[ri,ϕi], ϕₒ[ri,ϕi] = NaN, NaN, NaN
                else
                    rSystem[ri,ϕi], ϕSystem[ri,ϕi], ϕₒ[ri,ϕi] = rt, ϕt, ϕₒt
                end
            end
        end

        rSystem = [rSystem[i,:] for i in 1:nr]; ϕSystem = [ϕSystem[i,:] for i in 1:nr]; ΔA = [ΔA[i,:] for i in 1:nr]; ϕₒ = [ϕₒ[i,:] for i in 1:nr] #reshape, correct ϕ for other functions (based on ϕ to observer with ϕ = 0 at camera)
        rings = [ring(r = ri, i = i, v = v, I = I, ϕ = ϕi, ϕₒ = ϕₒi, ΔA = ΔAi, rMin=rMin, rMax=rMax, rot=rot, θₒ=θₒ; kwargs...) for (ri,ϕi,ΔAi,ϕₒi) in zip(rSystem,ϕSystem,ΔA,ϕₒ)]
        m = new(rings,Dict{Symbol,profile}(),camera(stack(α,dims=1),stack(β,dims=1)))
    end

    function model(r̄::Float64, rFac::Float64, Sα::Float64, i::Float64, nr::Int, nϕ::Int, scale::Symbol; kwargs...)
        """constructor for model struct -- takes in r̄, Sα, i, nr, nϕ, scale, and kwargs and returns a model object (detailed above) while checking for errors
        this version of the constructor creates a DiskWind model of the BLR as detailed in Long+2023
        params:
            r̄: mean radius of model (in terms of rₛ) {Float64}
            Sα: radius factor {Float64} 
            i: inclination angle (rad) {Float64}
            nr: number of radial bins {Int}
            nϕ: number of azimuthal bins {Int}
            scale: radial binning scale (:log or :linear)
            kwargs: extra keyword arguments for model constructor (see examples)
        returns:
            model object {model}
        """
        @assert (Sα != 1/2) && (Sα != 3/2) "Sα cannot be 1/2 or 3/2 as this divides by zero"
        @assert r̄ > 0 "r̄ must be greater than 0"

        rMin, rMax = get_rMinMaxDiskWind(r̄,rFac,Sα)
        # println("r̄ = $r̄, rMin = $rMin, rMax = $rMax")
        kwargs = values(kwargs); kwargs = merge(kwargs, (; α = Sα))
        # println("kwargs = $kwargs")
        model(rMin, rMax, i, nr, nϕ, DiskWindIntensity, vCircularDisk, scale; kwargs...)
    end
end

function DiskWindModel(rMin::Float64, rMax::Float64, i::Float64; nr::Int=128, nϕ::Int=256, I::Function=DiskWindIntensity, v::Function=vCircularDisk, scale::Symbol=:log, kwargs...)
    return model(rMin, rMax, i, nr, nϕ, I, v, scale; kwargs...)
end

function DiskWindModel(r̄::Float64, rFac::Float64, α::Float64, i::Float64; rot::Float64=0.0, nr::Int=128, nϕ::Int=256, scale::Symbol=:log, kwargs...)
    return model(r̄, rFac, α, i, nr, nϕ, scale; kwargs...)
end

function cloudModel(ϕₒ::Vector{Float64}, i::Vector{Float64}, rot::Vector{Float64}, θₒ::Vector{Float64}, θₒSystem::Float64, ξ::Float64; rₛ::Float64=1.0, μ::Float64=500., β::Float64=1.0, F::Float64=0.5,
    I::Union{Function,Float64}=IsotropicIntensity,v::Union{Function,Float64}=vCircularCloud,kwargs...)

    @assert length(ϕₒ) == length(i) == length(rot) == length(θₒ) "ϕ, i, rot, and θₒ must be the same length -- got $(length(ϕ)), $(length(i)), $(length(rot)), and $(length(θₒ))"
    rings = [drawCloud(i=i[j],θₒ=θₒ[j],rot=rot[j],ϕₒ=ϕₒ[j],μ=μ,F=F,β=β,rₛ=rₛ,θₒSystem=θₒSystem,I=I,v=v,ξ=ξ;kwargs...) for j=1:length(ϕₒ)]
    return model(rings)
end

function cloudModel(nClouds::Int64; μ::Float64=500., β::Float64=1.0, F::Float64=0.5, rₛ::Float64=1.0, θₒ::Float64=π/2, γ::Float64=1.0, ξ::Float64=1.0, i::Float64=0.0, 
    I::Union{Function,Float64}=IsotropicIntensity, v::Union{Function,Float64}=vCircularCloud, rng::AbstractRNG=Random.GLOBAL_RNG, kwargs...)
    ϕₒ = rand(rng,nClouds).*2π
    #θₒ = rand(nClouds).*θₒ #note: need to implement equation 13 here -- should add a system θₒ parameter to ring (or model?) struct, then each ring can have a different θ
    θ = acos.(cos(θₒ).+(1-cos(θₒ)).*rand(rng,nClouds).^γ) #θₒ for each cloud, from eqn 14
    #idea: models should have set i, θₒ, but define + operator to "add" two models together that may have rings with different i, θₒ
    rot = rand(rng,nClouds).*2π
    i = ones(nClouds).*i
    return cloudModel(ϕₒ,i,rot,θ,θₒ,ξ, rₛ=rₛ,μ=μ,β=β,F=F,I=I,v=v,rng=rng;kwargs...)
end

function cloudModel(nClouds::Int64,μ::Float64,β::Float64,F::Float64,θₒ::Float64,i::Float64,rₛ::Float64=1.0,I=IsotropicIntensity,v=vCircularCloud,rng::AbstractRNG=Random.GLOBAL_RNG)
    γ = getGamma(μ=μ,β=β,F=F)
    rings = Array{ring{Vector{Float64},Float64},1}(undef,nClouds)
    for n=1:nClouds
        rCam = getR(rₛ,γ,rng)
        ϕCam = rand(rng)*2π
        θₒn = rand(rng)*θₒ
        rot = rand(rng)*2π
        α,β = rCam*cos(ϕCam),rCam*sin(ϕCam)
        r,ϕ = raytrace(α,β,i,rot,θₒn)
        rings[n] = ring(r=r,i=i,rot=rot,θₒ=θₒn,v=v,I=I,ϕ=ϕ,ΔA=1.0)
    end
    return model(rings)
end



Base.show(io::IO, m::model) = begin 
    println(io, "model struct with $(length(m.rings)) rings:")
    if isdefined(m, :profiles) && length(m.profiles) > 0
        println(io, "\t-profiles: $(keys(m.profiles))")
    else
        println(io, "\t-no profiles set")
    end
    if isdefined(m, :camera)
        println(io, "\t-$(m.camera)")
    else
        println(io, "\t-no camera set")
    end
end
