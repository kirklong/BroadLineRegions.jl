#!/usr/bin/env julia

struct ring{V,F}
    """structure to hold parameters of model ring
    attributes:
        r: distance from central mass (in terms of rₛ) Union{Vector{Float64},Float64
            - in the case of circular orbits (e=0) or for a single angle r is a Float64, 
              otherwise it is a Vector{Float64} of distances corresponding to azimuthal angles ϕ
            - user can (optionally) supply a function to calculate r from other parameters (see constructor) that returns a Vector{Float64} or Float64
        i: inclination angle (rad) {Float64}
            - must be between 0 and π/2, with 0 being face-on and π/2 being edge-on
        e: eccentricity {Float64}
            - must be between 0 and 1, defaults to 0.0, if non-zero, r is a function of ϕ (must be supplied to constructor)
        v: line of sight velocity {Union{Vector{Float64},Float64}}
            - Float64 for a single angle, otherwise a Vector{Float64} of velocities corresponding to azimuthal angles ϕ
            - user can (optionally) supply a function to calculate v from other parameters (see constructor) that returns a Vector{Float64}
        I: intensity {Union{Vector{Float64},Float64,Matrix{Float64}}
            - Float64 for a single angle, otherwise a Vector{Float64} of intensities corresponding to azimuthal angles ϕ
            - user can (optionally) supply a function to calculate I from other parameters (see constructor) that returns a Vector{Float64}
        ϕ: azimuthal angle (rad) {Union{Vector{Float64},Float64}}
            - Float64 for a single angle, otherwise a Vector{Float64} of azimuthal angles between 0 and 2π
    
    constructor:
        ring(;r, i, e, v, I, ϕ, kwargs...)
        r: Float64, Vector{Float64} or Function
            - if function, must return a Vector{Float64} or Float64 corresponding to ϕ
        i: Float64
            - must be between 0 and π/2
        e: Float64
            - optional
            - must be between 0 and 1, defaults to 0.0
            - if non-zero, r should be a function (of at least ϕ) -- but can still be passed (with warning) as a Float64/Vector{Float64} if calculated elsewhere
        v: Float64, Vector{Float64} or Function
            - if function, must return a Vector{Float64} corresponding to ϕ
        I: Float64, Vector{Float64}, Matrix{Float64}, or Function
            - if function, must return a Vector{Float64} corresponding to ϕ for each r (if r is a Vector{Float64} must return a Matrix{Float64} of intensities for each r)
        ϕ: Float64 or Vector{Float64}
            - must be between 0 and 2π
            - if Float64, r, v, and I must be Float64 (or functions that return Float64)
        ΔA: Float64 or Vector{Float64}
            - projected area of ring in image (used in calculating profiles)
            - i.e. for log-spherical ring ΔA = r^2*Δr*Δϕ, for linear-spherical ring ΔA = r*Δr*Δϕ, for cloud ΔA could be size of cloud (note here r would be the image r not the physical r)
            - defaults to 1.0 if not provided
        kwargs: contain extra keyword arguments for v, I, and/or r if they are functions (see examples)
    """

    r::Union{V,F,Function}
    i::F
    e::F
    v::Union{V,F,Function}
    I::Union{V,F,Matrix{F},Function}
    ϕ::Union{V,F}
    ΔA::Union{V,F}

    function ring(;kwargs...) #could re-write this to use multiple dispatch? i.e. ring(;r::Float64, i::Float64, e::Float64, v::Float64, I::Float64, ϕ::Float64) etc.
        """
        constructor for ring struct -- takes in kwargs (detailed above) and returns a ring object (detailed above) while checking for errors
        """
        r = nothing; i = nothing; v = nothing; I = nothing; ϕ = nothing; e = 0.0; ΔA = 1.0
        try; r = kwargs[:r]; catch; error("r must be provided as kwarg"); end
        try; i = kwargs[:i]; catch; error("i must be provided as kwarg"); end
        try; e = kwargs[:e]; catch; println("e not provided: defaulting to 0.0"); end
        try; v = kwargs[:v]; catch; error("v must be provided as kwarg"); end
        try; I = kwargs[:I]; catch; error("I must be provided as kwarg"); end
        try; ϕ = kwargs[:ϕ]; catch; error("ϕ must be provided as kwarg"); end
        try; ΔA = kwargs[:ΔA]; catch; error("ΔA not provided, defaulting to 1.0"); end
        kwargs = values(kwargs)

        @assert typeof(i) == Float64 "i must be Float64, got $(typeof(i))"
        @assert (i >= 0) && (i <= π/2) "i must be between 0 and π/2, got $i"
        @assert typeof(e) == Float64 "e must be Float64, got $(typeof(e))"
        @assert (e >= 0) && (e < 1) "e must be between 0 and 1, got $e"
        @assert (typeof(ΔA) == Float64) || (typeof(ΔA) == Vector{Float64}) "ΔA must be Float64 or Vector{Float64}, got $(typeof(ΔA))"
        if typeof(ΔA) == Vector{Float64}
            @assert length(ΔA) == length(ϕ) "ΔA must be the same length as ϕ"
            @assert all(ΔA .> 0) "ΔA must be greater than 0"
        else
            @assert ΔA > 0 "ΔA must be greater than 0"
        end
        @assert (typeof(r) == Vector{Float64}) || (typeof(r) == Float64) || (isa(r,Function)) "r must be Float64, Vector{Float64} or Function, got $(typeof(v))"        
        if (e != 0.0)
            if isa(r,Function)
                try; r = r(;kwargs...); kwargs = merge(kwargs, (; r = r)); catch; error("error in function r -- check kwargs"); end
                @assert length(r) == length(ϕ) "r must return a vector of the same length as ϕ"
                @assert typeof(r) == typeof(ϕ) "r and ϕ must be the same type"
                @assert all(r .>= 0) "r must be greater than 0"
            else
                @warn "e is non-zero but r is not a function -- make sure r is calculated elsewhere correctly"
                if typeof(r) == Float64
                    @assert r > 0 "r must be greater than 0"
                else
                    @assert all(r .> 0) "r must be greater than 0"
                end
            end
        end

        @assert (typeof(ϕ) == Vector{Float64}) || (typeof(ϕ) == Float64) "ϕ must be Float64 or Vector{Float64}, got $(typeof(ϕ))"

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

        new{Vector{Float64},Float64}(r,i,e,v,I,ϕ,ΔA)
    end
end

Base.show(io::IO, r::ring) = begin
    println(io, "ring struct with inclination $(r.i) rad")
    println(io, "eccentricity: $(r.e)")
    if typeof(r.ϕ) == Float64
        println(io, "cloud with azimuthal angle: $(r.ϕ) rad")
    else
        println(io, "intensity distributed over $(length(r.ϕ)) azimuthal angles ($(minimum(r.ϕ)) < ϕ < $(maximum(r.ϕ)) rad)")
    end
    if typeof(r.r) == Float64
        println(io, "cloud at radius: $(r.r) rₛ")
    else
        println(io, "intensity distributed over $(length(r.r)) radii ($(minimum(r.r)) < r < $(maximum(r.r)) rₛ)")
    end
    if typeof(r.v) == Float64
        println(io, "cloud line of sight velocity: $(r.v) c")
    else
        println(io, "line of sight velocity: $(minimum(r.v)) < v < $(maximum(r.v)) c")
    end
    if typeof(r.I) == Float64
        println(io, "cloud intensity: $(r.I) arbitrary units")
    else
        println(io, "intensity: $(minimum(r.I)) < I < $(maximum(r.I)) arbitrary units")
    end
    if typeof(r.ΔA) == Float64
        println(io, "projected area: $(r.ΔA) rₛ²")
    else
        println(io, "projected area: $(minimum(r.ΔA)) < ΔA < $(maximum(r.ΔA)) rₛ²")
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

struct camera
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
    xRes, yRes = nothing, nothing
    if typeof(c.α) == Matrix{Float64}
        xRes, yRes = size(c.α)
    else
        xRes, yRes = length(c.α), length(c.β)
    end
    println(io, "camera with resolution ($xRes, $yRes) and range: $(minimum(c.α)) < α < $(maximum(c.α)) and $(minimum(c.β)) < β < $(maximum(c.β))")
end 

meshgrid(x,y) = (reshape(repeat(x,outer=length(y)),length(x),length(y)), reshape(repeat(y,inner=length(x)),length(x),length(y)))

mutable struct model #restructure to make mutable struct -- Δr doesn't need to be kept track of it only matters for binning profiles? profiles should be a dict i.e. LP => array, delayProfile => array, custom => array
    #make immutable to see if improves performance cost
    #add image as keyword arg to constructors -- only initialize if true to save on performance
    
    """structure to hold many rings and their parameters
    attributes:
        rings: Vector{ring}
            - list of ring objects

    constructor:
        ***TO DO***: want several different constructor methods based on whether user provides rings (and Δr, Δϕ) or not
        "ray tracing" should be in separate raytrace.jl file with helper methods for making these profiles from ensembles of rings
        should struct be mutable? 
    """
    rings::Vector{ring}
    profiles::Union{Nothing,Dict{Symbol,profile}}
    camera::Union{Nothing,camera}

    function model(rings::Vector{ring{Vector{Float64},Float64}})
        """
        constructor for model struct -- takes in rings and returns a model object (detailed above) while checking for errors
        """
        r = [r.r for r in rings]; ϕ = [r.ϕ for r in rings]; i = [r.i for r in rings]
        r_c = zeros(length(r)); ϕ_c = zeros(length(r))
        for (i,(ri,ϕi,ii)) in enumerate(zip(r,ϕ,i))
            r_c[i], ϕ_c[i] = photograph(ri,ϕi,ii) #get camera coordinates from physical 
        end
        α = r_c.*cos.(ϕ_c); β = r_c.*sin.(ϕ_c) #camera x and y
        new(rings,Dict{Symbol,profile}(),camera(stack(α,dims=1),stack(β,dims=1)))
    end

    function model(rMin::Float64, rMax::Float64, i::Float64, nr::Int, nϕ::Int, I::Function, v::Function, scale::Symbol; kwargs...)
        """constructor for model struct -- takes in rMin, rMax, i, nr, nϕ, I, v, and scale and returns a model object (detailed above) while checking for errors
        params:
            rMin: minimum radius of model (in terms of rₛ) {Float64}
            rMax: maximum radius of model (in terms of rₛ) {Float64}
            i: inclination angle (rad) {Float64} (all rings have the same inclination)
            nr: number of radial bins {Int}
            nϕ: number of azimuthal bins {Int}
            I: intensity function {Function} (defaults to DiskWindIntensity)
            v: velocity function {Function} (defaults to vCircular)
            scale: radial binning scale (:log or :linear)
            kwargs: extra keyword arguments for I and v if they are functions (see examples)
        returns:
            model object {model}
        """
        r = nothing; ΔA = nothing
        @assert rMin < rMax "rMin must be less than rMax"
        @assert nr > 0 "nr must be greater than 0"
        @assert nϕ > 0 "nϕ must be greater than 0"

        ϕ = collect(range(0,stop=2π,length=nϕ))
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

        rMesh, ϕMesh = meshgrid(r,ϕ)
        # ΔA = scale == :log ? @. rMesh^2*Δr*Δϕ : @. rMesh*Δr*Δϕ
        ΔA = scale == :log ? [rMesh[i,:].^2 .* (Δr * Δϕ) for i in 1:nr] : [rMesh[i,:] .* (Δr * Δϕ) for i in 1:nr]
        # α = r*cos(ϕ); β = @. r*sin(ϕ) #camera coordinates in x and y
        α = [ri.*cos.(ϕ) for ri in r]; β = [ri.*sin.(ϕ) for ri in r] #camera coordinates in x and y
        r = [sqrt.(αi.^2 .+ (βi./cos(i)).^2) for (αi,βi) in zip(α,β)] #converting camera r to physical r, raytraced back to disk
        ϕ = [atan.(βi./cos(i),αi) for (αi,βi) in zip(α,β)] #converting camera ϕ to physical ϕ, raytraced from camera to disk plane
        # r = @. sqrt(α^2 + (β/$cos(i))^2) #converting camera r to physical r, raytraced back to disk
        # ϕ = @. atan(β/$cos(i),α) #converting camera ϕ to physical ϕ, raytraced from camera to disk plane
        rings = [ring(r = ri, i = i, e = 0.0, v = v, I = I, ϕ = ϕi, ΔA = ΔAi, rMin=rMin, rMax=rMax; kwargs...) for (ri,ΔAi,ϕi) in zip(r,ΔA,ϕ)]
        """PROBLEM: need to recalculate r and ϕ for raytracing reasons (where ray hits disk)"""
        new(rings,Dict{Symbol,profile}(),camera(stack(α,dims=1),stack(β,dims=1)))
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
        model(rMin, rMax, i, nr, nϕ, DiskWindIntensity, vCircular, scale; kwargs...)
    end
end

function model(rMin::Float64, rMax::Float64, i::Float64; nr::Int=128, nϕ::Int=256, I::Function=DiskWindIntensity, v::Function=vCircular, scale::Symbol=:log, kwargs...)
    return model(rMin, rMax, i, nr, nϕ, I, v, scale; kwargs...)
end

function model(r̄::Float64, rFac::Float64, α::Float64, i::Float64; nr::Int=128, nϕ::Int=256, scale::Symbol=:log, kwargs...)
    return model(r̄, rFac, α, i, nr, nϕ, scale; kwargs...)
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
