#!/usr/bin/env julia

using Test 

include("structs.jl")
include("velocity.jl")
include("intensity.jl")

@testset "e = 0.0, initialized successfully" begin
    @test (ring(r=10.0, i=π/4, e=0.0, v=-1.0, I=1.0, ϕ=π/1)).e == 0.0
    @test (ring(r=10.0, i=π/4, v=-1.0, I=1.0, ϕ=π/1)).e == 0.0
    @test (ring(r=10.0, i=π/4, e=0.0, v=-1.0, I=1.0, ϕ=π/1)).v == -1.0
    @test (ring(r=10.0, i=π/4, e=0.0, v=-1.0, I=1.0, ϕ=π/1)).I == 1.0
    @test (ring(r=10.0, i=π/4, e=0.0, v=-1.0, I=1.0, ϕ=π/1)).ϕ == π/1
    @test (ring(r=10.0, i=π/4, e=0.0, v=vCircular, I=1.0, ϕ=π/1)).v == vCircular(r=10.0, i=π/4, ϕ=π/1)
    @test (ring(r=10.0, i=π/4, e=0.0, v=-1.0, I=IsotropicIntensity, ϕ=π/1)).I == IsotropicIntensity(r=10.,ϕ=π/1)
    ϕ = collect(range(0, stop=2π, length=10))
    @test (ring(r=10.0, i=π/4, e=0.0, v=vCircular, I=IsotropicIntensity, ϕ=ϕ)).I == IsotropicIntensity(r=10.,ϕ=ϕ)
    @test (ring(r=10.0, i=π/4, e=0.0, v=vCircular, I=IsotropicIntensity, ϕ=ϕ)).v == vCircular(r=10.0, i=π/4, ϕ=ϕ)
end

@testset "e != 0.0, initialized successfully" begin
    @test (ring(r=rElliptical, i=π/4, e=0.5, a=5., v=-1.0, I=1.0, ϕ=π/1)).e == 0.5
    @test (ring(r=rElliptical, i=π/4, e=0.5, a=5., v=-1.0, I=1.0, ϕ=π/1)).r == rElliptical(ϕ=π/1, a=5., e=0.5)
    @test (ring(r=rElliptical, i=π/4, e=0.5, a=5., v=vElliptical, I=1.0, ϕ=π/1)).v == vElliptical(a=5., i=π/4, ϕ=π/1, e=0.5)
    ϕ = collect(range(0, stop=2π, length=10)); r = rElliptical(ϕ=ϕ, a=5., e=0.5)
    @test (ring(r=rElliptical, i=π/4, e=0.5, a=5., v=vElliptical, I=IsotropicIntensity, ϕ=ϕ)).v == vElliptical(a=5., i=π/4, ϕ=ϕ, e=0.5)
    @test (ring(r=rElliptical, i=π/4, e=0.5, a=5., v=vElliptical, I=IsotropicIntensity, ϕ=ϕ)).I == IsotropicIntensity(r=r,ϕ=ϕ)
    @test (ring(r=rElliptical, i=π/4, e=0.5, a=5., v=vElliptical, I=IsotropicIntensity, ϕ=π/1)).r == rElliptical(ϕ=π/1, a=5., e=0.5)
end

@testset "e = 0.0, errors thrown correctly" begin
    @test_throws AssertionError ring(r=-1, i=π/4, e=0.0, v=-1.0, I=1.0, ϕ=π/1)
    @test_throws AssertionError ring(r=10.0, i=-π, e=0.0, v=-1.0, I=1.0, ϕ=π/1)
    @test_throws AssertionError ring(r=10.0, i=π/4, e=0.0, v="test", I=1.0, ϕ=π/1)
    @test_throws AssertionError ring(r=10.0, i=π/4, e=0.0, v=-1.0, I="test", ϕ=π/1)
    @test_throws AssertionError ring(r=10.0, i=π/4, e=0.0, v=-1.0, I=1.0, ϕ="test")
    @test_throws AssertionError ring(r=10.0, i=π/4, e=0.0, v=-1.0, I=1.0, ϕ=-π/1)
    @test_throws AssertionError ring(r=10.0, i=π/4, e=0.0, v=-1.0, I=1.0, ϕ=3π)
    @test_throws AssertionError ring(r=10.0, i=π/1, e=0.0, v=-1.0, I=1.0, ϕ=π/1)
    ϕ = collect(range(0, stop=2π, length=10)); r = rElliptical(ϕ=ϕ, a=5., e=0.5)
    @test_throws AssertionError (ring(r=10.0, i=π/4, e=0.0, v=-1.0, I=IsotropicIntensity, ϕ=ϕ)).I == IsotropicIntensity(ϕ=ϕ)
end


