using Adapt

@testset "ModelArrays flatten/cpu" begin
    mD = BLR.DiskWindModel(300.0, 900.0, 0.4, nr=16, nϕ=32, scale=:linear,
        I=BLR.IsotropicIntensity, v=BLR.vCircularDisk, τ=0.4, reflect=false)
    mC = BLR.cloudModel(250, μ=600.0, β=1.0, F=0.5, θₒ=0.4, i=0.4, γ=1.0, ξ=0.8,
        I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=0.1, rng=MersenneTwister(222))
    models = (mD, mC, mD + mC)

    for m in models
        ma = BLR.flatten(m)
        n = length(BLR.getVariable(m, :I, flatten=true))
        @test ma isa BLR.ModelArrays{Float64,Vector{Float64},Vector{Bool}}
        @test all(length(getfield(ma, f)) == n for f in fieldnames(typeof(ma)))
        @test ma.α == Float64.(vec(m.camera.α))
        @test ma.β == Float64.(vec(m.camera.β))
        @test isequal(ma.I, BLR.getVariable(m, :I, flatten=true))
        @test isequal(ma.ΔA, BLR.getVariable(m, :ΔA, flatten=true))
        @test isequal(ma.r, BLR.getVariable(m, :r, flatten=true))
        @test isequal(ma.x, BLR.getVariable(m, :x, flatten=true))
        @test isequal(ma.y, BLR.getVariable(m, :y, flatten=true))
        @test isequal(ma.z, BLR.getVariable(m, :z, flatten=true))
        @test isapprox(sum(ma.I .* ma.ΔA), sum(BLR.getVariable(m, :I, flatten=true) .* BLR.getVariable(m, :ΔA, flatten=true)), rtol=1e-12)
        @test eltype(ma.reflect) == Bool
        @test BLR.cpu(ma) == ma
        @test Adapt.adapt(Array, ma) == ma
    end

    ma32 = BLR.flatten(mD; T=Float32)
    @test ma32 isa BLR.ModelArrays{Float32,Vector{Float32},Vector{Bool}}
    @test eltype(ma32.I) == Float32
    @test ma32.α == Float32.(vec(mD.camera.α))

    @test BLR.cpu(mC) isa BLR.ModelArrays{Float64,Vector{Float64},Vector{Bool}}
    @test_throws ErrorException BLR.gpu(mD)
end
