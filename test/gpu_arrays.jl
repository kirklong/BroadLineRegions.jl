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

@testset "GPU backend selection + Float32-clean scalar kernels" begin
    if Base.get_extension(BLR, :BroadLineRegionsCUDAExt) === nothing &&
            Base.get_extension(BLR, :BroadLineRegionsMetalExt) === nothing
        err = try
            BLR.defaultGPUBackend()
        catch e
            e
        end
        @test err isa ErrorException && occursin("using CUDA", err.msg) && occursin("using Metal", err.msg)
        @test_throws ErrorException BLR.gpuCloudModel(10, 1)
        @test_throws ErrorException BLR.gpuDiskWindModel(300.0, 900.0, 0.4; f1=1.0, f2=1.0, f3=1.0, f4=1.0, α=1.0)
    end
    @test BLR._gpu_array_type(BLR.KernelAbstractions.CPU()) === Array
    @test BLR._check_gpu_eltype(BLR.KernelAbstractions.CPU(), Float64) === nothing
    # Metal has no Float64: the device scalar bodies must not promote Float32 inputs through untyped
    # Float64 literals (e.g. a bare sqrt(2) or π). Inferred return types expose any such promotion.
    noF64(rt) = !(Float64 in (rt isa DataType && rt <: Tuple ? fieldtypes(rt) : (rt,)))
    F = Float32
    @test only(Base.return_types(BLR._rt_disk_wind_i_scalar, NTuple{8,F})) === F
    @test only(Base.return_types(BLR._rt_v_circular_disk_scalar, NTuple{4,F})) === F
    @test noF64(only(Base.return_types(BLR._rt_disk_deproject_scalar, NTuple{21,F})))
    @test noF64(only(Base.return_types(BLR._rt_v_circular_radial_disk_scalar, NTuple{6,F})))
    @test only(Base.return_types(BLR._cloud_uniform, (Type{F}, UInt32, UInt32, UInt32, Int))) === F
    @test noF64(only(Base.return_types(BLR._rt_build_cloud_scalar,
        (Type{F}, UInt32, UInt32, UInt32, ntuple(_ -> F, 13)..., Bool, F, Bool, ntuple(_ -> F, 8)...))))
end
