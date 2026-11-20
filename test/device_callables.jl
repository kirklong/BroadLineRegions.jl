function test_device_callables(backend, ::Type{T}) where {T}
    @testset "On-device cloud callables ($T on $(typeof(backend)))" begin
        sampler = u -> u < T(0.5) ? T(150) + T(40) * u : T(900) + T(40) * (u - T(0.5))
        m = BLR.gpuCloudModel(1024, 42; i=0.4, radial_sampler=sampler,
            intensity=p -> p.r / T(150), velocity=p -> p.reflect ? -p.rₛ : p.rₛ, backend, T)
        r = Array(m.ma.r)
        inner = r .< T(500)
        @test 400 < count(inner) < 624
        @test all(T(150) .<= r[inner] .<= T(170))
        @test all(T(900) .<= r[.!inner] .<= T(920))
        repeat_model = BLR.residentCloudModel(1024, 42; i=0.4, radial_sampler=sampler, backend, T)
        @test Array(repeat_model.ma.r) == r
        @test Array(m.ma.I) ≈ r ./ T(150) rtol=8eps(T)
        @test Array(m.ma.v) == ifelse.(Array(m.ma.reflect), -one(T), one(T))
        @test_throws ArgumentError BLR.residentCloudModel(2, 1; radial_sampler=17, backend, T)
        @test_throws ArgumentError BLR.residentCloudModel(2, 1; intensity=17, backend, T)
    end
end

test_device_callables(BLR.KernelAbstractions.CPU(), Float32)
test_device_callables(BLR.KernelAbstractions.CPU(), Float64)
