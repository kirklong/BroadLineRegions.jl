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

        intensity = (r, ϕ, inc) -> inv(r)
        velocity = (r, ϕ, inc) -> ϕ * inc / r
        cloud = BLR.gpuCloudModel(128, 42; i=0.4, intensity, velocity, backend, T)
        disk = BLR.gpuDiskWindModel(100., 300., 0.4; nr=8, nϕ=16, intensity, velocity, backend, T)
        for model in (cloud, disk)
            r, ϕ, inc = Array(model.ma.r), Array(model.ma.ϕ), Array(model.ma.i)
            valid = isfinite.(r)
            @test any(valid)
            @test Array(model.ma.I)[valid] ≈ intensity.(r[valid], ϕ[valid], inc[valid]) rtol=8eps(T)
            @test Array(model.ma.v)[valid] ≈ velocity.(r[valid], ϕ[valid], inc[valid]) rtol=8eps(T)
        end

        # One adapter can request any number of fields in a chosen order.
        emissivity = (r, ϕ, inc, scale_radius) -> scale_radius / r
        flow = (r, ϕ, inc, x, opening, z) -> (abs(x) + opening * abs(z)) / r
        extra_kw = (; θₒ=0.2, rₛ=3.0, intensity=BLR.deviceCallback(emissivity, :rₛ),
            velocity=BLR.deviceCallback(flow, :x, :θₒSystem, :z), backend, T)
        cloud = BLR.gpuCloudModel(128, 42; i=0.4, extra_kw...)
        disk = BLR.gpuDiskWindModel(100., 300., 0.4; nr=8, nϕ=16, extra_kw...)
        for model in (cloud, disk)
            r, x, z = Array(model.ma.r), Array(model.ma.x), Array(model.ma.z)
            valid = isfinite.(r)
            @test any(valid)
            @test Array(model.ma.I)[valid] ≈ T(3) ./ r[valid] rtol=8eps(T)
            @test Array(model.ma.v)[valid] ≈ (abs.(x[valid]) .+ T(0.2) .* abs.(z[valid])) ./ r[valid] rtol=8eps(T)
        end
        @test_throws ArgumentError BLR.residentCloudModel(2, 1;
            intensity=BLR.deviceCallback(emissivity, :unknown_field), backend, T)
        @test_throws ArgumentError BLR.residentCloudModel(2, 1;
            intensity=BLR.deviceCallback(emissivity, :x, :z), backend, T)
    end
end

test_device_callables(BLR.KernelAbstractions.CPU(), Float32)
test_device_callables(BLR.KernelAbstractions.CPU(), Float64)
