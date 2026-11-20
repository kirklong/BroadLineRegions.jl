function test_device_callables(backend, ::Type{T}) where {T}
    @testset "On-device callables ($T on $(typeof(backend)))" begin
        # A callable struct or closure can own its prescription without extending BLR.
        radius = draw -> draw(0) < T(0.5) ? T(150) + T(20) * draw(1) : T(900) + T(20) * draw(1)
        opening = r -> r < T(500) ? T(0.3) : T(0.8)
        clouds = BLR.residentCloudModel(1024, 42; i=0.4, radial_sampler=radius, θₒ=opening,
            γ=r -> one(r), intensity=p -> p.opening_angle, velocity=(p, draw) -> draw(0),
            responsivity=p -> p.r < T(500) ? one(T) : T(2), backend, T)
        repeat_clouds = BLR.residentCloudModel(1024, 42; i=0.4, radial_sampler=radius,
            velocity=(p, draw) -> draw(0), backend, T)
        r = Array(clouds.ma.r)
        inner = r .< T(500)
        @test 400 < count(inner) < 624
        @test all(T(150) .<= r[inner] .<= T(170))
        @test all(T(900) .<= r[.!inner] .<= T(920))
        @test all(Array(clouds.ma.θₒ) .<= opening.(r) .+ 8eps(T))
        @test Array(clouds.ma.I) == opening.(r)
        @test Array(clouds.ma.η) == ifelse.(inner, one(T), T(2))
        @test Array(repeat_clouds.ma.r) == r
        @test Array(repeat_clouds.ma.v) == Array(clouds.ma.v)
        @test all(zero(T) .< Array(clouds.ma.v) .< one(T))

        # Callables returning the built-in constants preserve the turbulent velocity.
        turbulent = (; i=0.4, velocity=:turbulent, fFlow=1.0, σρᵣ=0.1, σΘᵣ=0.1, backend, T)
        fixed = BLR.residentCloudModel(128, 42; turbulent..., θₑ=0.5, σₜ=0.125)
        callable = BLR.residentCloudModel(128, 42; turbulent..., θₑ=r -> T(0.5),
            σₜ=(r, opening_angle, gamma) -> T(0.125))
        @test Array(callable.ma.v) == Array(fixed.ma.v)
        @test_throws ArgumentError BLR.residentCloudModel(2, 1; radial_sampler=17, backend, T)

        disk_kw = (; nr=8, nϕ=12, intensity=(r, ϕ, inc) -> one(T), backend, T)
        standard = BLR.residentDiskWindModel(100.0, 300.0, 0.4; disk_kw...)
        responsive = BLR.residentDiskWindModel(100.0, 300.0, 0.4; disk_kw...,
            responsivity=(r, ϕ, inc) -> T(2))
        inside = isfinite.(Array(standard.ma.r))
        @test any(inside)
        @test all(==(T(2)), Array(responsive.ma.η)[inside])
        @test isequal(Array(responsive.ma.r), Array(standard.ma.r))
        @test isequal(Array(responsive.ma.I), Array(standard.ma.I))
        t_edges = range(zero(T), T(650); length=50)
        @test Array(BLR.getΨt(responsive, t_edges)) ≈
            2 .* Array(BLR.getΨt(standard, t_edges)) rtol=32eps(T)
    end
end

test_device_callables(BLR.KernelAbstractions.CPU(), Float32)
test_device_callables(BLR.KernelAbstractions.CPU(), Float64)
