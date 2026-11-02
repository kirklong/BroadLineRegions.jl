@testset "Camera-to-ring coordinate round trip" begin
    for (i, rot, θₒ) in ((π/6, 0.0, 0.0), (0.4, 0.3, 0.2), (0.7, -0.4, 0.1))
        r3D = BLR.get_r3D(i, rot, θₒ)
        M = [sin(i) 0.0 -cos(i); 0.0 1.0 0.0; cos(i) 0.0 sin(i)] * r3D
        for (α, β) in ((50.0, 25.0), (-30.0, 40.0), (20.0, -60.0))
            r, ϕ, ϕ₀ = BLR.raytrace(α, β, i, rot, θₒ)
            αCheck, βCheck = BLR.photograph(r, ϕ₀, i, rot, θₒ)
            @test αCheck ≈ α
            @test βCheck ≈ β
            expected = BLR.raytrace(α, β, i, rot, θₒ, M)
            @test r ≈ expected[1]
            @test ϕ ≈ expected[2]
            @test ϕ₀ ≈ expected[3]
        end
    end
end
