using KernelAbstractions

@testset "raytrace bin-assign kernels" begin
    disk(; r1=300.0, r2=900.0, nr=8, nϕ=16, inc=0.4, τ=0.4) =
        BLR.DiskWindModel(r1, r2, inc, nr=nr, nϕ=nϕ, scale=:linear,
            I=BLR.IsotropicIntensity, v=BLR.vCircularDisk, τ=τ, reflect=false)

    clouds(n, seed; μ=600.0, inc=0.4, τ=0.1) =
        BLR.cloudModel(n, μ=μ, β=1.0, F=0.5, θₒ=0.4, i=inc, γ=1.0, ξ=1.0,
            I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=τ, rng=MersenneTwister(seed))

    models = (
        disk() + clouds(80, 101),
        clouds(80, 101) + disk(),
        disk(r1=250.0, r2=700.0, nr=6, nϕ=12) + disk(r1=500.0, r2=1000.0, nr=6, nϕ=12),
        disk() + clouds(50, 102) + clouds(50, 103, μ=850.0),
    )

    for m in models
        camStartInds = BLR.getFlattenedCameraIndices(m)
        points = BLR._rt_flatten_points(m, camStartInds)
        grids, pixels = BLR._rt_build_output(m, camStartInds)[1:2]
        gridArrays = BLR._rt_grid_arrays(grids)
        ma = BLR.flatten(m)

        keys = BLR._rt_bin_assign(ma, gridArrays; backend=KernelAbstractions.CPU())
        ref = BLR._rt_reference_pixel_keys(points, grids)

        @test keys == ref
        @test all(0 .<= keys .<= length(pixels))
        @test count(==(0), keys) == count(p -> p.discrete && BLR._rt_find_pixel(p, grids) == 0, points)
    end
end
