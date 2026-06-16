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

        perm = BLR._rt_sortperm_by_key_depth(ma, keys)
        refPerm = sortperm(1:length(points); by=i -> (ref[i], -BLR._rt_sort_depth(points[i].x)), alg=MergeSort)
        @test perm == refPerm
        @test BLR._rt_sorted_key_depth_pairs(keys, ma.x, perm) ==
              BLR._rt_sorted_key_depth_pairs(ref, [p.x for p in points], refPerm)

        IR = [1.0 + 0.1*(s-1) for s in 1:length(m.subModelStartInds)]
        scan = BLR._rt_scan_arrays(m, keys, perm, pixels, IR)
        out = BLR._rt_scan_output(length(pixels))
        BLR._rt_segmented_scan!(out, ma, keys, perm, scan; τCutOff=1.0, backend=KernelAbstractions.CPU())

        for pixInd in eachindex(pixels)
            bucket = findall(==(pixInd), ref)
            res = BLR._rt_scan_bucket(points, bucket, IR, pixels[pixInd].ΔA, 1.0)
            @test out.active[pixInd] == !isnothing(res)
            if !isnothing(res)
                @test isapprox(out.I[pixInd], res.I, rtol=1e-10)
                @test isapprox(out.v[pixInd], res.v, rtol=1e-10)
                @test isapprox(out.r[pixInd], res.r, rtol=1e-10)
                @test isapprox(out.ϕ[pixInd], res.ϕ, rtol=1e-10)
                @test isapprox(out.ϕ₀[pixInd], res.ϕ₀, rtol=1e-10)
                @test isapprox(out.i[pixInd], res.i, rtol=1e-10)
                @test isapprox(out.rot[pixInd], res.rot, rtol=1e-10)
                @test isapprox(out.θₒ[pixInd], res.θₒ, rtol=1e-10)
                @test isapprox(out.τ[pixInd], res.τ, rtol=1e-10)
                @test isapprox(out.η[pixInd], res.η, rtol=1e-10)
                @test isapprox(out.x[pixInd], res.x, rtol=1e-10)
                @test isapprox(out.y[pixInd], res.y, rtol=1e-10)
                @test isapprox(out.z[pixInd], res.z, rtol=1e-10)
                @test out.reflect[pixInd] == res.reflect
            else
                @test isnan(out.I[pixInd])
            end
        end
    end

    m32 = models[1]
    camStartInds32 = BLR.getFlattenedCameraIndices(m32)
    points32 = BLR._rt_flatten_points(m32, camStartInds32)
    grids32, pixels32 = BLR._rt_build_output(m32, camStartInds32)[1:2]
    ma32 = BLR.flatten(m32; T=Float32)
    keys32 = BLR._rt_bin_assign(ma32, BLR._rt_grid_arrays(grids32; T=Float32); backend=KernelAbstractions.CPU())
    perm32 = BLR._rt_sortperm_by_key_depth(ma32, keys32)
    IR32 = Float32[1.0, 0.25]
    scan32 = BLR._rt_scan_arrays(m32, keys32, perm32, pixels32, IR32; T=Float32)
    out32 = BLR._rt_scan_output(length(pixels32); T=Float32)
    BLR._rt_segmented_scan!(out32, ma32, keys32, perm32, scan32; τCutOff=Float32(1.0), backend=KernelAbstractions.CPU())
    ref32 = BLR._rt_reference_pixel_keys(points32, grids32)
    for pixInd in eachindex(pixels32)
        res = BLR._rt_scan_bucket(points32, findall(==(pixInd), ref32), Float64.(IR32), pixels32[pixInd].ΔA, 1.0)
        if !isnothing(res)
            @test isapprox(out32.I[pixInd], res.I, rtol=1e-4)
            @test isapprox(out32.v[pixInd], res.v, rtol=1e-4)
            @test isapprox(out32.r[pixInd], res.r, rtol=1e-4)
        end
    end

    keys = [2, 1, 2, 1, 0, 2, 0]
    x = [10.0, 4.0, 11.0, 6.0, 3.0, NaN, 9.0]
    perm = BLR._rt_sortperm_by_key_depth(keys, x)
    @test keys[perm] == [0, 0, 1, 1, 2, 2, 2]
    @test x[perm[1:2]] == [9.0, 3.0]
    @test x[perm[3:4]] == [6.0, 4.0]
    @test isequal(x[perm[7]], NaN)
end
