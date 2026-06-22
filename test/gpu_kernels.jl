using KernelAbstractions

function _test_same_nan_approx(actual, expected; atol=1e-12, rtol=0.0)
    @test size(actual) == size(expected)
    for idx in eachindex(actual, expected)
        if isnan(expected[idx])
            @test isnan(actual[idx])
        else
            @test isapprox(actual[idx], expected[idx]; atol=atol, rtol=rtol)
        end
    end
end

function _disk_field_values(m, sym)
    out = Matrix{Float64}(undef, length(m.rings), length(m.rings[1].r))
    for (idx, ring) in enumerate(m.rings)
        out[idx, :] .= getfield(ring, sym)
    end
    return vec(out)
end

@testset "profile histogram kernels" begin
    edges = collect(range(-2.0, stop=3.0, length=38))
    rng = MersenneTwister(99)
    x = vcat(-2.0 .+ 5.0 .* rand(rng, 10_000), copy(edges), [NaN, Inf, -Inf, -3.0, 4.0])
    y = vcat(rand(rng, 10_000), ones(length(edges)), ones(5))
    for overflow in (false, true)
        out = zeros(length(edges)-1)
        BLR._rt_weighted_histogram!(out, x, y, edges; overflow=overflow, backend=KernelAbstractions.CPU())
        ref = BLR.binnedSum(x, y, bins=edges, overflow=overflow)[3]
        @test all(isapprox.(out, ref; rtol=1e-12, atol=1e-12))
    end

    disk = BLR.DiskWindModel(300.0, 900.0, 0.4, nr=16, nϕ=32, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.7, f3=0.2, f4=0.9,
        α=1.2, reflect=false)
    clouds = BLR.cloudModel(300, μ=600.0, β=1.0, F=0.5, θₒ=0.4, i=0.4, γ=1.0, ξ=0.8,
        I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=0.0, rng=MersenneTwister(12))
    m = disk + clouds
    lineEdges = collect(range(-0.08, stop=0.08, length=41))
    ma = BLR.flatten(m)
    lp = zeros(length(lineEdges)-1)
    BLR._rt_line_profile!(lp, ma, lineEdges; overflow=true, backend=KernelAbstractions.CPU())
    refLP = BLR.getProfile(m, :line, bins=lineEdges, centered=false, overflow=true).binSums
    @test all(isapprox.(lp, refLP; rtol=1e-12, atol=1e-12))

    θ = ma.α .* 1e-10
    w = ma.I .* ma.ΔA
    σ² = zeros(length(lineEdges)-1)
    sumW = similar(σ²)
    sumWθ = similar(σ²)
    μ = similar(σ²)
    sumWδ² = similar(σ²)
    BLR._rt_weighted_variance!(σ², sumW, sumWθ, μ, sumWδ², ma.v, θ, w, lineEdges;
        overflow=true, backend=KernelAbstractions.CPU())
    refσ² = BLR.binnedVariance(ma.v, θ, w, bins=lineEdges, overflow=true)[3]
    for idx in eachindex(σ², refσ²)
        if isnan(refσ²[idx])
            @test isnan(σ²[idx])
        else
            @test isapprox(σ²[idx], refσ²[idx]; rtol=1e-10, atol=1e-30)
        end
    end
end

@testset "transfer function kernels" begin
    disk = BLR.DiskWindModel(300.0, 900.0, 0.4, nr=12, nϕ=24, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.7, f3=0.2, f4=0.9,
        α=1.2, τ=0.4, reflect=false)
    clouds = BLR.cloudModel(200, μ=600.0, β=1.0, F=0.5, θₒ=0.4, i=0.4, γ=1.0, ξ=0.8,
        I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=0.1, rng=MersenneTwister(31))
    m = BLR.raytrace!(disk + clouds)   # realistic usage: transfer functions run on a raytraced model
    ma = BLR.flatten(m)

    # delay kernel τ = η(r − x); equals getVariable(m, BLR.t) to rounding (the unified W1 convention)
    delays = BLR._rt_transfer_delays(ma; backend=KernelAbstractions.CPU())
    @test delays ≈ ma.η .* (ma.r .- ma.x)
    @test all(isapprox.(delays, BLR.getVariable(m, BLR.t, flatten=true); rtol=1e-9, atol=1e-9))

    finV = filter(isfinite, ma.v); finT = filter(isfinite, delays)
    vEdges = collect(range(minimum(finV), maximum(finV), length=21))
    tEdges = collect(range(minimum(finT), maximum(finT), length=16))

    # 2D Ψ(v,t): kernel vs the CPU accumulator on identical inputs (binning correctness, pre-floor)
    Ψ = zeros(length(vEdges)-1, length(tEdges)-1)
    BLR._rt_psi2d!(Ψ, ma.v, delays, ma.I, ma.ΔA, vEdges, tEdges; backend=KernelAbstractions.CPU())
    Ψref = zeros(size(Ψ))
    BLR.ΨAccumulate!(Ψref, ma.v, delays, ma.I, ma.ΔA, vEdges, tEdges)
    @test Ψ ≈ Ψref

    # 1D Ψ(t) + overflow buckets vs ΨtAccumulate!
    Ψt = zeros(length(tEdges)-1); under = zeros(1); over = zeros(1)
    BLR._rt_psit!(Ψt, under, over, delays, ma.I, ma.ΔA, tEdges; backend=KernelAbstractions.CPU())
    Ψtref = zeros(length(tEdges)-1)
    sU, sO = BLR.ΨtAccumulate!(Ψtref, delays, ma.I, ma.ΔA, tEdges)
    @test Ψt ≈ Ψtref
    @test under[1] ≈ sU
    @test over[1] ≈ sO

    # unified convention: the transfer functions now bin EXACTLY like binnedSum, including values
    # sitting on edges (left-exclusive interior; boundary points fold into overflow). Construct
    # delays on every edge + interior + out-of-range and compare the kernel AND ΨtAccumulate! to
    # binnedSum directly.
    tE = collect(0.0:1.0:6.0)
    dExact = vcat(copy(tE), [0.5, 5.5, -1.0, 7.0])
    wExact = collect(1.0:length(dExact)); Ae = ones(length(dExact))
    refNoOvf = BLR.binnedSum(dExact, wExact .* Ae; bins=tE, overflow=false)[3]
    refOvf = BLR.binnedSum(dExact, wExact .* Ae; bins=tE, overflow=true)[3]

    ΨtK = zeros(length(tE)-1); uK = zeros(1); oK = zeros(1)
    BLR._rt_psit!(ΨtK, uK, oK, dExact, wExact, Ae, tE; backend=KernelAbstractions.CPU())
    @test ΨtK ≈ refNoOvf                                   # interior bins (kernel)
    assembled = copy(ΨtK); assembled[1] += uK[1]; assembled[end] += oK[1]
    @test assembled ≈ refOvf                               # overflow folding (kernel)

    ΨtA = zeros(length(tE)-1)
    sU, sO = BLR.ΨtAccumulate!(ΨtA, dExact, wExact, Ae, tE)
    @test ΨtA ≈ refNoOvf                                   # interior bins (CPU transfer fn)
    asmA = copy(ΨtA); asmA[1] += sU; asmA[end] += sO
    @test asmA ≈ refOvf                                    # overflow folding (CPU transfer fn)
end

@testset "resident-model observables (CPU backend)" begin
    nanapprox(a, b; rtol=1e-8, atol=1e-10) = length(a) == length(b) && all(
        (isnan(b[i]) ? isnan(a[i]) : isapprox(a[i], b[i]; rtol=rtol, atol=atol)) for i in eachindex(a, b))

    disk = BLR.DiskWindModel(300.0, 900.0, 0.4, nr=12, nϕ=24, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.7, f3=0.2, f4=0.9,
        α=1.2, ηₒ=0.4, η₁=0.6, αRM=0.1, rNorm=700.0, τ=0.4, reflect=false)
    clouds = BLR.cloudModel(200, μ=600.0, β=1.0, F=0.5, θₒ=0.4, i=0.4, γ=1.0, ξ=0.8,
        I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=0.1, rng=MersenneTwister(41))
    m = BLR.raytrace!(disk + clouds)
    rm = BLR.resident(m; backend=KernelAbstractions.CPU())
    @test rm isa BLR.ResidentModel

    U = [40.0, -12.0]; V = [5.0, 33.0]; PA = 0.6; BLRAng = 1e-11
    for bins in (60, collect(range(-0.06, 0.06, length=41)))
        for sym in (:line, :r, :ϕ)
            ref = BLR.getProfile(m, sym; bins=bins)
            got = BLR.getProfile(rm, sym; bins=bins)
            @test got.binEdges ≈ ref.binEdges
            @test nanapprox(got.binSums, ref.binSums; rtol=1e-10)
        end
        # :delay uses the unified η(r−x) delay (matches getΨt); ≈ host getProfile(:delay) to ~1e-8
        @test nanapprox(BLR.getProfile(rm, :delay; bins=bins).binSums, BLR.getProfile(m, :delay; bins=bins).binSums)
        for sym in (:phase, :moment2)
            ref = BLR.getProfile(m, sym; bins=bins, U=U, V=V, PA=PA, BLRAng=BLRAng)
            got = BLR.getProfile(rm, sym; bins=bins, U=U, V=V, PA=PA, BLRAng=BLRAng)
            @test nanapprox(got.binSums, ref.binSums; rtol=1e-8)
        end
    end

    fv = filter(isfinite, BLR.getVariable(m, :v, flatten=true))
    fd = filter(isfinite, BLR.getVariable(m, BLR.t, flatten=true))
    vEdges = collect(range(minimum(fv), maximum(fv), length=16))
    tEdges = collect(range(minimum(fd), maximum(fd), length=21))
    for overflow in (false, true)
        @test nanapprox(BLR.getΨt(rm, tEdges, overflow), BLR.getΨt(m, tEdges, overflow); rtol=1e-8)
    end
    @test nanapprox(vec(BLR.getΨ(rm, vEdges, tEdges)), vec(BLR.getΨ(m, vEdges, tEdges)); rtol=1e-8)

    refS = BLR.secondMoment(m; U=U, V=V, PA=PA, BLRAng=BLRAng, returnAvg=true, bins=60)
    gotS = BLR.secondMoment(rm; U=U, V=V, PA=PA, BLRAng=BLRAng, returnAvg=true, bins=60)
    @test isapprox(refS[4], gotS[4]; rtol=1e-8)           # line-integrated size scalar
    @test nanapprox(gotS[3], refS[3]; rtol=1e-8)          # per-channel σ²(v)

    @test_throws ArgumentError BLR.getProfile(rm, :line; bins=60, dx=ones(2, 2))
end

@testset "disk construction kernels" begin
    rMin = 300.0
    rMax = 900.0
    inc = 0.4
    nr = 16
    nϕ = 32
    f1 = 1.0
    f2 = 0.7
    f3 = 0.2
    f4 = 0.9
    αsrc = 1.2
    ηₒ = 0.4
    η₁ = 0.6
    αRM = 0.1
    rNorm = 700.0
    m = BLR.DiskWindModel(rMin, rMax, inc, nr=nr, nϕ=nϕ, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=f1, f2=f2, f3=f3, f4=f4,
        α=αsrc, ηₒ=ηₒ, η₁=η₁, αRM=αRM, rNorm=rNorm, reflect=false)

    r3D = BLR.get_r3D(inc, 0.0, 0.0)
    undoTilt = [sin(inc) 0.0 -cos(inc); 0.0 1.0 0.0; cos(inc) 0.0 sin(inc)]
    M = undoTilt * r3D
    rSystem = similar(m.camera.α)
    ϕSystem = similar(m.camera.α)
    ϕ₀ = similar(m.camera.α)
    η = similar(m.camera.α)
    xSystem = similar(m.camera.α)
    ySystem = similar(m.camera.α)
    zSystem = similar(m.camera.α)
    BLR._rt_disk_deproject!(rSystem, ϕSystem, ϕ₀, η, xSystem, ySystem, zSystem,
        m.camera.α, m.camera.β, inc, 0.0, 0.0, M, r3D, rMin, rMax, ηₒ, η₁, αRM, rNorm;
        backend=KernelAbstractions.CPU())

    _test_same_nan_approx(rSystem, _disk_field_values(m, :r))
    _test_same_nan_approx(ϕSystem, _disk_field_values(m, :ϕ))
    _test_same_nan_approx(ϕ₀, _disk_field_values(m, :ϕ₀))
    _test_same_nan_approx(η, _disk_field_values(m, :η))
    _test_same_nan_approx(xSystem, _disk_field_values(m, :x))
    _test_same_nan_approx(ySystem, _disk_field_values(m, :y))
    _test_same_nan_approx(zSystem, _disk_field_values(m, :z))

    v = similar(rSystem)
    BLR._rt_v_circular_disk!(v, rSystem, ϕSystem, inc; backend=KernelAbstractions.CPU())
    _test_same_nan_approx(v, _disk_field_values(m, :v))

    I = similar(rSystem)
    BLR._rt_disk_wind_i!(I, rSystem, ϕSystem, inc, f1, f2, f3, f4, αsrc, rMin, rMax;
        backend=KernelAbstractions.CPU())
    _test_same_nan_approx(I, _disk_field_values(m, :I))
end

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

    model_builders = (
        () -> disk() + clouds(80, 101),
        () -> clouds(80, 101) + disk(),
        () -> disk(r1=250.0, r2=700.0, nr=6, nϕ=12) + disk(r1=500.0, r2=1000.0, nr=6, nϕ=12),
        () -> disk() + clouds(50, 102) + clouds(50, 103, μ=850.0),
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

    for build in model_builders
        refModel = BLR.raytrace!(build(); τCutOff=1.0)
        backendModel = BLR.raytrace!(build(); τCutOff=1.0, backend=KernelAbstractions.CPU())
        for sym in (:I, :v, :r, :ϕ, :ϕ₀, :τ)
            @test BLR.getVariable(backendModel, sym, flatten=true) == BLR.getVariable(refModel, sym, flatten=true)
        end
        @test backendModel.subModelStartInds == refModel.subModelStartInds
        @test backendModel.camera.α == refModel.camera.α
        @test backendModel.camera.β == refModel.camera.β
    end

    ref32Model = BLR.raytrace!(model_builders[1](); IRatios=[1.0, 0.25], τCutOff=1.0)
    backend32Model = BLR.raytrace!(model_builders[1](); IRatios=[1.0, 0.25], τCutOff=1.0,
        backend=KernelAbstractions.CPU(), T=Float32)
    for sym in (:I, :v, :r)
        @test all(isapprox.(BLR.getVariable(backend32Model, sym, flatten=true),
            BLR.getVariable(ref32Model, sym, flatten=true), rtol=1e-4))
    end

    # velocity-dependent optical depth (per-point τ vector) is unimplemented and must error the
    # same way on the default CPU scan and the backend scan -- not silently produce wrong numbers.
    velModel() = begin
        d = disk()
        for r in d.rings
            r.τ = fill(0.4, length(r.I))
        end
        d + clouds(40, 109)
    end
    @test_throws ErrorException BLR.raytrace!(velModel(); τCutOff=1.0)
    @test_throws ErrorException BLR.raytrace!(velModel(); τCutOff=1.0, backend=KernelAbstractions.CPU())
end
