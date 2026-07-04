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

@testset "on-device DiskWind construction (CPU backend)" begin
    fkw = (f1=1.0, f2=0.7, f3=0.2, f4=0.9, α=1.2, ηₒ=0.4, η₁=0.6, αRM=0.1, rNorm=700.0)
    rMin, rMax, inc, nr, nϕ = 311.7, 887.3, 0.4, 20, 40
    cols = (:r, :ϕ, :ϕ₀, :i, :rot, :θₒ, :v, :I, :ΔA, :τ, :η, :x, :y, :z, :α, :β)
    # The host masks out-of-range pixels with a per-point round(r,sigdigits=9); the construction
    # kernel uses plain compares against host-rounded bounds (the plan's B2 rule). They can therefore
    # disagree only on pixels whose deprojected r sits *exactly* on rMin/rMax (the camera grid starts
    # at rMin*cos(i), so the +x meridian of the inner/outer ring lands on the boundary) — the
    # documented single-boundary-ring difference. Allow a NaN-mask mismatch only there; require
    # bit-for-bit agreement everywhere else.
    isboundary(rv, lo, hi) = isfinite(rv) && (abs(rv - lo) / lo < 1e-9 || abs(rv - hi) / hi < 1e-9)
    function _compare_cols(rm, ref, lo, hi; checkcols=cols)
        rhost = ref.r
        for c in checkcols
            a = getfield(rm.ma, c); b = getfield(ref, c)
            @test length(a) == length(b)
            for k in eachindex(a, b)
                if isnan(a[k]) != isnan(b[k])
                    @test isboundary(rhost[k], lo, hi)
                elseif isfinite(b[k])
                    @test isapprox(a[k], b[k]; atol=0.0, rtol=1e-12)
                end
            end
        end
    end
    for scale in (:linear, :log)
        m = BLR.DiskWindModel(rMin, rMax, inc; nr=nr, nϕ=nϕ, scale=scale,
            I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, fkw...)
        ref = BLR.flatten(m; T=Float64)
        rm = BLR.residentDiskWindModel(rMin, rMax, inc; nr=nr, nϕ=nϕ, scale=scale, fkw...)
        @test rm isa BLR.ResidentModel
        @test rm.nSubModels == 1
        @test length(rm.ma.I) == nr * nϕ
        _compare_cols(rm, ref, rMin, rMax)
        # reflect (per-ring scalar false) is exact everywhere — no boundary subtlety
        @test rm.ma.reflect == ref.reflect
        # drop-in for the resident observables pipeline: line profile matches host-resident
        a = BLR.getProfile(BLR.resident(m), :line; bins=50).binSums
        b = BLR.getProfile(rm, :line; bins=50).binSums
        @test all((isnan(a[i]) ? isnan(b[i]) : isapprox(a[i], b[i]; rtol=1e-8, atol=1e-12)) for i in eachindex(a, b))
    end

    # r̄/rFac/Sα parameterization matches the host DiskWindModel(r̄, rFac, Sα, i; ...) geometry
    # (both set the intensity power-law α = Sα internally, so α is not passed here)
    fkwNoα = (f1=1.0, f2=0.7, f3=0.2, f4=0.9, ηₒ=0.4, η₁=0.6, αRM=0.1, rNorm=700.0)
    r̄, rFac, Sα = 500.0, 2.0, 1.2
    loh, hih = BLR.get_rMinMaxDiskWind(r̄, rFac, Sα)
    mh = BLR.DiskWindModel(r̄, rFac, Sα, inc; nr=nr, nϕ=nϕ, scale=:log, fkwNoα...)
    rmh = BLR.residentDiskWindModel(r̄, rFac, Sα, inc; nr=nr, nϕ=nϕ, scale=:log, fkwNoα...)
    refh = BLR.flatten(mh; T=Float64)
    _compare_cols(rmh, refh, loh, hih; checkcols=(:r, :v, :I, :η, :x, :ΔA))

    # user-supplied GPU-safe scalar intensity/velocity hook runs in the same fused kernel
    Icustom = (r, ϕ, i) -> 2.0 * cos(ϕ)^2
    vcustom = (r, ϕ, i) -> 0.001 * sin(ϕ)
    rmc = BLR.residentDiskWindModel(rMin, rMax, inc; nr=nr, nϕ=nϕ, scale=:linear,
        intensity=Icustom, velocity=vcustom)
    refr = BLR.residentDiskWindModel(rMin, rMax, inc; nr=nr, nϕ=nϕ, scale=:linear, fkw...).ma.r
    for k in eachindex(rmc.ma.r)
        if isnan(refr[k])
            @test isnan(rmc.ma.I[k]) && isnan(rmc.ma.v[k])
        else
            @test rmc.ma.I[k] ≈ 2.0 * cos(rmc.ma.ϕ[k])^2
            @test rmc.ma.v[k] ≈ 0.001 * sin(rmc.ma.ϕ[k])
        end
    end

    # vCircularRadialDisk: bit-matches the host DiskWindModel(...; v=vCircularRadialDisk, vᵣFrac, inflow)
    for (vf, inf) in ((0.0, true), (0.33, true), (0.5, false))
        mr = BLR.DiskWindModel(rMin, rMax, inc; nr=nr, nϕ=nϕ, scale=:log,
            I=BLR.DiskWindIntensity, v=BLR.vCircularRadialDisk, vᵣFrac=vf, inflow=inf, fkw...)
        rmr = BLR.residentDiskWindModel(rMin, rMax, inc; nr=nr, nϕ=nϕ, scale=:log,
            vᵣFrac=vf, inflow=inf, fkw...)
        vref = BLR.flatten(mr).v
        fin = isfinite.(vref) .& isfinite.(rmr.ma.v)
        @test maximum(abs.(vref[fin] .- rmr.ma.v[fin])) < 1e-12
    end

    # error paths
    @test_throws ArgumentError BLR.residentDiskWindModel(rMin, rMax, inc; nr=nr, nϕ=nϕ)  # missing f1..f4
    @test_throws ArgumentError BLR.residentDiskWindModel(rMax, rMin, inc; nr=nr, nϕ=nϕ, fkw...)  # rMin>=rMax
    @test_throws ArgumentError BLR.residentDiskWindModel(rMin, rMax, inc; nr=1, nϕ=nϕ, fkw...)   # nr<=1
    @test_throws ArgumentError BLR.residentDiskWindModel(rMin, rMax, inc; nr=nr, nϕ=nϕ, scale=:bad, fkw...)
end

@testset "on-device cloud construction (CPU backend)" begin
    params = (μ=600.0, β=1.0, F=0.5, rₛ=1.0, θₒ=0.5, γ=1.0, ξ=0.8, i=0.4,
              ηₒ=0.4, η₁=0.6, αRM=0.1, rNorm=700.0)
    N, seed = 4000, 12345
    rm = BLR.residentCloudModel(N, seed; params...)
    @test rm isa BLR.ResidentModel
    @test rm.nSubModels == 1
    @test length(rm.ma.r) == N
    ma = rm.ma

    # cloud columns: ΔA=1, τ=0, isotropic I=1, camera α=y/β=z, i scalar = inc
    @test all(ma.ΔA .== 1.0)
    @test all(ma.τ .== 0.0)
    @test all(ma.I .== 1.0)
    @test ma.α == ma.y && ma.β == ma.z
    @test all(ma.i .== params.i)
    # geometry self-consistency: rotation + reflection are isometries -> |xyz| == r
    @test maximum(abs.(sqrt.(ma.x .^ 2 .+ ma.y .^ 2 .+ ma.z .^ 2) .- ma.r)) < 1e-9

    # deterministic transforms are BIT-EXACT vs the host package scalars given identical draws
    # (only the RNG draws differ from host :philox — validated by the KS test below).
    s = UInt64(seed); k0 = UInt32(s & 0xffffffff); k1 = UInt32((s >> 32) & 0xffffffff)
    for p in 1:150
        cidx = UInt32(p)
        ϕ₀ = BLR._cloud_uniform(Float64, k0, k1, cidx, 0) * 2π
        uθ = BLR._cloud_uniform(Float64, k0, k1, cidx, 1)
        θ = acos(cos(params.θₒ) + (1 - cos(params.θₒ)) * uθ^params.γ)
        rot = BLR._cloud_uniform(Float64, k0, k1, cidx, 2) * 2π
        g, nn = BLR._cloud_gamma(Float64, k0, k1, cidx, 1 / params.β^2, 3)
        r = params.rₛ + params.μ * params.F + params.μ * params.β^2 * (1 - params.F) * g
        xyz = BLR.rotate3D_scalar(r, ϕ₀, params.i, rot, θ, false)
        uref = BLR._cloud_uniform(Float64, k0, k1, cidx, nn)
        reflect = (xyz[3] < BLR.midPlaneXZ(xyz[1], params.i)) && (uref > params.ξ)
        xyz = reflect ? BLR.reflect_scalar(xyz[1], xyz[2], xyz[3], params.i) : xyz
        sini, cosi = sincos(params.i)
        ϕref = atan(xyz[2], sini * xyz[1] - cosi * xyz[3])
        vref = BLR.vCircularCloud(r=r, ϕ₀=ϕ₀, i=params.i, rot=rot, θₒ=θ, rₛ=params.rₛ, reflect=reflect)
        ηref = BLR.response(r; ηₒ=params.ηₒ, η₁=params.η₁, αRM=params.αRM, rNorm=params.rNorm)
        @test ma.r[p] == r
        @test ma.x[p] == xyz[1] && ma.y[p] == xyz[2] && ma.z[p] == xyz[3]
        @test ma.ϕ[p] == ϕref
        @test ma.v[p] == vref
        @test ma.η[p] == ηref
        @test ma.reflect[p] == reflect
    end

    # determinism: same seed -> identical realization; different seed -> different
    rm2 = BLR.residentCloudModel(N, seed; params...)
    @test rm2.ma.r == ma.r && rm2.ma.v == ma.v
    rm3 = BLR.residentCloudModel(N, seed + 1; params...)
    @test rm3.ma.r != ma.r

    # statistical equivalence to host :philox radii (two-sample KS), both Gamma paths
    ks2(a, b) = (A = sort(a); B = sort(b); na = length(A); nb = length(B);
        maximum(abs(searchsortedlast(A, x) / na - searchsortedlast(B, x) / nb) for x in vcat(A, B)))
    for βv in (1.0, 1.2)                     # shape=1 (direct) and shape<1 (boost)
        Nks = 20000
        rd = BLR.residentCloudModel(Nks, 77; μ=600.0, β=βv, F=0.5, θₒ=0.5, γ=1.0, ξ=0.8, i=0.4, rₛ=1.0).ma.r
        mh = BLR.cloudModel(Nks; μ=600.0, β=βv, F=0.5, θₒ=0.5, γ=1.0, ξ=0.8, i=0.4, rₛ=1.0,
            I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, rng=:philox, seed=77)
        D = ks2(rd, BLR.getVariable(mh, :r, flatten=true))
        @test D < 1.358 * sqrt(2 / Nks)      # ~95% two-sample KS critical value
    end

    # drop-in for the resident observables pipeline
    pl = BLR.getProfile(rm, :line; bins=40)
    @test count(isfinite, pl.binSums) == length(pl.binSums)
    @test isapprox(sum(filter(isfinite, pl.binSums)), Float64(N); rtol=1e-10)  # I=ΔA=1 -> Σ flux = N

    # cloudIntensity (intensity=:cloud, κ): bit-matches host given identical draws
    κ = 0.4
    rmI = BLR.residentCloudModel(N, seed; params..., intensity=:cloud, κ=κ)
    for p in 1:120
        cidx = UInt32(p)
        ϕ₀ = BLR._cloud_uniform(Float64, k0, k1, cidx, 0) * 2π
        uθ = BLR._cloud_uniform(Float64, k0, k1, cidx, 1)
        θ = acos(cos(params.θₒ) + (1 - cos(params.θₒ)) * uθ^params.γ)
        rot = BLR._cloud_uniform(Float64, k0, k1, cidx, 2) * 2π
        g, _ = BLR._cloud_gamma(Float64, k0, k1, cidx, 1 / params.β^2, 3)
        r = params.rₛ + params.μ * params.F + params.μ * params.β^2 * (1 - params.F) * g
        Iref = BLR.cloudIntensity(r=r, ϕ=0.0, θₒ=θ, ϕ₀=ϕ₀, rot=rot, i=params.i, κ=κ)
        @test rmI.ma.I[p] == Iref
    end

    # vCloudTurbulentEllipticalFlow (velocity=:turbulent): statistically equal to the host turbulent
    # model (two-sample KS on the LOS velocity); geometry still self-consistent.
    tp = (σρᵣ=0.2, σρc=0.04, σΘᵣ=0.4, σΘc=0.1, θₑ=35 / 180 * π, fEllipse=0.8, fFlow=0.0, σₜ=0.05)
    Nt = 20000
    rmT = BLR.residentCloudModel(Nt, 31; μ=600.0, β=1.0, F=0.5, rₛ=1.0, θₒ=0.5, γ=5.0, ξ=0.3, i=0.35,
        intensity=:cloud, κ=κ, velocity=:turbulent, tp...)
    @test all(isfinite, rmT.ma.v)
    @test maximum(abs.(sqrt.(rmT.ma.x .^ 2 .+ rmT.ma.y .^ 2 .+ rmT.ma.z .^ 2) .- rmT.ma.r)) < 1e-9
    mhT = BLR.cloudModel(Nt; μ=600.0, β=1.0, F=0.5, rₛ=1.0, θₒ=0.5, γ=5.0, ξ=0.3, i=0.35, κ=κ,
        I=BLR.cloudIntensity, v=BLR.vCloudTurbulentEllipticalFlow, tp..., rng=:philox, seed=31)
    @test ks2(rmT.ma.v, BLR.getVariable(mhT, :v, flatten=true)) < 1.358 * sqrt(2 / Nt)

    @test_throws ArgumentError BLR.residentCloudModel(0, seed; params...)
    @test_throws ArgumentError BLR.residentCloudModel(N, seed; params..., intensity=:bogus)
    @test_throws ArgumentError BLR.residentCloudModel(N, seed; params..., velocity=:bogus)
end

@testset "device-resident raytrace! (CPU backend)" begin
    dk(r1, r2; nr=10, nϕ=20, τ=5.0) = BLR.DiskWindModel(r1, r2, 0.4; nr=nr, nϕ=nϕ, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.5, f3=0.2, f4=0.3, α=1.0, τ=τ)
    cl(n, seed; μ=600.0, τ=0.1) = BLR.cloudModel(n; μ=μ, β=1.0, F=0.5, θₒ=0.4, i=0.4, γ=1.0, ξ=0.8,
        I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=τ, rng=:philox, seed=seed)
    fluxsum(rm) = sum(filter(isfinite, rm.ma.I .* rm.ma.ΔA))

    # the device-resident pipeline (bin→sort→segment→scan→compact, all on the backend) reproduces the
    # host raytrace! combine for every topology. Output ordering differs (active pixels then free
    # clouds), so compare order-independent aggregates: surviving point count, total flux, line profile.
    function check(builder; rfc=false, IR=1.0)
        href = BLR.resident(BLR.raytrace!(builder(); IRatios=IR, τCutOff=1.0, raytraceFreeClouds=rfc))
        rm = BLR.resident(builder(); raytrace=true)
        @test rm.rt isa BLR.RaytraceMeta
        rrt = BLR.raytrace!(rm; IRatios=IR, τCutOff=1.0, raytraceFreeClouds=rfc)
        @test rrt isa BLR.ResidentModel
        @test length(rrt.ma.I) == length(href.ma.I)
        @test isapprox(fluxsum(rrt), fluxsum(href); rtol=1e-10)
        a = BLR.getProfile(href, :line; bins=40).binSums
        b = BLR.getProfile(rrt, :line; bins=40).binSums
        @test all((isnan(a[i]) ? isnan(b[i]) : isapprox(a[i], b[i]; rtol=1e-8, atol=1e-10)) for i in eachindex(a, b))
    end
    check(() -> dk(300., 900.) + cl(300, 1))                       # disk + clouds (binned)
    check(() -> cl(300, 1) + dk(300., 900.))                       # order swapped
    check(() -> dk(250., 700.) + dk(500., 1000.))                  # N-disk overlap (grid union)
    check(() -> dk(300., 900.) + cl(200, 2) + cl(200, 3, μ=800.))  # 3 submodels
    check(() -> dk(300., 900.) + cl(300, 4); IR=[1.0, 0.25])       # non-uniform IRatios
    check(() -> cl(150, 1, μ=300., τ=2.0) + cl(150, 2, μ=320., τ=2.0); rfc=true)  # free clouds (attenuate)
    check(() -> cl(400, 5, μ=50., τ=3.0) + cl(400, 6, μ=55., τ=3.0); rfc=true)    # dense free-cloud overlaps

    # a handle with no submodels carries no metadata and cannot be device-raytraced
    @test BLR.resident(dk(300., 900.); raytrace=true).rt === nothing
    @test_throws ErrorException BLR.raytrace!(BLR.resident(dk(300., 900.)))
end

@testset "on-device-built model raytrace! (Phase 2, CPU backend)" begin
    # residentDiskWindModel/residentCloudModel carry raytrace metadata built ON-DEVICE, and `+` merges
    # it, so a model assembled with no host `ring`s raytraces end-to-end.
    rmDisk = BLR.residentDiskWindModel(300.0, 900.0, 0.4; nr=16, nϕ=32, scale=:linear,
        f1=1.0, f2=0.5, f3=0.2, f4=0.3, α=1.0)
    rmCl = BLR.residentCloudModel(400, 7; μ=600.0, β=1.0, F=0.5, θₒ=0.4, i=0.4, γ=1.0, ξ=0.8)
    @test rmDisk.rt isa BLR.RaytraceMeta && BLR._rt_meta_has_grid(rmDisk.rt)
    @test length(rmDisk.rt.grids) == 1 && rmDisk.rt.grids[1].nPixels == 16 * 32
    @test rmCl.rt isa BLR.RaytraceMeta && !BLR._rt_meta_has_grid(rmCl.rt)

    rm = rmDisk + rmCl
    @test rm.nSubModels == 2
    @test length(rm.rt.grids) == 1 && rm.rt.grids[1].nPixels == 16 * 32
    @test unique(rm.rt.submodel) == [1, 2]
    @test sum(.!rm.rt.discrete) == 16 * 32 && sum(rm.rt.discrete) == 400   # disk grid pts vs clouds

    rrt = BLR.raytrace!(rm)
    @test rrt isa BLR.ResidentModel
    @test all(isfinite, filter(!isnan, rrt.ma.I))
    @test isfinite(sum(filter(isfinite, rrt.ma.I .* rrt.ma.ΔA)))

    # the on-device disk grid edges match the host DiskWindModel grid exactly (deterministic geometry)
    mh = BLR.DiskWindModel(300.0, 900.0, 0.4; nr=16, nϕ=32, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.5, f3=0.2, f4=0.3, α=1.0) +
        BLR.cloudModel(50; μ=600.0, i=0.4, rng=:philox, seed=1)
    hmeta = BLR._rt_build_meta(mh)
    @test sort(rmDisk.rt.grids[1].grid.rMin) ≈ sort(hmeta.grids[1].grid.rMin)
    @test sort(rmDisk.rt.grids[1].grid.rMax) ≈ sort(hmeta.grids[1].grid.rMax)

    # determinism: same seed -> identical raytrace result
    rm2 = BLR.residentDiskWindModel(300.0, 900.0, 0.4; nr=16, nϕ=32, scale=:linear, f1=1.0, f2=0.5, f3=0.2, f4=0.3, α=1.0) +
        BLR.residentCloudModel(400, 7; μ=600.0, β=1.0, F=0.5, θₒ=0.4, i=0.4, γ=1.0, ξ=0.8)
    rrt2 = BLR.raytrace!(rm2)
    @test sort(rrt2.ma.I) == sort(rrt.ma.I)

    # generic N-continuous combine on-device: multiple disks merge into a multi-block grid that
    # `raytrace!` assembles innermost-first (like the host union). With matching τ, bit-reproducible
    # disks match host raytrace! to machine precision for non-overlapping / nested stacks, and to the
    # documented boundary/depth-tie tolerance where grids partially overlap (identical pixel set).
    ddk(r1, r2; τ=5.0) = BLR.residentDiskWindModel(r1, r2, 0.4; nr=10, nϕ=20, scale=:linear,
        f1=1.0, f2=0.5, f3=0.2, f4=0.3, α=1.0, τ=τ)
    hdk(r1, r2; τ=5.0) = BLR.DiskWindModel(r1, r2, 0.4; nr=10, nϕ=20, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.5, f3=0.2, f4=0.3, α=1.0, τ=τ)
    twoDisks = ddk(250.0, 700.0) + ddk(500.0, 1000.0)
    @test length(twoDisks.rt.grids) == 2                          # two grid blocks, not an error
    rdd = BLR.raytrace!(twoDisks)
    href = BLR.resident(BLR.raytrace!(hdk(250.0, 700.0) + hdk(500.0, 1000.0)))
    @test sort(filter(isfinite, rdd.ma.ΔA)) ≈ sort(filter(isfinite, href.ma.ΔA))  # identical pixel set
    fluxR = sum(filter(isfinite, rdd.ma.I .* rdd.ma.ΔA))
    fluxH = sum(filter(isfinite, href.ma.I .* href.ma.ΔA))
    @test isapprox(fluxR, fluxH; rtol=1e-3)                       # ≤ boundary/depth-tie tolerance

    # non-overlapping / nested stacks are exact
    for (db, hb) in ((() -> ddk(200., 400.) + ddk(600., 1000.), () -> hdk(200., 400.) + hdk(600., 1000.)),
                     (() -> ddk(200., 1100.) + ddk(400., 900.) + ddk(600., 800.),
                      () -> hdk(200., 1100.) + hdk(400., 900.) + hdk(600., 800.)))
        r = BLR.raytrace!(db()); h = BLR.resident(BLR.raytrace!(hb()))
        @test length(r.ma.I) == length(h.ma.I)
        @test isapprox(sum(filter(isfinite, r.ma.I .* r.ma.ΔA)), sum(filter(isfinite, h.ma.I .* h.ma.ΔA)); rtol=1e-10)
    end
end

@testset "ResidentModel device combine (+)" begin
    disk = BLR.DiskWindModel(300.0, 900.0, 0.4; nr=12, nϕ=24, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.7, f3=0.2, f4=0.9,
        α=1.2, ηₒ=0.4, η₁=0.6, αRM=0.1, rNorm=700.0)
    clouds = BLR.cloudModel(300; μ=600.0, β=1.0, F=0.5, θₒ=0.4, i=0.4, γ=1.0, ξ=0.8,
        I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, rng=:philox, seed=7)

    # resident(m1) + resident(m2) concatenates columns in the same order as resident(m1 + m2),
    # so it is bit-for-bit identical (and stays on-backend — vcat of device arrays never leaves the GPU)
    combined = BLR.resident(disk) + BLR.resident(clouds)
    reference = BLR.resident(disk + clouds)
    @test combined.nSubModels == 2
    @test combined.nSubModels == reference.nSubModels
    @test length(combined.ma.I) == length(disk.rings) * length(disk.rings[1].I) + 300
    for f in fieldnames(BLR.ModelArrays)
        @test isequal(getfield(combined.ma, f), getfield(reference.ma, f))
    end
    # observable equivalence
    @test isequal(BLR.getProfile(combined, :line; bins=40).binSums,
        BLR.getProfile(reference, :line; bins=40).binSums)

    # order matters (disk+clouds vs clouds+disk are different concatenations, same multiset of points)
    flipped = BLR.resident(clouds) + BLR.resident(disk)
    @test sort(filter(isfinite, Array(flipped.ma.r))) ≈ sort(filter(isfinite, Array(combined.ma.r)))

    # different element types cannot be combined
    @test_throws ArgumentError BLR.resident(disk; T=Float64) + BLR.resident(clouds; T=Float32)

    # mixing a host model with a resident handle errors with an actionable message (not a MethodError)
    @test_throws ArgumentError disk + BLR.resident(clouds)
    @test_throws ArgumentError BLR.resident(disk) + clouds
    err = try
        disk + BLR.resident(clouds)
    catch e
        e
    end
    @test err isa ArgumentError && occursin("resident", err.msg) && occursin("cpu(", err.msg)
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

@testset "resident composite model (W4-G1/G2, CPU backend)" begin
    # Mirrors the "composite forwarding + spectrum (W4-T4/T5/T6)" fixture in runtests.jl --
    # deliberately tiny (seconds, not minutes). Full-GPU agreement tests are W4-G4 (test/gpu_cuda.jl).
    mDisk = BLR.DiskWindModel(300.0, 900.0, 0.4; nr=8, nϕ=16, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.0, f3=0.0, f4=0.0, α=1.0,
        τ=0.4, reflect=false)
    cm = BLR.CompositeModel(mDisk; line="Ha", lineCenter=6562.8)
    BLR.addLine!(cm, mDisk; line="Hb", lineCenter=4861.3, fluxRatio=0.35)

    # --- W4-G1: resident forwarding (Float64 on the CPU backend, like the other resident CPU tests) ---
    rcm = BLR.resident(cm; backend=KernelAbstractions.CPU())
    @test rcm isa BLR.ResidentCompositeModel
    @test length(rcm) == 2
    @test rcm.lines == cm.lines
    @test rcm.lineCenters == cm.lineCenters
    @test rcm.fluxRatios == cm.fluxRatios
    for line in cm.lines
        rm = rcm[line]
        @test rm isa BLR.ResidentModel
        @test eltype(rm.ma.I) == Float64
        @test isequal(rm.ma.v, BLR.getVariable(cm[line], :v, flatten=true))
        @test isequal(rm.ma.I, BLR.getVariable(cm[line], :I, flatten=true))
    end
    # per-line kwargs forward (e.g. T=Float32)
    rcm32 = BLR.resident(cm; T=Float32)
    @test eltype(rcm32["Ha"].ma.I) == Float32 && eltype(rcm32["Hb"].ma.I) == Float32
    # getindex miss errors listing the known lines (same message style as the CPU CompositeModel)
    errIdx = try; rcm["nope"]; nothing; catch e; e; end
    @test errIdx isa ErrorException && occursin("nope", errIdx.msg) && occursin("Ha", errIdx.msg)
    # show prints one row per line
    shown = sprint(show, rcm)
    @test occursin("2 line(s)", shown) && occursin("Ha", shown) && occursin("Hb", shown)
    @test occursin("fluxRatio=0.35", shown)
    # gpu(cm) delegates per line -> without CUDA loaded the gpu(::Any) fallback error still fires
    @test_throws ErrorException BLR.gpu(cm)
end
