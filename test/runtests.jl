using BroadLineRegions
using Test
using Random

# write tests here
@testset "DiskWind model initialized successfully" begin
    mLAll = BLR.DiskWindModel(8.5e3,50.,1.,45/180*π,nr=4096,nϕ=1024,
        I=BLR.DiskWindIntensity,v=BLR.vCircularDisk,f1=1.0,f2=1.0,
        f3=1.0,f4=1.0,τ=5.,reflect=false)
    mLf1 = BLR.DiskWindModel(8.5e3,50.,1.,45/180*π,nr=4096,nϕ=1024,
        I=BLR.DiskWindIntensity,v=BLR.vCircularDisk,f1=1.0,f2=0.0,
        f3=0.0,f4=0.0,τ=5.,reflect=false)
    mLf2 = BLR.DiskWindModel(8.5e3,50.,1.,45/180*π,nr=4096,nϕ=1024,
        I=BLR.DiskWindIntensity,v=BLR.vCircularDisk,f1=0.0,f2=1.0,
        f3=0.0,f4=0.0,τ=5.,reflect=false)
    mLf3 = BLR.DiskWindModel(8.5e3,50.,1.,45/180*π,nr=4096,nϕ=1024,
        I=BLR.DiskWindIntensity,v=BLR.vCircularDisk,f1=0.0,f2=0.0,
        f3=1.0,f4=0.0,τ=5.,reflect=false)
    mLf4 = BLR.DiskWindModel(8.5e3,50.,1.,45/180*π,nr=4096,nϕ=1024,
        I=BLR.DiskWindIntensity,v=BLR.vCircularDisk,f1=0.0,f2=0.0,
        f3=0.0,f4=1.0,τ=5.,reflect=false)
    @test typeof(mLAll) == BLR.model
    LPAll = BLR.getProfile(mLAll,:line,bins=101,centered=true)
    LPf1 = BLR.getProfile(mLf1,:line,bins=101,centered=true)
    LPf2 = BLR.getProfile(mLf2,:line,bins=101,centered=true)
    LPf3 = BLR.getProfile(mLf3,:line,bins=101,centered=true)
    LPf4 = BLR.getProfile(mLf4,:line,bins=101,centered=true)
    @test maximum(LPAll.binSums) > 0.0
    @test maximum(LPf1.binSums) > 0.0
    @test maximum(LPf2.binSums) > 0.0
    @test maximum(LPf3.binSums) > 0.0
    @test maximum(LPf4.binSums) > 0.0
    @test minimum(LPAll.binSums) >= 0.0
    @test minimum(LPf1.binSums) >= 0.0
    @test minimum(LPf2.binSums) >= 0.0
    @test minimum(LPf3.binSums) >= 0.0
    @test minimum(LPf4.binSums) >= 0.0
    @test isapprox(LPAll.binCenters[findmax(LPAll.binSums)[2]],0.0,atol = 5e-3)
    @test isapprox(LPf1.binCenters[findmax(LPf1.binSums)[2]], 0.0, atol = 5e-3)
    firstMax = findmax(LPf2.binSums)[2]
    @test isapprox(LPf2.binCenters[firstMax],(0.002*sign(LPf2.binCenters[firstMax])), atol = 5e-3)
    mask = [i!=firstMax for i=1:length(LPf2.binCenters)]
    secondMax = findmax(LPf2.binSums[mask])[2]+1
    @test isapprox(LPf2.binCenters[secondMax], (0.002*sign(LPf2.binCenters[secondMax])), atol = 5e-3)
    @test isapprox(LPf3.binCenters[findmax(LPf3.binSums)[2]], 0.0, atol = 5e-3)
    firstMax = findmax(LPf4.binSums)[2]
    @test isapprox(LPf4.binCenters[firstMax], (0.002*sign(LPf2.binCenters[firstMax])), atol = 5e-3)
    mask = [i!=firstMax for i=1:length(LPf4.binCenters)]
    secondMax = findmax(LPf4.binSums[mask])[2]+1
    @test isapprox(LPf4.binCenters[secondMax], (0.002*sign(LPf2.binCenters[secondMax])), atol = 5e-3)
    M = 1.7e8*2e30; rs = 2*6.67e-11*M/9e16; rsDay = rs/3e8/24/3600
    tCenters,Ψt = BLR.getΨt(mLAll,501,100/rsDay)
    @test isapprox(tCenters[findmax(Ψt)[2]]*rsDay, 40.0, atol = 5e-1)
end

@testset "Cloud model initialized successfully" begin 
    M = 10^(6.5)*2e30; rs = 2*6.67e-11*M/9e16; rsDay = rs/3e8/24/3600
    mP1 = BLR.cloudModel(1_000_000,μ=4/rsDay,F=0.25,β=1.0,θₒ=40/180*π,i=20/180*π,
               κ=-0.4,γ=5.0,ξ=0.3,fEllipse=0.0,fFlow=0.0,θₑ=0.0,σₜ=0.0,
               σρc=0.0,σΘᵣ=0.0,σΘc=0.0,σρᵣ=0.0,
               I=BLR.cloudIntensity,v=BLR.vCloudTurbulentEllipticalFlow,τ=0.0)
    @test typeof(mP1) == BLR.model
    LP1 = BLR.getProfile(mP1,:line,bins=101,centered=true)
    @test maximum(LP1.binSums) > 0.0
    @test minimum(LP1.binSums) >= 0.0
    @test isapprox(LP1.binCenters[findmax(LP1.binSums)[2]],0.003,atol = 5e-3)
    mP2 = BLR.cloudModel(1_000_000,μ=4/rsDay,F=0.25,β=0.8,θₒ=30/180*π,i=20/180*π,
        κ=-0.4,γ=5.0,ξ=0.1,fEllipse=1.0,σₜ=0.0,
        fFlow=0.0,θₑ=0.0,σρc=0.0,σΘᵣ=0.0,σΘc=0.0,σρᵣ=0.0,
        I=BLR.cloudIntensity,v=BLR.vCloudTurbulentEllipticalFlow,τ=0.0)
    LP2 = BLR.getProfile(mP2,:line,bins=101,centered=true)
    @test maximum(LP2.binSums) > 0.0
    @test minimum(LP2.binSums) >= 0.0
    firstMax = findmax(LP2.binSums)[2]
    @test isapprox(LP2.binCenters[firstMax], (0.001*sign(LP2.binCenters[firstMax])), atol = 5e-3)
    mask = sign(LP2.binCenters[firstMax]) == 1 ? [i < 51 for i=1:length(LP2.binCenters)] : [i > 50 for i=1:length(LP2.binCenters)]
    secondMax = findmax(LP2.binSums[mask])[2]
    @test isapprox(LP2.binCenters[mask][secondMax],(0.001*sign(LP2.binCenters[mask][secondMax])),atol = 5e-3)
    tCenters,Ψt = BLR.getΨt(mP1,501,10/rsDay)
    @test isapprox(tCenters[findmax(Ψt)[2]]*rsDay, 1.3, atol = 5e-1)
    tCenters,Ψt = BLR.getΨt(mP2,501,10/rsDay)
    @test isapprox(tCenters[findmax(Ψt)[2]]*rsDay, 1.8, atol = 5e-1)
end

function _cloud_field_tuple(m)
    return (
        BLR.getVariable(m, :r, flatten=true),
        BLR.getVariable(m, :ϕ, flatten=true),
        BLR.getVariable(m, :ϕ₀, flatten=true),
        BLR.getVariable(m, :rot, flatten=true),
        BLR.getVariable(m, :θₒ, flatten=true),
        BLR.getVariable(m, :x, flatten=true),
        BLR.getVariable(m, :y, flatten=true),
        BLR.getVariable(m, :z, flatten=true),
        BLR.getVariable(m, :I, flatten=true),
        BLR.getVariable(m, :v, flatten=true),
    )
end

function _ks_statistic(x, y)
    xs = sort(x)
    ys = sort(y)
    i = 1
    j = 1
    nx = length(xs)
    ny = length(ys)
    d = 0.0
    while i <= nx || j <= ny
        if j > ny || (i <= nx && xs[i] <= ys[j])
            v = xs[i]
        else
            v = ys[j]
        end
        while i <= nx && xs[i] <= v
            i += 1
        end
        while j <= ny && ys[j] <= v
            j += 1
        end
        d = max(d, abs((i - 1) / nx - (j - 1) / ny))
    end
    return d
end

@testset "Philox cloud generator" begin
    args = (; μ=500.0, β=0.8, F=0.25, θₒ=40/180*π, i=20/180*π,
        γ=2.0, ξ=0.3, I=BLR.cloudIntensity, v=BLR.vCircularCloud, κ=-0.2, τ=0.0)

    rng1 = MersenneTwister(123)
    mt = BLR.cloudModel(256; rng=rng1, args...)
    rng2 = MersenneTwister(123)
    ϕ₀ = rand(rng2, 256) .* 2π
    θ = acos.(cos(args.θₒ) .+ (1 - cos(args.θₒ)) .* rand(rng2, 256).^args.γ)
    rot = rand(rng2, 256) .* 2π
    mtRef = BLR.cloudModel(ϕ₀, ones(256) .* args.i, rot, θ, args.θₒ, args.ξ;
        rₛ=1.0, μ=args.μ, β=args.β, F=args.F, I=args.I, v=args.v,
        rng=rng2, κ=args.κ, τ=args.τ)
    @test _cloud_field_tuple(mt) == _cloud_field_tuple(mtRef)

    philox1 = BLR.cloudModel(512; rng=:philox, seed=9876, parallel=true, args...)
    philox2 = BLR.cloudModel(512; rng=:philox, seed=9876, parallel=true, args...)
    philoxSeq = BLR.cloudModel(512; rng=:philox, seed=9876, parallel=false, args...)
    philoxOther = BLR.cloudModel(512; rng=:philox, seed=9877, parallel=false, args...)
    @test _cloud_field_tuple(philox1) == _cloud_field_tuple(philox2)
    @test _cloud_field_tuple(philox1) == _cloud_field_tuple(philoxSeq)
    @test _cloud_field_tuple(philox1)[1] != _cloud_field_tuple(philoxOther)[1]

    ringsRev = [BLR._drawPhiloxCloud(j, 9876; args...) for j in 512:-1:1]
    philoxRev = BLR.model(reverse(ringsRev))
    @test _cloud_field_tuple(philoxSeq) == _cloud_field_tuple(philoxRev)

    mtR = BLR.getVariable(BLR.cloudModel(4096; rng=MersenneTwister(4321), args...), :r, flatten=true)
    philoxR = BLR.getVariable(BLR.cloudModel(4096; rng=:philox, seed=4321, parallel=false, args...), :r, flatten=true)
    d = _ks_statistic(mtR, philoxR)
    @test d < 1.63 * sqrt((length(mtR) + length(philoxR)) / (length(mtR) * length(philoxR)))

    @test_throws ArgumentError BLR.cloudModel(4; rng=:philox, args...)
    @test_throws ArgumentError BLR.cloudModel(4; rng=:unknown, args...)
    @test_throws ArgumentError BLR.cloudModel(4; rng=MersenneTwister(1), seed=1, args...)
end

@testset "second moment profiles" begin
    #binnedVariance machinery on hand-built arrays
    edges = [0.0,0.5,1.0]
    x = [0.1,0.1,0.9]; θ = [1.0,3.0,5.0]; w = [1.0,1.0,2.0]
    _,_,σ² = BLR.binnedVariance(x,θ,w,bins=edges)
    @test σ²[1] ≈ 1.0
    @test σ²[2] ≈ 0.0 #single point in bin -> point source has no size
    _,_,σ² = BLR.binnedVariance(x,θ,w,bins=[0.0,0.4,0.6,1.0])
    @test isnan(σ²[2]) #empty bins are NaN, not 0
    #pass 2 follows binnedSum's edge rules: interior edge point goes to the bin below
    _,_,σ² = BLR.binnedVariance([0.1,0.5],[0.0,2.0],[1.0,1.0],bins=edges)
    @test σ²[1] ≈ 1.0 && isnan(σ²[2])
    #points exactly on the first/last edges (or beyond) only counted with overflow=true
    _,_,σ² = BLR.binnedVariance([0.0,0.1,1.0,0.9],[0.0,2.0,1.0,3.0],[1.0,1.0,1.0,1.0],bins=edges)
    @test σ²[1] ≈ 0.0 && σ²[2] ≈ 0.0
    _,_,σ² = BLR.binnedVariance([-1.0,0.1,2.0,0.9],[0.0,2.0,1.0,3.0],[1.0,1.0,1.0,1.0],bins=edges,overflow=true)
    @test σ²[1] ≈ 1.0 && σ²[2] ≈ 1.0

    #model-level checks
    m = BLR.DiskWindModel(8.5e3,50.,1.,45/180*π,nr=128,nϕ=256,
        I=BLR.DiskWindIntensity,v=BLR.vCircularDisk,f1=1.0,f2=1.0,
        f3=1.0,f4=1.0,τ=5.,reflect=false)
    U = [10.0,30.0]; V = [20.0,-40.0]; PA = 0.7; BLRAng = 1e-10
    res = BLR.secondMoment(m,U=U,V=V,PA=PA,BLRAng=BLRAng,bins=101,centered=true)
    @test length(res) == 2
    edges,centers,σ²,σ²tot = res[1]
    finite = isfinite.(σ²)
    @test any(finite)
    @test all(σ²[finite] .>= 0.0)
    @test σ²tot > 0.0
    #angular size has the baseline length divided out -- scaling both U and V leaves σ² unchanged
    res2 = BLR.secondMoment(m,U=2 .*U,V=2 .*V,PA=PA,BLRAng=BLRAng,bins=101,centered=true)
    @test all(isapprox.(res2[1][3][finite],σ²[finite],rtol=1e-12))
    @test isapprox(res2[1][4],σ²tot,rtol=1e-12)
    #σ² scales as BLRAng²
    res4 = BLR.secondMoment(m,U=U,V=V,PA=PA,BLRAng=2*BLRAng,bins=101,centered=true)
    @test all(isapprox.(res4[1][3][finite],4 .*σ²[finite],rtol=1e-10))
    @test isapprox(res4[1][4],4*σ²tot,rtol=1e-10)
    #line-integrated size matches an independent direct computation over all points
    U′ = cos(PA)*U[1]+sin(PA)*V[1]; V′ = -sin(PA)*U[1]+cos(PA)*V[1]
    s = vec(m.camera.α.*BLRAng.*U′ .+ m.camera.β.*BLRAng.*V′)./hypot(U′,V′)
    wts = BLR.getVariable(m,:I,flatten=true).*BLR.getVariable(m,:ΔA,flatten=true)
    vflat = BLR.getVariable(m,:v,flatten=true)
    mask = isfinite.(vflat) .& isfinite.(s.*wts)
    s̄ = sum(wts[mask].*s[mask])/sum(wts[mask])
    σ²direct = sum(wts[mask].*(s[mask].-s̄).^2)/sum(wts[mask])
    @test isapprox(σ²tot,σ²direct,rtol=1e-10)
    #law of total variance: line-integrated size exceeds the flux-weighted mean of per-channel sizes
    resOF = BLR.secondMoment(m,U=U,V=V,PA=PA,BLRAng=BLRAng,bins=101,centered=true,overflow=true)
    LP = BLR.getProfile(m,:line,bins=101,centered=true,overflow=true)
    fOF = isfinite.(resOF[1][3])
    meanChannelSize = sum(LP.binSums[fOF].*resOF[1][3][fOF])/sum(LP.binSums[fOF])
    @test resOF[1][4] > meanChannelSize #photocenter spread term is the rotation signature
    #returnAvg mirrors phase conventions
    edgesAvg,centersAvg,σ²Avg,σ²totAvg = BLR.secondMoment(m,U=U,V=V,PA=PA,BLRAng=BLRAng,bins=101,centered=true,returnAvg=true)
    @test length(σ²Avg) == length(σ²)
    @test σ²totAvg ≈ (res[1][4]+res[2][4])/2
    #getProfile hook
    p = BLR.getProfile(m,:moment2,U=U,V=V,PA=PA,BLRAng=BLRAng,bins=101,centered=true)
    @test typeof(p) == BLR.profile
    @test p.name == :moment2
    @test all(isapprox.(p.binSums[finite],σ²Avg[finite],rtol=1e-12))
    BLR.setProfile!(m,p)
    @test :moment2 ∈ keys(m.profiles)
    #getProfile forwards its bins argument to phase (centered=false so the bin count is not recomputed)
    pPhase = BLR.getProfile(m,:phase,U=U,V=V,PA=PA,BLRAng=BLRAng,bins=51,centered=false)
    @test length(pPhase.binSums) == 51
end

@testset "second moment analytic ring checks" begin
    #a thin uniform ring of radius r₀ at inclination i has closed-form projected second moments:
    #σ²(ψ) = ⟨r²⟩BLRAng²(cos²ψ + cos²i·sin²ψ)/2 with ψ the sky angle from the node line.
    #⟨r²⟩ comes from the model's intrinsic radii/weights, so these checks are independent of the
    #camera coordinates and binning machinery that secondMoment uses
    i60 = 60/180*π; BLRAng = 1e-10
    mRing = BLR.DiskWindModel(990.,1010.,i60,nr=256,nϕ=512,I=BLR.IsotropicIntensity,v=BLR.vCircularDisk)
    rFlat = BLR.getVariable(mRing,:r,flatten=true)
    wFlat = BLR.getVariable(mRing,:I,flatten=true).*BLR.getVariable(mRing,:ΔA,flatten=true)
    ringMask = isfinite.(rFlat.*wFlat)
    r̄² = sum(wFlat[ringMask].*rFlat[ringMask].^2)/sum(wFlat[ringMask])
    A = r̄²*BLRAng^2/2 #variance along the node line (sky major axis)
    χ = collect(range(0,π,length=61))[1:end-1] #baseline sky angles
    σ²s = [BLR.secondMoment(mRing,U=[50*cos(c)],V=[50*sin(c)],PA=0.0,BLRAng=BLRAng,bins=11)[1][4] for c in χ]
    @test isapprox(maximum(σ²s), A, rtol=2e-2)
    @test isapprox(minimum(σ²s), A*cos(i60)^2, rtol=2e-2)
    #trace invariance: any two orthogonal baselines sum to ⟨r²⟩BLRAng²(1+cos²i)/2
    @test all(isapprox.(σ²s[1:30].+σ²s[31:60], A*(1+cos(i60)^2), rtol=1e-2))
    #rotating the model by PA equals rotating the baselines by -PA
    σ²PA = BLR.secondMoment(mRing,U=[50*cos(1.0)],V=[50*sin(1.0)],PA=0.7,BLRAng=BLRAng,bins=11)[1][4]
    σ²rot = BLR.secondMoment(mRing,U=[50*cos(0.3)],V=[50*sin(0.3)],PA=0.0,BLRAng=BLRAng,bins=11)[1][4]
    @test isapprox(σ²PA, σ²rot, rtol=1e-10)
    #per-channel: with v ∝ sinφ the isovelocity pair (φ, π-φ) shares its node-line coordinate, so
    #channel images are two points split ⊥ to the node line: σ⊥(v) = √⟨r²⟩cos(i)√(1-(v/vmax)²)
    #while σ∥ ≈ 0 (bin-width smearing only)
    θ₀ = χ[argmax(σ²s)] #node-line direction
    resP = BLR.secondMoment(mRing,U=[50*cos(θ₀+π/2)],V=[50*sin(θ₀+π/2)],PA=0.0,BLRAng=BLRAng,bins=41,centered=true)[1]
    resN = BLR.secondMoment(mRing,U=[50*cos(θ₀)],V=[50*sin(θ₀)],PA=0.0,BLRAng=BLRAng,bins=41,centered=true)[1]
    k0 = argmin(abs.(resP[2])) #channel at line center
    @test isapprox(sqrt(resP[3][k0]), sqrt(r̄²)*cos(i60)*BLRAng, rtol=2e-2)
    @test resN[3][k0] < 1e-2*resP[3][k0]
end

@testset "cached xyz coordinates (getXYZ)" begin
    #disk model: cache populated at construction, NaN-masked consistently by removeNaN!
    mD = BLR.DiskWindModel(500., 5., 1., 30/180*π, nr=32, nϕ=64, scale=:log,
        f1=1.0, f2=1.0, f3=1.0, f4=1.0, reflect=false, τ=5.)
    for r in mD.rings
        @test !isnothing(r.x) && !isnothing(r.y) && !isnothing(r.z)
        @test length(r.x) == length(r.I)
        @test isnan.(r.x) == isnan.(r.I) #NaN sentinels must line up with intensity
    end
    BLR.removeNaN!(mD)
    for r in mD.rings
        @test length(r.x) == length(r.I) && length(r.y) == length(r.I) && length(r.z) == length(r.I)
        @test !any(isnan, r.x)
    end
    #spot-check: cached values match a fresh rotate3D computation
    #(after removeNaN! per-ring scalars may have been expanded to vectors -- handle both, see removeNaN!)
    r1 = mD.rings[1]
    k = 1
    i1 = typeof(r1.i) == Float64 ? r1.i : r1.i[k]
    rot1 = typeof(r1.rot) == Float64 ? r1.rot : r1.rot[k]
    θₒ1 = typeof(r1.θₒ) == Float64 ? r1.θₒ : r1.θₒ[k]
    reflect1 = typeof(r1.reflect) == Bool ? r1.reflect : r1.reflect[k]
    fresh = BLR.rotate3D_scalar(r1.r[k], r1.ϕ₀[k], i1, rot1, θₒ1, reflect1)
    @test isapprox(r1.x[k], fresh[1], rtol=1e-12) && isapprox(r1.y[k], fresh[2], rtol=1e-12) && isapprox(r1.z[k], fresh[3], rtol=1e-12)

    #cloud model: cache populated by drawCloud (post-reflection), camera built from it
    mC = BLR.cloudModel(2_000, μ=500., β=1.0, F=0.5, θₒ=40/180*π, i=20/180*π, γ=1.0, ξ=0.5,
        I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=0.0, rng=MersenneTwister(7))
    for r in mC.rings[1:100]
        @test typeof(r.x) == Float64
        fresh = BLR.rotate3D_scalar(r.r, r.ϕ₀, r.i, r.rot, r.θₒ, r.reflect)
        @test isapprox(r.x, fresh[1], rtol=1e-12) && isapprox(r.y, fresh[2], rtol=1e-12) && isapprox(r.z, fresh[3], rtol=1e-12)
    end
    #camera coordinates must equal cached (y, z)
    @test all(mC.camera.α[j] == mC.rings[j].y for j in 1:length(mC.rings))
    @test all(mC.camera.β[j] == mC.rings[j].z for j in 1:length(mC.rings))
    #delay spot-check: tCloud must agree with direct recomputation
    rc = mC.rings[1]
    @test isapprox(BLR.tCloud(rc), rc.η*(rc.r - BLR.rotate3D_scalar(rc.r, rc.ϕ₀, rc.i, rc.rot, rc.θₒ, rc.reflect)[1]), rtol=1e-12)

    #lazy path: a hand-built ring without x/y/z computes and caches on first access
    rl = BLR.ring(r=100.0, i=0.3, v=0.001, I=1.0, ϕ=0.5, ϕ₀=0.5, ΔA=1.0, rot=0.2, θₒ=0.1, reflect=true)
    @test isnothing(rl.x)
    xyz = BLR.getXYZ(rl)
    @test !isnothing(rl.x) && rl.x == xyz[1] && rl.y == xyz[2] && rl.z == xyz[3]
    @test all(collect(xyz) .== BLR.rotate3D(100.0, 0.5, 0.3, 0.2, 0.1, true)) #reflection already applied exactly once
end

@testset "getVariable cache" begin
    m = BLR.DiskWindModel(500., 5., 1., 30/180*π, nr=16, nϕ=32, scale=:log,
        f1=1.0, f2=1.0, f3=1.0, f4=1.0, reflect=false, τ=5.)
    v1 = BLR.getVariable(m, :I)
    @test BLR.getVariable(m, :I) === v1 #cache hit returns the memoized array
    vf = BLR.getVariable(m, :I, flatten=true)
    @test isequal(vf, vec(v1)) #flattened variant cached under its own key but consistent
    #invalidation contract: reset! after direct ring mutation
    m.rings[1].I[1] = 12345.6789
    BLR.reset!(m)
    v2 = BLR.getVariable(m, :I)
    @test v2 !== v1
    @test v2[1,1] == 12345.6789
    #package delay functions are memoized; user functions are always re-evaluated
    d1 = BLR.getVariable(m, BLR.t)
    @test BLR.getVariable(m, BLR.t) === d1
    myfun(ring) = ring.r
    @test BLR.getVariable(m, myfun) !== BLR.getVariable(m, myfun)
    @test isequal(BLR.getVariable(m, myfun), BLR.getVariable(m, :r))
    #combining models starts with a fresh cache
    mc = m + m
    @test isempty(mc.cache)
    #cache can be disabled entirely
    m.cache = nothing
    @test BLR.getVariable(m, :I) !== BLR.getVariable(m, :I)
    #removeNaN! re-enables and repopulates consistently
    m2 = BLR.DiskWindModel(500., 5., 1., 30/180*π, nr=16, nϕ=32, scale=:log,
        f1=1.0, f2=1.0, f3=1.0, f4=1.0, reflect=false, τ=5.)
    BLR.getVariable(m2, :I) #populate
    BLR.removeNaN!(m2)
    @test !any(isnan, BLR.getVariable(m2, :I, flatten=true)) #post-removal gather is fresh, not stale
end

@testset "combined model profiles + gathers" begin
    mD = BLR.DiskWindModel(500., 5., 1., 45/180*π, nr=32, nϕ=64, scale=:log,
        f1=1.0, f2=1.0, f3=1.0, f4=1.0, reflect=false, τ=5.)
    mC = BLR.cloudModel(5_000, μ=500., β=1., F=0.5, θₒ=30/180*π, γ=1., ξ=1., i=45/180*π,
        I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, rescale=1e-5, τ=0.0, rng=MersenneTwister(11))
    mc = mD + mC
    #fast paths must match the per-ring closure computation exactly (combined delays = tCloud everywhere)
    pFast = BLR.getProfile(mc, :delay, bins=31, centered=true)
    d(ring) = BLR.tCloud(ring).*ring.I
    den(ring) = (typeof(ring.r) == Float64 && typeof(ring.ϕ) == Float64) ? ring.I*ring.η : ring.I.*ring.η
    pNum = BLR.binModel(31, m=mc, yVariable=d, xVariable=:v, centered=true)
    pDen = BLR.binModel(31, m=mc, yVariable=den, xVariable=:v, centered=true)
    @test isequal(pFast.binSums, pNum[3]./pDen[3])
    pR = BLR.getProfile(mc, :r, bins=31, centered=true)
    rfun(ring) = ring.r.*ring.I
    @test isequal(pR.binSums, BLR.binModel(31,m=mc,yVariable=rfun,xVariable=:v,centered=true)[3] ./
                              BLR.binModel(31,m=mc,yVariable=:I,xVariable=:v,centered=true)[3])
    #delay arrays are memoized for combined models too
    @test BLR.getVariable(mc, BLR.t) === BLR.getVariable(mc, BLR.t)
    #tCloud (general r - x) reduces to the analytic tDisk formula on flat rings
    for k in (1, 16, 32)
        td = BLR.tDisk(mD.rings[k]); tg = BLR.tCloud(mD.rings[k])
        fin = .!isnan.(td)
        @test maximum(abs.((td[fin] .- tg[fin]) ./ td[fin])) < 1e-10
    end
    #combined transfer functions are order-independent and work clouds-first (errored before)
    mc2 = mC + mD #clouds first
    tEdges = collect(range(0.0, stop=2000.0, length=101))
    Ψ1 = BLR.getΨt(mc, tEdges); Ψ2 = BLR.getΨt(mc2, tEdges)
    @test all(isapprox.(Ψ1, Ψ2, rtol=1e-10))
    #delay profile and transfer function use identical delays for combined models
    delays = BLR.getVariable(mc, BLR.t)
    dPer = vcat([BLR.tCloud(r) for r in mc.rings[mc.subModelStartInds[2]:end]]...)
    nDisk = length(delays) - length(dPer)
    @test isequal(delays[nDisk+1:end], dPer)
    #expandPerPoint aligns per-ring scalar fields (cloud η here) with the flattened I
    ηpp = BLR.expandPerPoint(mc, :η)
    @test length(ηpp) == length(BLR.getVariable(mc, :I, flatten=true))
    #Function gathers on >2 submodels match Symbol gathers (regression: chunk offsets were wrong and errored before)
    m3 = mc + BLR.cloudModel(1_000, μ=500., β=1., F=0.5, θₒ=30/180*π, γ=1., ξ=1., i=45/180*π,
        I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, rescale=1e-5, τ=0.0, rng=MersenneTwister(12))
    getr(ring) = ring.r
    @test isequal(BLR.getVariable(m3, getr), BLR.getVariable(m3, :r))

    #DEFINITION (combined models): Ψ(t) is the histogram of I⋅ΔA over the general tCloud delays
    #t = η(r − x), per point, for every submodel. Reference computed ring-by-ring from first
    #principles -- seed- and version-independent, so this pins the convention without magic numbers.
    dRef = Float64[]; iRef = Float64[]; aRef = Float64[]; vRef = Float64[]; ηRef = Float64[]
    for r in mc.rings
        td = BLR.tCloud(r)
        n = length(r.I)
        append!(dRef, td isa Number ? (td,) : td)
        append!(iRef, r.I isa Number ? (r.I,) : r.I)
        append!(vRef, r.v isa Number ? (r.v,) : r.v)
        append!(aRef, r.ΔA isa Number ? (r.ΔA for _ in 1:n) : r.ΔA)
        append!(ηRef, r.η isa Number ? (r.η for _ in 1:n) : r.η)
    end
    ΨRef = map(1:length(tEdges)-1) do j
        mask = (dRef .> tEdges[j]) .& (dRef .<= tEdges[j+1]) #binnedSum convention (left-exclusive)
        s = sum(iRef[mask].*aRef[mask])
        s > 0 ? s : 1e-30 #replicate getΨt's empty-bin floor
    end
    @test all(isapprox.(Ψ1, ΨRef, rtol=1e-10))
    #and the :delay profile is the intensity(×η)-weighted mean of those same tCloud delays per velocity bin
    vEdges = collect(range(-0.05, stop=0.05, length=22))
    pD = BLR.getProfile(mc, :delay, bins=vEdges, centered=false)
    dProfRef = map(1:length(vEdges)-1) do j
        mask = (vRef .>= vEdges[j]) .& (vRef .< vEdges[j+1]) .& isfinite.(dRef) .& isfinite.(iRef)
        sum(dRef[mask].*iRef[mask].*aRef[mask]) / sum(iRef[mask].*ηRef[mask].*aRef[mask])
    end
    fin = isfinite.(pD.binSums) .& isfinite.(dProfRef)
    @test count(fin) > 10 #most velocity bins should be populated
    @test all(isapprox.(pD.binSums[fin], dProfRef[fin], rtol=1e-8))
end

@testset "binnedSum direct indexing" begin
    rng = MersenneTwister(99)
    nb = 37; lo = -2.0; hi = 3.0
    edges = collect(range(lo, stop=hi, length=nb+1))
    #random values + values exactly on every edge + NaN/Inf + out-of-range
    xs = vcat(lo .+ (hi-lo).*rand(rng, 10_000), copy(edges), [NaN, Inf, -Inf, lo-1.0, hi+1.0])
    ys = vcat(rand(rng, 10_000), ones(length(edges)), ones(5))
    for ovf in (false, true)
        #Int bins (uniform fast path) vs the same edges as an explicit vector (searchsortedfirst path)
        e1, c1, r1 = BLR.binnedSum(xs, ys, bins=nb, centered=false, minX=lo, maxX=hi, overflow=ovf)
        e2, c2, r2 = BLR.binnedSum(xs, ys, bins=edges, overflow=ovf)
        @test e1 == e2
        @test r1 == r2 #bit-identical, including exact-edge assignment convention
    end
end

@testset "transfer functions single-pass" begin
    #reference = the old masked algorithm, written out explicitly; new code must match
    #(bin assignment exactly; sums to summation-order rounding)
    mD = BLR.DiskWindModel(500., 5., 1., 30/180*π, nr=16, nϕ=32, scale=:log,
        f1=1.0, f2=1.0, f3=1.0, f4=1.0, reflect=false, τ=5.) #keeps NaN sentinel pixels -- exercises NaN exclusion
    mC = BLR.cloudModel(2_000, μ=500., β=1.0, F=0.5, θₒ=40/180*π, i=20/180*π, γ=1.0, ξ=0.5,
        I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=0.0, rng=MersenneTwister(21))
    for m in (mD, mC, mD + mC)
        I = BLR.getVariable(m,:I); ΔA = BLR.getVariable(m,:ΔA); v = BLR.getVariable(m,:v); d = BLR.getVariable(m,BLR.t)
        vEdges = collect(range(-0.05, stop=0.05, length=12))
        tEdges = collect(range(0.0, stop=1500.0, length=9))
        Ψ = BLR.getΨ(m, vEdges, tEdges)
        for i in 1:length(vEdges)-1, j in 1:length(tEdges)-1
            #binnedSum convention: left-exclusive / right-inclusive interior edges (> left, <= right)
            mask = (v .> vEdges[i]) .& (v .<= vEdges[i+1]) .& (d .> tEdges[j]) .& (d .<= tEdges[j+1])
            s = sum(I[mask].*ΔA[mask])
            expected = s > 0 ? s : 1e-30
            @test isapprox(Ψ[i,j], expected, rtol=1e-12)
            @test (Ψ[i,j] == 1e-30) == (expected == 1e-30) #empty/poisoned bins identical
        end
        for ovf in (false, true)
            Ψt = BLR.getΨt(m, tEdges, ovf)
            ref = map(1:length(tEdges)-1) do j
                mask = (d .> tEdges[j]) .& (d .<= tEdges[j+1])
                s = sum(I[mask].*ΔA[mask])
                s > 0 ? s : 1e-30
            end
            if ovf
                sU = sum(I[(d .<= tEdges[1])].*ΔA[(d .<= tEdges[1])])
                sO = sum(I[(d .>= tEdges[end])].*ΔA[(d .>= tEdges[end])])
                ref[1] += sU > 0 ? sU : 1e-30
                ref[end] += sO > 0 ? sO : 1e-30
            end
            @test all(isapprox.(Ψt, ref, rtol=1e-12))
        end
    end
end

@testset "+ associativity (W4-T0)" begin
    # Bug being fixed: Base.:+ used to push! a single boundary for the whole
    # right-hand operand, discarding its internal subModelStartInds -- so
    # a+(b+c) collapsed to 2 slots vs 3 for (a+b)+c. Models here are kept
    # deliberately tiny; this testset should add seconds, not minutes.
    d(; r1=300.0, r2=900.0, nr=3, nϕ=6, inc=0.4, τ=0.4) =
        BLR.DiskWindModel(r1, r2, inc, nr=nr, nϕ=nϕ, scale=:linear,
            I=BLR.IsotropicIntensity, v=BLR.vCircularDisk, τ=τ, reflect=false)
    c(n, seed; μ=600.0, inc=0.4, τ=0.1) =
        BLR.cloudModel(n, μ=μ, β=1.0, F=0.5, θₒ=0.4, i=inc, γ=1.0, ξ=1.0,
            I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=τ, rng=:philox, seed=seed)

    a = d(r1=300.0, r2=500.0)
    b = c(20, 201)
    e = d(r1=600.0, r2=800.0)

    left = (a + b) + e
    right = a + (b + e)

    # (a) associativity: identical submodel boundaries (3 slots, not 2) and
    # identical flattened camera indices.
    @test length(left.subModelStartInds) == 3
    @test left.subModelStartInds == right.subModelStartInds
    @test BLR.getFlattenedCameraIndices(left) == BLR.getFlattenedCameraIndices(right)

    # raytrace! returns a NEW model rather than mutating -- capture the return.
    rtLeft = BLR.raytrace!(left)
    rtRight = BLR.raytrace!(right)
    LPLeft = BLR.getProfile(rtLeft, :line, bins=11, centered=true)
    LPRight = BLR.getProfile(rtRight, :line, bins=11, centered=true)
    @test LPLeft.binSums == LPRight.binSums
    @test LPLeft.binCenters == LPRight.binCenters

    # (b) mixed right-nest disk+(cloud+disk) matches its left-associated
    # equivalent (a,b,e above already form exactly this disk/cloud/disk mix).
    mixedRight = a + (b + e)
    mixedLeft = (a + b) + e
    @test mixedRight.subModelStartInds == mixedLeft.subModelStartInds
    rtMixedRight = BLR.raytrace!(mixedRight)
    rtMixedLeft = BLR.raytrace!(mixedLeft)
    LPmr = BLR.getProfile(rtMixedRight, :line, bins=11, centered=true)
    LPml = BLR.getProfile(rtMixedLeft, :line, bins=11, centered=true)
    @test LPmr.binSums == LPml.binSums
    @test LPmr.binCenters == LPml.binCenters
end

@testset "model params + rebuild (W4-T1)" begin
    # Models here are deliberately tiny -- this testset should add seconds, not minutes.
    cloudArgs = (; μ=600.0, β=1.0, F=0.5, θₒ=0.4, i=0.4, γ=1.0, ξ=1.0,
        I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=0.1)

    # (a) rMin/rMax-form DiskWindModel records :DiskWindModel; rebuild -> bit-identical :line binSums
    mDW = BLR.DiskWindModel(300.0, 900.0, 0.4; nr=8, nϕ=16, scale=:linear,
        I=BLR.IsotropicIntensity, v=BLR.vCircularDisk, τ=0.4, reflect=false)
    @test mDW.params.constructor == :DiskWindModel
    @test mDW.params.rMin == 300.0 && mDW.params.rMax == 900.0
    mDWr = BLR.rebuild(mDW)
    @test BLR.getProfile(mDW, :line, bins=21, centered=true).binSums ==
          BLR.getProfile(mDWr, :line, bins=21, centered=true).binSums

    # mean-form (:DiskWindModelMean) also records provenance and round-trips
    mMean = BLR.DiskWindModel(500.0, 5.0, 1.0, 0.4; nr=8, nϕ=16, scale=:log,
        I=BLR.IsotropicIntensity, v=BLR.vCircularDisk, τ=0.4)
    @test mMean.params.constructor == :DiskWindModelMean
    @test BLR.getProfile(mMean, :line, bins=21, centered=true).binSums ==
          BLR.getProfile(BLR.rebuild(mMean), :line, bins=21, centered=true).binSums

    # (b) philox cloudModel rebuilt from params is bit-identical
    mC = BLR.cloudModel(64; rng=:philox, seed=42, cloudArgs...)
    @test mC.params.constructor == :cloudModel
    @test mC.params.rng == :philox && mC.params.seed == 42
    @test _cloud_field_tuple(mC) == _cloud_field_tuple(BLR.rebuild(mC))

    # (c) legacy-rng cloud model records rng=:legacy
    mLeg = BLR.cloudModel(32; rng=MersenneTwister(7), cloudArgs...)
    @test mLeg.params.constructor == :cloudModel
    @test mLeg.params.rng == :legacy

    # vectors-form cloudModel records :cloudModelVectors; rebuild reproduces the (input) ϕ₀
    ϕ₀v = collect(range(0.0, 2π, length=9))[1:8]
    mVec = BLR.cloudModel(ϕ₀v, fill(0.4, 8), fill(0.1, 8), fill(0.4, 8), 0.4, 1.0;
        μ=600.0, I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, rng=MersenneTwister(3), τ=0.0)
    @test mVec.params.constructor == :cloudModelVectors
    @test isequal(BLR.getVariable(BLR.rebuild(mVec), :ϕ₀, flatten=true),
                  BLR.getVariable(mVec, :ϕ₀, flatten=true))

    # (d) combined diskwind + philox-cloud records the nested :+ node; rebuild reproduces both
    # submodels bit-identically INCLUDING subModelStartInds (use isequal -- disk has NaN sentinels)
    mCombined = mDW + mC
    @test mCombined.params.constructor == :+
    @test mCombined.params.left.constructor == :DiskWindModel
    @test mCombined.params.right.constructor == :cloudModel
    mCombinedR = BLR.rebuild(mCombined)
    @test mCombinedR.subModelStartInds == mCombined.subModelStartInds
    @test isequal(BLR.getVariable(mCombinedR, :v, flatten=true),
                  BLR.getVariable(mCombined, :v, flatten=true))
    @test isequal(BLR.getVariable(mCombinedR, :I, flatten=true),
                  BLR.getVariable(mCombined, :I, flatten=true))
    @test BLR.getProfile(mCombinedR, :line, bins=21, centered=true).binSums ==
          BLR.getProfile(mCombined, :line, bins=21, centered=true).binSums

    # (d′) right-nested a+(b+e) (3 slots post-W4-T0) records the nested tree and rebuilds with
    # identical subModelStartInds and getFlattenedCameraIndices
    a = BLR.DiskWindModel(300.0, 500.0, 0.4; nr=6, nϕ=12, scale=:linear,
        I=BLR.IsotropicIntensity, v=BLR.vCircularDisk, τ=0.4)
    b = BLR.cloudModel(16; rng=:philox, seed=201, cloudArgs...)
    e = BLR.DiskWindModel(600.0, 800.0, 0.4; nr=6, nϕ=12, scale=:linear,
        I=BLR.IsotropicIntensity, v=BLR.vCircularDisk, τ=0.4)
    rn = a + (b + e)
    @test rn.params.constructor == :+
    @test rn.params.right.constructor == :+ #nested, not flattened
    @test length(rn.subModelStartInds) == 3
    rnR = BLR.rebuild(rn)
    @test rnR.subModelStartInds == rn.subModelStartInds
    @test BLR.getFlattenedCameraIndices(rnR) == BLR.getFlattenedCameraIndices(rn)

    # (e) raytrace! of the (d) model with a non-uniform IRatios vector records the :raytrace! wrapper
    # and rebuild reproduces the raytraced :line profile bit-identically (raytrace! RETURNS a new model)
    rt = BLR.raytrace!(mDW + mC; IRatios=[1.0, 2.0])
    @test rt.params.constructor == :raytrace!
    @test rt.params.IRatios == [1.0, 2.0]
    @test rt.params.parent.constructor == :+
    rtR = BLR.rebuild(rt)
    @test BLR.getProfile(rtR, :line, bins=21, centered=true).binSums ==
          BLR.getProfile(rt, :line, bins=21, centered=true).binSums

    # (f) an unmatched override key throws
    @test_throws ErrorException BLR.rebuild(mDW; notARealKey=1.0)
    # a valid override is accepted and broadcasts to the leaf
    @test typeof(BLR.rebuild(mDW; τ=0.9)) == BLR.model

    # (g) getindex propagates the matching :+ params subtree to extracted submodels,
    # for both association orders (a, b, e reused from (d′) above)
    ln = (a + b) + e
    for m3 in (ln, rn)
        @test m3[1].params == a.params
        @test m3[2].params == b.params
        @test m3[3].params == e.params
    end
    # the extracted record is rebuild-able (philox cloud -> bit-identical gas)
    @test _cloud_field_tuple(BLR.rebuild(rn[2])) == _cloud_field_tuple(b)
    # single-slot extraction of a leaf model is the whole model -> the record transfers
    @test mDW[1].params == mDW.params
    # a side built without a public constructor stays params-less without breaking the other side
    # (the known side's slot count pins the unknown side's arithmetically)
    aRaw = deepcopy(a); aRaw.params = nothing
    mMix = aRaw + e
    @test mMix[1].params === nothing
    @test mMix[2].params == e.params
    # submodels sliced out of a multi-slot raytraced model carry no construction record (raytrace!
    # compacts/regroups slots, and re-raytracing a slice alone is not equivalent to slicing the
    # raytraced whole -- occlusion couples the submodels); rt from (e) above is multi-slot here
    @test length(rt.subModelStartInds) > 1
    @test all(rt[j].params === nothing for j in 1:length(rt.subModelStartInds))
    # ...but a leaf added on top of a raytraced model still resolves
    mRtPlus = rt + e
    nSlots = length(mRtPlus.subModelStartInds)
    @test mRtPlus[nSlots].params == e.params
    @test mRtPlus[1].params === nothing
end

@testset "CompositeModel + addLine! (W4-T2/T3)" begin
    # Models here are deliberately tiny -- this testset should add seconds, not minutes.
    cloudArgs = (; μ=600.0, β=1.0, F=0.5, θₒ=0.4, i=0.4, γ=1.0, ξ=1.0,
        I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=0.1)

    # --- T2: construction, duplicate-name error, indexing error, show ---
    mHa = BLR.DiskWindModel(300.0, 900.0, 0.4; nr=6, nϕ=12, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.0, f3=0.0, f4=0.0, α=1.0,
        τ=0.4, reflect=false)
    cm = BLR.CompositeModel(mHa; line="Ha", lineCenter=6562.8)
    @test cm.lines == ["Ha"]
    @test cm.models["Ha"] === mHa
    @test cm.lineCenters["Ha"] == 6562.8
    @test cm.fluxRatios["Ha"] == 1.0
    @test length(cm) == 1
    @test_throws ErrorException BLR.CompositeModel(mHa; line="Ha", lineCenter=-1.0) #lineCenter must be positive

    # duplicate name error
    @test_throws ErrorException BLR.addLine!(cm, mHa; line="Ha", lineCenter=1.0)

    # indexing error lists known lines
    err = try; cm["nope"]; nothing; catch e; e; end
    @test err isa ErrorException
    @test occursin("nope", err.msg) && occursin("Ha", err.msg)

    # show output smoke test
    s = sprint(show, cm)
    @test occursin("Ha", s)
    @test occursin("CompositeModel", s)

    # --- T3(a): rMin/rMax-form DiskWindModel reuse with an intensity-only override (α) ---
    # identical r/ϕ/v, different I (DiskWindIntensity's α is a source-function power law that
    # reaches only the intensity, since rMin/rMax are given explicitly in this parameterization)
    BLR.addLine!(cm; line="Hb", lineCenter=4861.3, from="Ha", α=2.5)
    mHb = cm["Hb"]
    @test isequal(BLR.getVariable(mHa, :r, flatten=true), BLR.getVariable(mHb, :r, flatten=true))
    @test isequal(BLR.getVariable(mHa, :ϕ, flatten=true), BLR.getVariable(mHb, :ϕ, flatten=true))
    @test isequal(BLR.getVariable(mHa, :v, flatten=true), BLR.getVariable(mHb, :v, flatten=true))
    @test !isequal(BLR.getVariable(mHa, :I, flatten=true), BLR.getVariable(mHb, :I, flatten=true))

    # --- T3(a'): the SAME kind of α override on a mean-form (:DiskWindModelMean) model
    # MOVES the grid -- DIFFERENT r arrays. Pins the honest (non-invariant) constructor
    # semantics -- do not "fix" this.
    mMeanA = BLR.DiskWindModel(500.0, 5.0, 1.0, 0.4; nr=6, nϕ=12, scale=:log,
        I=BLR.IsotropicIntensity, v=BLR.vCircularDisk, τ=0.4)
    cmMean = BLR.CompositeModel(mMeanA; line="A", lineCenter=1000.0)
    BLR.addLine!(cmMean; line="B", lineCenter=2000.0, α=2.5)
    @test !isequal(BLR.getVariable(cmMean["A"], :r, flatten=true), BLR.getVariable(cmMean["B"], :r, flatten=true))

    # --- T3(b): philox cloud reuse without seed override -> identical r; new seed -> different ---
    mC1 = BLR.cloudModel(64; rng=:philox, seed=42, cloudArgs...)
    cmC = BLR.CompositeModel(mC1; line="C1", lineCenter=1.0)
    BLR.addLine!(cmC; line="C2", lineCenter=2.0) #no seed override -> reuse stored seed
    @test isequal(BLR.getVariable(cmC["C1"], :r, flatten=true), BLR.getVariable(cmC["C2"], :r, flatten=true))
    BLR.addLine!(cmC; line="C3", lineCenter=3.0, seed=99) #new seed -> different clouds
    @test !isequal(BLR.getVariable(cmC["C1"], :r, flatten=true), BLR.getVariable(cmC["C3"], :r, flatten=true))
    # seed=nothing override is rejected (philox requires an explicit seed)
    @test_throws ArgumentError BLR.addLine!(cmC; line="C4", lineCenter=4.0, seed=nothing)

    # legacy-rng reuse warns once and rebuilds with GLOBAL_RNG (statistically different, not identical)
    mLeg = BLR.cloudModel(16; rng=MersenneTwister(7), cloudArgs...)
    cmLeg = BLR.CompositeModel(mLeg; line="L1", lineCenter=1.0)
    @test_logs (:warn,) BLR.addLine!(cmLeg; line="L2", lineCenter=2.0)

    # --- T3(c): reuse from a params-less model throws the documented error ---
    mRaw = BLR.model(mHa.rings[1:2], nothing, nothing, [1]) #low-level constructor: params === nothing
    cmRaw = BLR.CompositeModel(mRaw; line="R", lineCenter=1.0)
    @test_throws ErrorException BLR.addLine!(cmRaw; line="R2", lineCenter=2.0)
    # sanity: the explicit-model method doesn't care about params at all
    BLR.addLine!(cmRaw, mHa; line="R3", lineCenter=3.0)
    @test cmRaw.lines == ["R", "R3"]
end

@testset "composite forwarding + spectrum (W4-T4/T5/T6)" begin
    # Models here are deliberately tiny -- this testset should add seconds, not minutes.
    # rMin/rMax comfortably above 200 rₛ so |v| = √(rₛ/2r) stays well below 0.05c (T5 needs this).
    mDisk = BLR.DiskWindModel(300.0, 900.0, 0.4; nr=8, nϕ=16, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.0, f3=0.0, f4=0.0, α=1.0,
        τ=0.4, reflect=false)

    # a two-line composite sharing the geometry: Hα (ref, fluxRatio 1.0) + Hβ (fluxRatio 0.35)
    cm = BLR.CompositeModel(mDisk; line="Ha", lineCenter=6562.8)
    BLR.addLine!(cm, mDisk; line="Hb", lineCenter=4861.3, fluxRatio=0.35)

    # --- T4: getProfile forwarding equals the direct call FIELD-WISE (profile has no ==) ---
    for line in cm.lines
        p = BLR.getProfile(cm, :line; line=line, bins=21, centered=true)
        q = BLR.getProfile(cm[line], :line; bins=21, centered=true)
        @test p.binSums == q.binSums
        @test p.binCenters == q.binCenters
        @test p.binEdges == q.binEdges
    end

    # --- T4: getVariable forwarding ---
    @test isequal(BLR.getVariable(cm, :v; line="Ha", flatten=true),
                  BLR.getVariable(cm["Ha"], :v, flatten=true))
    @test isequal(BLR.getVariable(cm, :I; line="Hb", flatten=false),
                  BLR.getVariable(cm["Hb"], :I, flatten=false))

    # --- T4: :ratio is a W5 placeholder; a missing `line` kwarg errors naming the kwarg ---
    @test_throws ErrorException BLR.getProfile(cm, :ratio; lines=("Ha", "Hb"))
    errNoLine = try; BLR.getProfile(cm, :line); nothing; catch e; e; end
    @test errNoLine isa ErrorException && occursin("line", errNoLine.msg)

    # --- T4: _fluxWeights identity -- default (centered=true) :line profile integrates to fluxRatio ---
    w = BLR._fluxWeights(cm)
    for line in cm.lines
        lp = BLR.getProfile(cm, :line; line=line, centered=true) #default centered edges pad past data
        @test isapprox(sum(lp.binSums) * w[line], cm.fluxRatios[line], rtol=1e-12)
    end

    # --- T4: raytrace! forwarding REASSIGNS the returned model; single-submodel lines are skipped ---
    mCl = BLR.cloudModel(20; rng=:philox, seed=7, μ=600.0, β=1.0, F=0.5, θₒ=0.4, i=0.4,
        γ=1.0, ξ=1.0, I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=0.1)
    mComb = mDisk + mCl #two submodels -> raytrace! does real work
    cmRT = BLR.CompositeModel(mComb; line="X", lineCenter=6562.8)
    BLR.addLine!(cmRT, mDisk; line="Y", lineCenter=4861.3) #single submodel -> skipped by raytrace!
    @test cmRT["X"].camera.raytraced == false
    origX = cmRT["X"]
    BLR.raytrace!(cmRT)
    @test cmRT["X"].camera.raytraced == true    #flipped
    @test cmRT["X"] !== origX                    #reassigned to the NEW returned model
    @test cmRT["Y"].camera.raytraced == false    #single-submodel line left unaltered (no warn-spam)

    # --- T5: wavelength() scalar/vector round-trip ---
    @test BLR.wavelength(0.0, 6562.8) == 6562.8
    vvec = [-0.01, 0.0, 0.02]
    λvec = BLR.wavelength(vvec, 6562.8)
    @test λvec == 6562.8 .* (1.0 .+ vvec)
    @test λvec ./ 6562.8 .- 1.0 ≈ vvec #recover v

    # --- T5: overlap check. Hα/Hβ (widely separated centers) do NOT overlap ---
    @test isempty(BLR.lineOverlap(cm))

    # push centers together (same geometry) until the intervals DO overlap; check the reported interval
    vmin, vmax = BLR._finiteVRange(mDisk)
    @test vmin < 0.0 && vmax > 0.0 && abs(vmin) < 0.05 && abs(vmax) < 0.05
    cA = 6563.0; cB = 6600.0 #cB/cA = 1.0056 < (1+vmax)/(1+vmin), so they overlap
    cmO = BLR.CompositeModel(mDisk; line="A", lineCenter=cA)
    BLR.addLine!(cmO, mDisk; line="B", lineCenter=cB)
    ov = BLR.lineOverlap(cmO)
    @test length(ov) == 1
    @test ov[1].lineA == "A" && ov[1].lineB == "B"
    @test isapprox(ov[1].λlo, cB*(1.0+vmin), rtol=1e-12) #hand computation: max(cA,cB)*(1+vmin)=cB*(1+vmin)
    @test isapprox(ov[1].λhi, cA*(1.0+vmax), rtol=1e-12) #min(cA,cB)*(1+vmax)=cA*(1+vmax)

    # --- T6: combined spectrum ---
    edges, centers, flux, total = BLR.getSpectrum(cm; bins=64)
    @test length(centers) == 64 && length(edges) == 65
    # integrated flux per line equals its fluxRatio (overflow keeps the boundary points)
    @test isapprox(sum(flux["Ha"]), 1.0, rtol=1e-12)
    @test isapprox(sum(flux["Hb"]), 0.35, rtol=1e-12)
    # total is the elementwise sum of the parts
    @test total == flux["Ha"] .+ flux["Hb"]
    # two non-overlapping lines -> two disjoint nonzero regions with a zero gap between them
    haNZ = findall(>(0.0), flux["Ha"]); hbNZ = findall(>(0.0), flux["Hb"])
    @test !isempty(haNZ) && !isempty(hbNZ)
    @test isempty(intersect(haNZ, hbNZ))
    @test any(==(0.0), total) #a gap of empty bins between the two lines

    # redshift shifts the grid by (1+z); integrated fluxes are unchanged
    e2, c2, flux2, total2 = BLR.getSpectrum(cm; bins=64, z=0.1)
    @test isapprox(sum(flux2["Ha"]), 1.0, rtol=1e-12)
    @test isapprox(edges[1]*1.1, e2[1], rtol=1e-12) && isapprox(edges[end]*1.1, e2[end], rtol=1e-12)

    # a zero-flux line errors loudly by name instead of silently vanishing from the spectrum
    # (adversarial review 2026-07-03: fluxRatio/0 = Inf weight would otherwise make binAccumulate!
    # silently drop every weighted point of that line, breaking sum(flux[line]) == fluxRatio)
    mZero = BLR.DiskWindModel(300.0, 900.0, 0.4; nr=6, nϕ=12, scale=:linear,
        I=BLR.IsotropicIntensity, rescale=0.0, v=BLR.vCircularDisk, τ=0.4, reflect=false)
    cmZero = BLR.CompositeModel(mDisk; line="Ha", lineCenter=6563.0)
    BLR.addLine!(cmZero, mZero; line="dark", lineCenter=4861.0)
    errZ = try BLR._fluxWeights(cmZero); nothing; catch err; err; end
    @test errZ isa ErrorException && occursin("dark", errZ.msg)
    @test_throws ErrorException BLR.getSpectrum(cmZero; bins=32)
end

@testset "composite multiline models (W4-T8)" begin
    # End-to-end scenario from the W4 plan. Models are deliberately tiny -- seconds, not minutes.
    # rMin/rMax-form (explicit) DiskWindModel as Hα (reference line; in this parameterization α is a
    # source-function power law that reaches ONLY the intensity function -- see the W4-T3 testset) +
    # Hβ added via addLine! parameter reuse with a different α and fluxRatio=0.35.
    # Radii matter for (c): |v| = √(rₛ/2r) ≲ 0.05 needs r ≳ 200 rₛ, so rMin=300 keeps the Hα/Hβ
    # wavelength intervals physically disjoint (at smaller radii the overlap would be REAL).
    rMin, rMax, inc = 300.0, 900.0, 0.4
    # f2-only (Keplerian shear) intensity gives a double-horned profile whose peaks sit at nonzero
    # ±|v| -- a meaningful peak position to pin, in the style of the existing DiskWind peak tests.
    mHα = BLR.DiskWindModel(rMin, rMax, inc; nr=6, nϕ=64, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=0.0, f2=1.0, f3=0.0, f4=0.0, α=1.0,
        τ=0.4, reflect=false)
    cm = BLR.CompositeModel(mHα; line="Hα", lineCenter=6563.0)
    BLR.addLine!(cm; line="Hβ", lineCenter=4861.0, fluxRatio=0.35, α=2.5) #param reuse, intensity-only override

    # (a) identical geometry arrays between the two lines (α override moved only I, not the grid)
    for var in (:r, :ϕ, :v)
        @test isequal(BLR.getVariable(cm, var; line="Hα", flatten=true),
                      BLR.getVariable(cm, var; line="Hβ", flatten=true))
    end

    # (b) getSpectrum integrated fluxes = the fluxRatios (integrated-flux semantics, decision 3)
    edges, centers, flux, total = BLR.getSpectrum(cm; bins=64)
    @test isapprox(sum(flux["Hα"]), 1.0, rtol=1e-12)
    @test isapprox(sum(flux["Hβ"]), 0.35, rtol=1e-12)

    # (c) Hα/Hβ do not overlap in wavelength at these radii
    @test isempty(BLR.lineOverlap(cm))

    # (d) each line's profile peaks at the same |v| positions as the equivalent standalone model
    # (standalones constructed independently -- disk-wind construction is deterministic)
    mHαAlone = BLR.DiskWindModel(rMin, rMax, inc; nr=6, nϕ=64, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=0.0, f2=1.0, f3=0.0, f4=0.0, α=1.0,
        τ=0.4, reflect=false)
    mHβAlone = BLR.DiskWindModel(rMin, rMax, inc; nr=6, nϕ=64, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=0.0, f2=1.0, f3=0.0, f4=0.0, α=2.5,
        τ=0.4, reflect=false)
    for (line, standalone) in (("Hα", mHαAlone), ("Hβ", mHβAlone))
        p = BLR.getProfile(cm, :line; line=line, bins=41, centered=true)
        q = BLR.getProfile(standalone, :line; bins=41, centered=true)
        pPeak = p.binCenters[findmax(p.binSums)[2]]
        qPeak = q.binCenters[findmax(q.binSums)[2]]
        @test abs(pPeak) == abs(qPeak)
        @test abs(pPeak) > 0.0 #double-horned: the peak is NOT at line center
        # second horn mirrors the first (existing peak-pinning style)
        firstMax = findmax(p.binSums)[2]
        mask = [j != firstMax for j in eachindex(p.binSums)]
        secondPeak = p.binCenters[mask][findmax(p.binSums[mask])[2]]
        @test isapprox(secondPeak, -pPeak, atol=1e-12)
    end

    # (e) SIGN-CONVENTION PIN (locked decision 4: stored v is redshift-positive -- near-side inflow
    # stores +v). A FULL disk with pure radial inflow is exactly front/back symmetric in v (the far
    # side's blueshift mirrors the near side's redshift), so mask the emission to the NEAR side
    # (ϕ ∈ (-π/2, π/2), the half tilted toward the camera -- IϕDiskWindMask's ϕ=0 is at the camera)
    # with isotropic (f4-only) emission: pure near-side inflow must put essentially ALL the flux
    # REDWARD of lineCenter. The symmetric peak tests above cannot catch a mirrored spectrum; this
    # asymmetric case can. If this test fails, the velocity sign convention broke -- do NOT "fix"
    # signs in src/ to make it pass (the convention was verified against the source 2026-07-01).
    mIn = BLR.DiskWindModel(300.0, 900.0, 0.4; nr=6, nϕ=64, scale=:linear, τ=0.1, reflect=false,
        I=BLR.IϕDiskWindMask, v=BLR.vCircularRadialDisk, vᵣFrac=1.0, inflow=true,
        f1=0.0, f2=0.0, f3=0.0, f4=1.0, α=0.0, ϕMin=-π/2+0.02, ϕMax=π/2-0.02)
    cmIn = BLR.CompositeModel(mIn; line="inflow", lineCenter=6563.0)
    _, cIn, fluxIn, _ = BLR.getSpectrum(cmIn; bins=101)
    lc = cmIn.lineCenters["inflow"]
    redFlux = sum(fluxIn["inflow"][cIn .> lc])
    blueFlux = sum(fluxIn["inflow"][cIn .< lc])
    @test redFlux > blueFlux
    @test redFlux > 0.99 && blueFlux < 0.01 #near-side pure inflow: all redshifted (unit integral)
end

include("raytrace_reference.jl")
include("gpu_arrays.jl")
include("gpu_kernels.jl")
include("recipes.jl")

# GPU correctness tests are opt-in: they need a CUDA device and the CUDA.jl weak dependency in the
# active environment. Run with `BLR_TEST_CUDA=1` from an env that has CUDA.jl available.
if get(ENV, "BLR_TEST_CUDA", "0") == "1"
    include("gpu_cuda.jl")
else
    @info "skipping GPU correctness tests (set BLR_TEST_CUDA=1 with CUDA.jl available to run them)"
end

## NOTE add JET to the test environment, then uncomment
# using JET
# @testset "static analysis with JET.jl" begin
#     JET.test_package(BLR, target_modules=(BLR,))
# end

## NOTE add Aqua to the test environment, then uncomment
# @testset "QA with Aqua" begin
#     import Aqua
#     Aqua.test_all(BLR)
# end
