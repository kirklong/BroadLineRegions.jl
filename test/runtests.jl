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
        mask = (dRef .>= tEdges[j]) .& (dRef .< tEdges[j+1])
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
            mask = (v .>= vEdges[i]) .& (v .< vEdges[i+1]) .& (d .>= tEdges[j]) .& (d .< tEdges[j+1])
            s = sum(I[mask].*ΔA[mask])
            expected = s > 0 ? s : 1e-30
            @test isapprox(Ψ[i,j], expected, rtol=1e-12)
            @test (Ψ[i,j] == 1e-30) == (expected == 1e-30) #empty/poisoned bins identical
        end
        for ovf in (false, true)
            Ψt = BLR.getΨt(m, tEdges, ovf)
            ref = map(1:length(tEdges)-1) do j
                mask = (d .>= tEdges[j]) .& (d .< tEdges[j+1])
                s = sum(I[mask].*ΔA[mask])
                s > 0 ? s : 1e-30
            end
            if ovf
                sU = sum(I[(d .< tEdges[1])].*ΔA[(d .< tEdges[1])])
                sO = sum(I[(d .>= tEdges[end])].*ΔA[(d .>= tEdges[end])])
                ref[1] += sU > 0 ? sU : 1e-30
                ref[end] += sO > 0 ? sO : 1e-30
            end
            @test all(isapprox.(Ψt, ref, rtol=1e-12))
        end
    end
end

include("raytrace_reference.jl")
include("gpu_arrays.jl")

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
