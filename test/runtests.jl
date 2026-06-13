using BroadLineRegions
using Test

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
