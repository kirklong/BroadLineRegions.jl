using BLR
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
    @test minimum(LPAll.binSums) > 0.0
    @test minimum(LPf1.binSums) > 0.0
    @test minimum(LPf2.binSums) > 0.0
    @test minimum(LPf3.binSums) > 0.0
    @test minimum(LPf4.binSums) > 0.0
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
    @test isapprox(tCenters[findmax(Ψt)[2]]*rsDay, 10.0, atol = 5e-1)
end

@testset "Cloud model initialized successfully" begin 
    mP1 = BLR.cloudModel(1_000_000,μ=4/rsDay,F=0.25,β=1.0,θₒ=40/180*π,i=20/180*π,
               κ=-0.4,γ=5.0,ξ=0.3,fEllipse=0.0,fFlow=0.0,θₑ=0.0,σₜ=0.0,
               σρc=0.0,σΘᵣ=0.0,σΘc=0.0,σρᵣ=0.0,
               I=BLR.cloudIntensity,v=BLR.vCloudTurbulentEllipticalFlow,τ=0.0)
    @test typeof(mP1) == BLR.model
    LP1 = BLR.getProfile(mP1,:line,bins=101,centered=true)
    @test maximum(LP1.binSums) > 0.0
    @test minimum(LP1.binSums) > 0.0
    @test isapprox(LP1.binCenters[findmax(LP1.binSums)[2]],0.025,atol = 5e-3)
    mP2 = BLR.cloudModel(1_000_000,μ=4/rsDay,F=0.25,β=0.8,θₒ=30/180*π,i=20/180*π,
        κ=-0.4,γ=5.0,ξ=0.1,fEllipse=1.0,σₜ=0.0,
        fFlow=0.0,θₑ=0.0,σρc=0.0,σΘᵣ=0.0,σΘc=0.0,σρᵣ=0.0,
        I=BLR.cloudIntensity,v=BLR.vCloudTurbulentEllipticalFlow,τ=0.0)
    LP2 = BLR.getProfile(mP2,:line,bins=101,centered=true)
    @test maximum(LP2.binSums) > 0.0
    @test minimum(LP2.binSums) > 0.0
    firstMax = findmax(LP2.binSums)[2]
    @test isapprox(LP2.binCenters[firstMax], (0.007*sign(LP2.binCenters[firstMax])), atol = 5e-3)
    mask = sign(LP2.binCenters[firstMax]) == 1 ? [i < 51 for i=1:length(LP2.binCenters)] :
        [i > 50 for i=1:length(LP2.binCenters)]
    secondMax = findmax(LP2.binSums[mask])[2]
    @test isapprox(LP2.binCenters[mask][secondMax],(0.007*sign(LP2.binCenters[mask][secondMax])),atol = 5e-3)
    M = 10^(6.5)*2e30; rs = 2*6.67e-11*M/9e16; rsDay = rs/3e8/24/3600
    tCenters,Ψt = BLR.getΨt(mP1,501,0.5/rsDay)
    @test isapprox(tCenters[findmax(Ψt)[2]]*rsDay, 0.025, atol = 5e-3)
    tCenters,Ψt = BLR.getΨt(mP2,501,0.5/rsDay)
    @test isapprox(tCenters[findmax(Ψt)[2]]*rsDay, 0.03, atol = 5e-2)
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
