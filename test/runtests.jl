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
