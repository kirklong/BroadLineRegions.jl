using StableRNGs: StableRNG

@testset "raytrace combiner" begin
    disk(; r1=300.0, r2=900.0, nr=8, nϕ=16, inc=0.4, τ=0.4) =
        BLR.DiskWindModel(r1, r2, inc, nr=nr, nϕ=nϕ, scale=:linear,
            I=BLR.IsotropicIntensity, v=BLR.vCircularDisk, τ=τ, reflect=false)

    clouds(n, seed; μ=600.0, inc=0.4, τ=0.1) =
        BLR.cloudModel(n, μ=μ, β=1.0, F=0.5, θₒ=0.4, i=inc, γ=1.0, ξ=1.0,
            I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=τ, rng=StableRNG(seed))

    raytrace_summary(m; kwargs...) = begin
        rt = BLR.raytrace!(m; kwargs...)
        I = BLR.getVariable(rt, :I, flatten=true)
        ΔA = BLR.getVariable(rt, :ΔA, flatten=true)
        LP = BLR.getProfile(rt, :line, bins=11, centered=true)
        (rt=rt, I=I, ΔA=ΔA, flux=sum(I .* ΔA), line=LP.binSums)
    end

    # Continuous-model baselines derive from raytrace! at commit a96fba6, with
    # partial obscuration accounted for. StableRNG keeps cloud fixtures identical
    # across Julia versions. Profile tolerances allow floating-point reassociation.
    sDC = raytrace_summary(disk() + clouds(80, 101), τCutOff=1.0)
    @test !any(isnan, sDC.I)
    @test isapprox(sDC.flux, 2.0579449066647638e6, rtol=1e-9)
    # The cloud-profile baseline uses the CCW orbital tangent in
    # vCircularCloud/vCloudTurbulentEllipticalFlow.
    @test all(isapprox.(sDC.line,
        [12788.5143410726, 114441.91252453177, 332069.3055598979,
         189578.46476604807, 246954.69609633373, 266290.68656407285,
         151530.12725557524, 202676.7106953462, 306816.39818805375,
         180537.36862948653, 54260.72204434513], rtol=1e-8))

    # Two overlapping continuous disks with uniform weights.
    ddOverlap = raytrace_summary(disk(r1=250.0, r2=700.0, nr=6, nϕ=12) +
                                 disk(r1=500.0, r2=1000.0, nr=6, nϕ=12);
                                 IRatios=[1.0, 1.0], τCutOff=1.0)
    @test !any(isnan, ddOverlap.I)
    @test isapprox(ddOverlap.flux, 3.1808156444130e6, rtol=1e-9)
    @test all(isapprox.(ddOverlap.line,
        [15948.1628068452, 170891.15213961044, 632870.6737718652,
         304022.78831852623, 189400.88227433746, 554548.3257905911,
         189400.88227433746, 304022.78831852623, 632870.6737718652,
         170891.15213961044, 15948.1628068452], rtol=1e-7))

    # Mixed disk/cloud results are independent of submodel order.
    sCD = raytrace_summary(clouds(80, 101) + disk(), τCutOff=1.0)
    @test all(isapprox.(sCD.line, sDC.line, rtol=1e-12))
    @test isapprox(sCD.flux, sDC.flux, rtol=1e-12)

    # One disk and two independent cloud submodels exercise multiple contributors
    # per sightline. Flux values here and below are combiner regression baselines.
    s3 = raytrace_summary(disk() + clouds(60, 102) + clouds(40, 103, μ=850.0), τCutOff=1.0)
    @test !any(isnan, s3.I)
    @test isfinite(s3.flux)
    @test isapprox(s3.flux, 2.0579620075784149e6, rtol=1e-8)

    # Three overlapping continuous disks form a chain; overlaps are resolved
    # iteratively while building the output grid innermost-first.
    s3disk = raytrace_summary(disk(r1=250.0, r2=700.0, nr=6, nϕ=12) +
                              disk(r1=500.0, r2=1000.0, nr=6, nϕ=12) +
                              disk(r1=800.0, r2=1400.0, nr=6, nϕ=12);
                              IRatios=[1.0, 1.0, 1.0], τCutOff=1.0)
    @test !any(isnan, s3disk.I)
    @test length(s3disk.I) == length(s3disk.rt.camera.α)
    @test isapprox(s3disk.flux, 6.7772238771869e6, rtol=1e-8)

    # IRatios is a global per-submodel weight, independent of submodel order.
    # Disk+cloud and cloud+disk with matched weights describe the same physics.
    sDC2 = raytrace_summary(disk() + clouds(80, 101), IRatios=[1.0, 0.25], τCutOff=1.0)
    sCD2 = raytrace_summary(clouds(80, 101) + disk(), IRatios=[0.25, 1.0], τCutOff=1.0)
    @test all(isapprox.(sDC2.line, sCD2.line, rtol=1e-12))
    @test isapprox(sDC2.flux, sCD2.flux, rtol=1e-12)

    # Flux conservation under area resampling: a cloud merged into a disk
    # pixel contributes the same flux I*ΔA it would alone. Adding one cloud of
    # weighted intensity 0.25 raises the total I*ΔA by exactly 0.25.
    dFluxModel = disk()
    diskFlux = sum(BLR.getVariable(dFluxModel, :I, flatten=true) .*
                   BLR.getVariable(dFluxModel, :ΔA, flatten=true))
    # Stable seed 3 places the cloud in the unobscured merged-pixel case.
    diskCloud = BLR.raytrace!(disk() + clouds(1, 3, μ=600.0), IRatios=[1.0, 0.25], τCutOff=1.0)
    diskCloudFlux = sum(BLR.getVariable(diskCloud, :I, flatten=true) .*
                        BLR.getVariable(diskCloud, :ΔA, flatten=true))
    @test isapprox(diskCloudFlux - diskFlux, 0.25, rtol=1e-8, atol=1e-8)

    # Free-cloud raytracing attenuates a rear cloud without merging it: it is
    # dimmed by exp(-τ) of the cloud in front but keeps its own velocity
    # (no merging into an intensity-weighted-mean point).
    front = BLR.ring(r=10.0, i=0.0, rot=0.0, θₒ=0.0, v=1.0, I=1.0,
        ϕ=0.0, ϕ₀=0.0, ΔA=Float64(pi), reflect=false, τ=0.6, η=1.0,
        Δr=1.0, Δϕ=1.0, scale=nothing, x=5.0, y=0.0, z=0.0)
    back = BLR.ring(r=10.0, i=0.0, rot=0.0, θₒ=0.0, v=-1.0, I=2.0,
        ϕ=0.0, ϕ₀=0.0, ΔA=Float64(pi), reflect=false, τ=0.1, η=1.0,
        Δr=1.0, Δϕ=1.0, scale=nothing, x=0.0, y=0.0, z=0.0)
    free = BLR.raytrace!(BLR.model([front]) + BLR.model([back]), raytraceFreeClouds=true, τCutOff=1.0)
    vFree = BLR.getVariable(free, :v, flatten=true)
    IFree = BLR.getVariable(free, :I, flatten=true)
    @test sort(vFree) == [-1.0, 1.0]
    @test isapprox(IFree[findfirst(==(-1.0), vFree)], 2.0 * exp(-0.6), rtol=1e-12)
    @test isapprox(IFree[findfirst(==(1.0), vFree)], 1.0, rtol=1e-12)

    # Partial obscuration (obscuredFrac > 1): a small front tile (area 1, τ=0.5)
    # sits in front of a larger rear tile (area 4) in one sightline. The rear tile
    # is 4× bigger, so 1/4 of it is
    # covered (attenuated by exp(-0.5)) and 3/4 peeks around at full strength:
    #   I = I_front + (exp(-0.5)·1 + 3) = 1 + exp(-0.5) + 3 = 4 + exp(-0.5).
    pf(x, ΔA, τ) = BLR._RaytracePoint(1, 1, 1, 0.0, 0.0, 0.0, 0.0, x, 1.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, ΔA, τ, 1.0, false, false, Inf)
    pts = [pf(10.0, 1.0, 0.5), pf(0.0, 4.0, 0.0)]   # front (small), rear (large)
    res = BLR._rt_scan_bucket(pts, [1, 2], [1.0], 1.0, 1.0)   # outputΔA=1, τCutOff=1
    @test isapprox(res.I, 4.0 + exp(-0.5), rtol=1e-12)

    # Disk-obscured-cloud filters classify rings by type across multiple cloud groups.
    obsc3() = disk(nr=6, nϕ=12) + clouds(30, 11) + clouds(30, 12, μ=750.0)
    z3 = BLR.zeroDiskObscuredClouds!(obsc3())
    @test !any(isnan, BLR.getVariable(z3, :I, flatten=true))
    r3 = BLR.removeDiskObscuredClouds!(obsc3())
    Ir3 = BLR.getVariable(r3, :I, flatten=true)
    @test !any(isnan, Ir3)
    @test length(Ir3) <= length(BLR.getVariable(obsc3(), :I, flatten=true))  # some clouds removed
end
