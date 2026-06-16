@testset "raytrace Phase A combiner" begin
    disk(; r1=300.0, r2=900.0, nr=8, nϕ=16, inc=0.4, τ=0.4) =
        BLR.DiskWindModel(r1, r2, inc, nr=nr, nϕ=nϕ, scale=:linear,
            I=BLR.IsotropicIntensity, v=BLR.vCircularDisk, τ=τ, reflect=false)

    clouds(n, seed; μ=600.0, inc=0.4, τ=0.1) =
        BLR.cloudModel(n, μ=μ, β=1.0, F=0.5, θₒ=0.4, i=inc, γ=1.0, ξ=1.0,
            I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=τ, rng=MersenneTwister(seed))

    raytrace_summary(m; kwargs...) = begin
        rt = BLR.raytrace!(m; kwargs...)
        I = BLR.getVariable(rt, :I, flatten=true)
        ΔA = BLR.getVariable(rt, :ΔA, flatten=true)
        LP = BLR.getProfile(rt, :line, bins=11, centered=true)
        (rt=rt, I=I, ΔA=ΔA, flux=sum(I .* ΔA), line=LP.binSums)
    end

    # ------------------------------------------------------------------
    # (1) EQUIVALENCE to the pre-rewrite raytrace!.
    # The flux/line values below were CAPTURED by running the original
    # raytrace! on origin/main (pre-rewrite, commit a96fba6) for the cases
    # that path could handle. The rewrite must reproduce them: this proves
    # the combiner did not change observable behaviour on the working paths.
    # (rtol 1e-8 absorbs floating-point re-association from the new order of
    # summation; old vs new agree to ~1e-13 in practice.)
    # ------------------------------------------------------------------
    # NOTE on provenance: the pre-rewrite raytrace! on origin/main gave
    #   flux = 2.0391087024702898e6
    #   line = [50951.8414, 228718.8705, 356006.2909, 246949.2266, 263973.4663,
    #           1.0, 303694.1658, 357323.0208, 196441.9867, 35047.8335, 1.0]
    # The rewrite reproduced those EXACTLY before the obscuredFrac fix. The values
    # below are origin/main + that one deliberate correction (partial-obscuration
    # formula; see the "obscuredFrac" unit test). Since the obscuredFrac edit is the
    # only delta from the verified-equivalent rewrite, these isolate that fix.
    sDC = raytrace_summary(disk() + clouds(80, 101), τCutOff=1.0)
    @test !any(isnan, sDC.I)
    @test isapprox(sDC.flux, 2.0579430970021e6, rtol=1e-9)
    @test all(isapprox.(sDC.line,
        [54265.913923740234, 234483.85358264975, 360127.72516325157,
         246949.22664271292, 266292.08547965763, 1.0,
         303694.77223813196, 357323.62727836345, 196442.59312292485,
         38361.2995706341, 1.0], rtol=1e-8))

    # Two overlapping continuous disks, uniform weights. origin/main handled a
    # single overlap pair through its main loop (only length(overlaps) > 1 hit
    # the broken recursive branch); the new unified combiner reproduces it.
    ddOverlap = raytrace_summary(disk(r1=250.0, r2=700.0, nr=6, nϕ=12) +
                                 disk(r1=500.0, r2=1000.0, nr=6, nϕ=12);
                                 IRatios=[1.0, 1.0], τCutOff=1.0)
    # origin/main pre-fix: flux=3.1843276387323737e6; the rewrite matched it exactly
    # before the obscuredFrac fix. Values below = origin/main + obscuredFrac fix.
    @test !any(isnan, ddOverlap.I)
    @test isapprox(ddOverlap.flux, 3.1808156444130e6, rtol=1e-9)
    @test all(isapprox.(ddOverlap.line,
        [15948.1628068452, 170891.15213961044, 632870.6737718652,
         304022.78831852623, 189400.88227433746, 554548.3257905911,
         189400.88227433746, 304022.78831852623, 632870.6737718652,
         170891.15213961044, 15948.1628068452], rtol=1e-7))

    # ------------------------------------------------------------------
    # (2) NEW CAPABILITIES — origin/main THROWS on these (so no golden
    # exists; values here are a forward regression baseline of the rewrite,
    # not an equivalence claim).
    # ------------------------------------------------------------------
    # Order independence. origin/main threw `setindex!(::Float64, ...)` when a
    # discrete (cloud) submodel came first; the rewrite is order-independent.
    sCD = raytrace_summary(clouds(80, 101) + disk(), τCutOff=1.0)
    @test all(isapprox.(sCD.line, sDC.line, rtol=1e-12))
    @test isapprox(sCD.flux, sDC.flux, rtol=1e-12)

    # One disk + TWO independent cloud submodels. origin/main threw the same
    # setindex! error (its main loop assumes scalar fields stay scalar). The
    # single scan core handles >2 contributors per sightline.
    s3 = raytrace_summary(disk() + clouds(60, 102) + clouds(40, 103, μ=850.0), τCutOff=1.0)
    @test !any(isnan, s3.I)
    @test isfinite(s3.flux)
    @test isapprox(s3.flux, 2.0579611782195e6, rtol=1e-8)  # rewrite baseline (origin/main throws on disk+2cloud)

    # THREE overlapping continuous disks (a chain). This is the case that hit the
    # old `length(overlaps) > 1` recursive branch, where origin/main threw
    # `non-boolean (Int64) used in boolean context` (the `endInd = i = length(...)`
    # `=` vs `==` bugs, raytrace.jl:674/687). The unified output-grid build
    # (innermost-first, overlaps resolved iteratively) replaces that recursion.
    s3disk = raytrace_summary(disk(r1=250.0, r2=700.0, nr=6, nϕ=12) +
                              disk(r1=500.0, r2=1000.0, nr=6, nϕ=12) +
                              disk(r1=800.0, r2=1400.0, nr=6, nϕ=12);
                              IRatios=[1.0, 1.0, 1.0], τCutOff=1.0)
    @test !any(isnan, s3disk.I)
    @test length(s3disk.I) == length(s3disk.rt.camera.α)
    @test isapprox(s3disk.flux, 6.7772238771869e6, rtol=1e-8)  # rewrite baseline (origin/main throws)

    # ------------------------------------------------------------------
    # (3) IRatios is a GLOBAL per-submodel weight (deliberate change from the
    # old "relative-to-base" scaling). It must be applied identically to a
    # component regardless of submodel order, so disk+cloud and cloud+disk
    # with matched weights are the same physics.
    # ------------------------------------------------------------------
    sDC2 = raytrace_summary(disk() + clouds(80, 101), IRatios=[1.0, 0.25], τCutOff=1.0)
    sCD2 = raytrace_summary(clouds(80, 101) + disk(), IRatios=[0.25, 1.0], τCutOff=1.0)
    @test all(isapprox.(sDC2.line, sCD2.line, rtol=1e-12))
    @test isapprox(sDC2.flux, sCD2.flux, rtol=1e-12)

    # ------------------------------------------------------------------
    # (4) Flux conservation under area resampling: a cloud merged into a disk
    # pixel contributes the same flux I*ΔA it would alone. Adding one cloud of
    # weighted intensity 0.25 raises the total I*ΔA by exactly 0.25.
    # ------------------------------------------------------------------
    dFluxModel = disk()
    diskFlux = sum(BLR.getVariable(dFluxModel, :I, flatten=true) .*
                   BLR.getVariable(dFluxModel, :ΔA, flatten=true))
    diskCloud = BLR.raytrace!(disk() + clouds(1, 501, μ=600.0), IRatios=[1.0, 0.25], τCutOff=1.0)
    diskCloudFlux = sum(BLR.getVariable(diskCloud, :I, flatten=true) .*
                        BLR.getVariable(diskCloud, :ΔA, flatten=true))
    @test isapprox(diskCloudFlux - diskFlux, 0.25, rtol=1e-8, atol=1e-8)

    # ------------------------------------------------------------------
    # (5) Free-cloud raytracing is attenuate-don't-merge: a rear cloud is
    # dimmed by exp(-τ) of the cloud in front but KEEPS its own velocity
    # (no merging into an intensity-weighted-mean point).
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # (6) Partial obscuration (obscuredFrac > 1), analytic check of the corrected
    # formula. A small front tile (area 1, τ=0.5) sits in front of a larger rear
    # tile (area 4) in one sightline. The rear tile is 4× bigger, so 1/4 of it is
    # covered (attenuated by exp(-0.5)) and 3/4 peeks around at full strength:
    #   I = I_front + (exp(-0.5)·1 + 3) = 1 + exp(-0.5) + 3 = 4 + exp(-0.5).
    # (origin/main computed `tmp_I += weights*unobscured_I` — intensity² — and only
    # evaluated obscuredFrac once; both fixed here.)
    pf(x, ΔA, τ) = BLR._RaytracePoint(1, 1, 1, 0.0, 0.0, 0.0, 0.0, x, 1.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, ΔA, τ, 1.0, false, false, Inf)
    pts = [pf(10.0, 1.0, 0.5), pf(0.0, 4.0, 0.0)]   # front (small), rear (large)
    res = BLR._rt_scan_bucket(pts, [1, 2], [1.0], 1.0, 1.0)   # outputΔA=1, τCutOff=1
    @test isapprox(res.I, 4.0 + exp(-0.5), rtol=1e-12)

    # ------------------------------------------------------------------
    # (7) zeroDiskObscuredClouds! / removeDiskObscuredClouds! with >2 submodels.
    # origin/main threw here (the `length(startInds) > 2` recursion sliced rings
    # without the camera/subModelStartInds and called zero...! positionally). The
    # logic classifies rings by type, so it now handles disk + N cloud groups.
    # ------------------------------------------------------------------
    obsc3() = disk(nr=6, nϕ=12) + clouds(30, 11) + clouds(30, 12, μ=750.0)
    z3 = BLR.zeroDiskObscuredClouds!(obsc3())
    @test !any(isnan, BLR.getVariable(z3, :I, flatten=true))
    r3 = BLR.removeDiskObscuredClouds!(obsc3())
    Ir3 = BLR.getVariable(r3, :I, flatten=true)
    @test !any(isnan, Ir3)
    @test length(Ir3) <= length(BLR.getVariable(obsc3(), :I, flatten=true))  # some clouds removed
end
