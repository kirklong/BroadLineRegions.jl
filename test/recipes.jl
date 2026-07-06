using RecipesBase
using KernelAbstractions
using Plots # the W4-T7 CompositeModel testset checks recipe entry points that return Plots.Plot objects

# Apply a @userplot recipe without a plotting backend and gather, per attribute, the concatenated &
# sorted values across all generated series. A scatter recipe is order-independent (it's a point
# cloud), so comparing the sorted multiset of each coordinate is a robust host-vs-resident check.
function _recipe_series(recipeObj)
    return RecipesBase.apply_recipe(Dict{Symbol,Any}(), recipeObj)
end
function _gather(rds, key)
    vals = Float64[]
    for rd in rds
        pa = rd.plotattributes
        haskey(pa, key) && append!(vals, Float64.(pa[key]))
    end
    return sort(vals)
end
# NaN-aware ≈ for the sorted value multisets (image plots unmasked points, so marker_z may carry NaN).
_nanapprox(a, b; rtol=1e-10) = length(a) == length(b) &&
    all(isnan(b[i]) ? isnan(a[i]) : isapprox(a[i], b[i]; rtol=rtol) for i in eachindex(a, b))

@testset "Plots recipes: ResidentModel parity with host model" begin
    disk = BLR.DiskWindModel(300.0, 900.0, 0.4, nr=12, nϕ=24, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.7, f3=0.2, f4=0.9,
        α=1.2, ηₒ=0.4, η₁=0.6, αRM=0.1, rNorm=700.0, τ=0.4, reflect=false)
    clouds = BLR.cloudModel(150, μ=600.0, β=1.0, F=0.5, θₒ=0.4, i=0.4, γ=1.0, ξ=0.8,
        I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=0.1, rng=MersenneTwister(41))

    @testset "image" begin
        for m in (disk, clouds, BLR.raytrace!(disk + clouds))
            rm = BLR.resident(m; backend=KernelAbstractions.CPU())
            for variable in (:I, :v, :r)
                host = _recipe_series(BLR.Image((m, variable)))
                res  = _recipe_series(BLR.Image((rm, variable)))
                @test _nanapprox(_gather(host, :x), _gather(res, :x))
                @test _nanapprox(_gather(host, :y), _gather(res, :y))
                @test _nanapprox(_gather(host, :marker_z), _gather(res, :marker_z))
            end
        end
        # default variable is :I when only the model is passed (1-arg path)
        rmDisk = BLR.resident(disk; backend=KernelAbstractions.CPU())
        @test _nanapprox(_gather(_recipe_series(BLR.Image((disk,))), :marker_z),
                         _gather(_recipe_series(BLR.Image((rmDisk,))), :marker_z))
    end

    @testset "plot3d" begin
        # single-submodel handles need no raytrace metadata; combined ones carry per-point submodel
        # info via rt (raytrace=true) so the recipe can split series by submodel.
        cases = ((disk, false), (clouds, false), (BLR.raytrace!(disk + clouds), true))
        for (m, rt) in cases
            rm = BLR.resident(m; backend=KernelAbstractions.CPU(), raytrace=rt)
            for variable in (:I, :r)
                host = _recipe_series(BLR.Plot3d((m, variable, false)))   # annotate=false: scatter only
                res  = _recipe_series(BLR.Plot3d((rm, variable, false)))
                @test _nanapprox(_gather(host, :x), _gather(res, :x))
                @test _nanapprox(_gather(host, :y), _gather(res, :y))
                @test _nanapprox(_gather(host, :z), _gather(res, :z))
                @test _nanapprox(_gather(host, :marker_z), _gather(res, :marker_z))
            end
        end
        # a combined ResidentModel without raytrace metadata cannot recover submodel boundaries
        rmNoMeta = BLR.resident(disk + clouds; backend=KernelAbstractions.CPU())
        @test_throws ArgumentError _recipe_series(BLR.Plot3d((rmNoMeta, :I, false)))
    end

    @testset "profile" begin
        rm = BLR.resident(disk; backend=KernelAbstractions.CPU())
        # recipe series should reproduce getProfile on the resident handle
        ref = BLR.getProfile(rm, :line)
        rds = _recipe_series(BLR.Profile((rm, :line)))
        @test length(rds) == 1
        @test rds[1].plotattributes[:x] ≈ ref.binCenters
        @test rds[1].plotattributes[:y] ≈ ref.binSums           # single profile is unnormalized
        # a ResidentModel stores no preset profiles, so a bare profile(rm) must error
        @test_throws ErrorException _recipe_series(BLR.Profile((rm,)))
    end
end

@testset "Plots recipes: CompositeModel (W4-T7)" begin
    # Tiny two-line composite sharing one geometry (mirrors the style of the W4-T4/T5/T6 testset in
    # runtests.jl) -- deliberately small, this testset should add seconds, not minutes.
    mDisk = BLR.DiskWindModel(300.0, 900.0, 0.4; nr=6, nϕ=12, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.0, f3=0.0, f4=0.0, α=1.0,
        τ=0.4, reflect=false)
    cm = BLR.CompositeModel(mDisk; line="Ha", lineCenter=6562.8)
    BLR.addLine!(cm, mDisk; line="Hb", lineCenter=4861.3, fluxRatio=0.35)

    @testset "profile: bare `profile` struct" begin
        p = BLR.getProfile(mDisk, :line, bins=21)
        rds = _recipe_series(BLR.Profile((p,)))
        @test length(rds) == 1
        @test rds[1].plotattributes[:x] == p.binCenters
        @test rds[1].plotattributes[:y] == p.binSums   # bare struct is plotted unnormalized
    end

    @testset "profile: CompositeModel" begin
        # default variable (:line): one series per line, labeled by line name, each normalized by its
        # own max |binSums| (so shapes are comparable regardless of fluxRatio).
        rds = _recipe_series(BLR.Profile((cm,)))
        @test length(rds) == length(cm.lines)
        for (i, line) in enumerate(cm.lines)
            pLine = BLR.getProfile(cm, :line; line=line)
            norm = maximum(abs(b) for b in pLine.binSums if !isnan(b))
            @test rds[i].plotattributes[:label] == line
            @test rds[i].plotattributes[:x] == pLine.binCenters
            @test rds[i].plotattributes[:y] ≈ pLine.binSums ./ norm
        end
        # an explicit profile name is honored too
        rdsR = _recipe_series(BLR.Profile((cm, :line)))
        @test length(rdsR) == length(cm.lines)
        # :ratio (W5) plots through the composite recipe with the lines pair passed POSITIONALLY
        # (`lines` is a Plots magic-attribute alias, so a keyword would never reach the recipe),
        # rendering the single unnormalized ratio series -- identical to plotting the bare struct
        pR = BLR.getProfile(cm, :ratio; lines=("Ha", "Hb"))
        rdsRatio = _recipe_series(BLR.Profile((cm, :ratio, ("Ha", "Hb"))))
        @test length(rdsRatio) == 1
        @test rdsRatio[1].plotattributes[:x] == pR.binCenters
        @test isequal(rdsRatio[1].plotattributes[:y], pR.binSums) #isequal: empty-denominator bins are NaN
        # the full Plots pipeline call works too (this is the user-facing entry point)
        pltRatio = BLR.profile(cm, :ratio, ("Ha", "Hb"))
        @test pltRatio isa Plots.Plot
        @test isequal(collect(pltRatio.series_list[1][:y]), pR.binSums)
        # without the lines pair the getProfile error still surfaces, naming the kwarg
        errRatio = try; _recipe_series(BLR.Profile((cm, :ratio))); nothing; catch e; e; end
        @test errRatio isa ErrorException && occursin("lines", errRatio.msg)
        # the pair is rejected for any other profile name (both arms of the guard)
        @test_throws ErrorException _recipe_series(BLR.Profile((cm, :line, ("Ha", "Hb"))))
        @test_throws ErrorException _recipe_series(BLR.Profile((mDisk, :line, ("Ha", "Hb"))))
    end

    @testset "profile: degenerate (zero-emission) line does not throw" begin
        # adversarial review 2026-07-03: a zero-emission line has all-zero :line bins (norm would be
        # 0 -> divide by zero) and all-NaN intensity-weighted :r bins (norm's NaN-filtered iterator
        # would be empty -> maximum throws). Both fall back to norm = 1 (plotted unnormalized) so the
        # other lines still render.
        mZero = BLR.DiskWindModel(300.0, 900.0, 0.4; nr=6, nϕ=12, scale=:linear,
            I=BLR.IsotropicIntensity, rescale=0.0, v=BLR.vCircularDisk, τ=0.4, reflect=false)
        cmZ = BLR.CompositeModel(mDisk; line="Ha", lineCenter=6562.8)
        BLR.addLine!(cmZ, mZero; line="dark", lineCenter=4861.3)
        rdsLine = _recipe_series(BLR.Profile((cmZ,))) #:line -- "dark" bins are all zero
        @test length(rdsLine) == 2
        @test all(iszero, rdsLine[2].plotattributes[:y]) #unnormalized zeros, not 0/0 = NaN
        rdsR = _recipe_series(BLR.Profile((cmZ, :r))) #:r -- "dark" bins are all NaN (0/0 weighted mean)
        @test length(rdsR) == 2
        @test all(isnan, rdsR[2].plotattributes[:y]) #plotted as-is (Plots skips NaN points)
    end

    @testset "spectrum" begin
        # no overlap: per-line series + one black "total" series
        sp = BLR.spectrum(cm; bins=40)
        @test sp isa Plots.Plot
        @test length(sp.series_list) == length(cm.lines) + 1
        labels = [s[:label] for s in sp.series_list]
        @test all(line -> line in labels, cm.lines)
        @test "total" in labels
        totalSeries = only(filter(s -> s[:label] == "total", sp.series_list))
        @test totalSeries[:linewidth] == 2

        # overlapping lines add one shaded band series per overlapping pair
        cA, cB = 6563.0, 6600.0
        cmO = BLR.CompositeModel(mDisk; line="A", lineCenter=cA)
        BLR.addLine!(cmO, mDisk; line="B", lineCenter=cB)
        ov = BLR.lineOverlap(cmO)
        @test length(ov) == 1
        spO = BLR.spectrum(cmO; bins=40)
        @test length(spO.series_list) == length(cmO.lines) + 1 + length(ov)
    end

    @testset "image forwarding" begin
        im = BLR.image(cm, :I; line="Ha")
        @test im isa Plots.Plot
        @test length(im.series_list) == 1
        # `line` has no default -- an intensity image mixing lines' arbitrary units is not meaningful
        @test_throws UndefKeywordError BLR.image(cm, :I)
    end

    @testset "plot3d forwarding" begin
        # line=<name> is identical to plotting that line's model directly
        pLine = BLR.plot3d(cm; line="Ha")
        pDirect = BLR.plot3d(mDisk)
        @test pLine isa Plots.Plot
        @test length(pLine.series_list) == length(pDirect.series_list)

        # line=nothing (default) overlays every line: one scatter3d series per line (single-submodel
        # models here) plus one camera-annotation series (annotate=true default)
        pAll = BLR.plot3d(cm)
        @test pAll isa Plots.Plot
        @test length(pAll.series_list) == length(cm.lines) + 1

        # annotate=false drops the camera series
        pNoAnnotate = BLR.plot3d(cm, false)
        @test length(pNoAnnotate.series_list) == length(cm.lines)
    end
end

@testset "Plots recipes: ResidentCompositeModel (W4-G3)" begin
    # Same tiny two-line fixture as the CompositeModel testset above; the resident composite (Float64
    # on the CPU backend) must drive the same recipes with matching series.
    mDisk = BLR.DiskWindModel(300.0, 900.0, 0.4; nr=6, nϕ=12, scale=:linear,
        I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.0, f3=0.0, f4=0.0, α=1.0,
        τ=0.4, reflect=false)
    cm = BLR.CompositeModel(mDisk; line="Ha", lineCenter=6562.8)
    BLR.addLine!(cm, mDisk; line="Hb", lineCenter=4861.3, fluxRatio=0.35)
    rcm = BLR.resident(cm; backend=KernelAbstractions.CPU())

    @testset "profile: series match the CompositeModel recipe" begin
        rdsC = _recipe_series(BLR.Profile((cm,)))
        rdsR = _recipe_series(BLR.Profile((rcm,)))
        @test length(rdsR) == length(cm.lines)
        for i in eachindex(rdsC)
            @test rdsR[i].plotattributes[:label] == rdsC[i].plotattributes[:label]
            # identical default edges (deterministic constructBinEdges from the same Float64 range)
            @test rdsR[i].plotattributes[:x] == rdsC[i].plotattributes[:x]
            @test isapprox(rdsR[i].plotattributes[:y], rdsC[i].plotattributes[:y], rtol=1e-10)
        end
        # an explicit profile name is honored too
        rdsR2 = _recipe_series(BLR.Profile((rcm, :line)))
        @test length(rdsR2) == length(cm.lines)
        # :ratio (W5) drives the same recipe branch as the host composite, via the resident getProfile
        rdsRatioC = _recipe_series(BLR.Profile((cm, :ratio, ("Ha", "Hb"))))
        rdsRatioR = _recipe_series(BLR.Profile((rcm, :ratio, ("Ha", "Hb"))))
        @test length(rdsRatioR) == 1
        @test rdsRatioR[1].plotattributes[:x] == rdsRatioC[1].plotattributes[:x]
        @test _nanapprox(collect(rdsRatioR[1].plotattributes[:y]), collect(rdsRatioC[1].plotattributes[:y]))
        # without the lines pair the resident getProfile error still surfaces
        @test_throws ErrorException _recipe_series(BLR.Profile((rcm, :ratio)))
    end

    @testset "spectrum" begin
        # same recipe body as the host composite: per-line series + black "total"
        sp = BLR.spectrum(rcm; bins=40)
        @test sp isa Plots.Plot
        @test length(sp.series_list) == length(cm.lines) + 1
        labels = [s[:label] for s in sp.series_list]
        @test all(line -> line in labels, cm.lines)
        @test "total" in labels
        # the total series matches the host composite spectrum (Float64 CPU backend)
        spC = BLR.spectrum(cm; bins=40)
        totR = only(filter(s -> s[:label] == "total", sp.series_list))
        totC = only(filter(s -> s[:label] == "total", spC.series_list))
        @test isapprox(collect(totR[:y]), collect(totC[:y]), rtol=1e-10)

        # overlapping lines shade a band, driven by lineOverlap(::ResidentCompositeModel)
        cmO = BLR.CompositeModel(mDisk; line="A", lineCenter=6563.0)
        BLR.addLine!(cmO, mDisk; line="B", lineCenter=6600.0)
        rcmO = BLR.resident(cmO; backend=KernelAbstractions.CPU())
        ov = BLR.lineOverlap(rcmO)
        @test ov == BLR.lineOverlap(cmO) && length(ov) == 1
        spO = BLR.spectrum(rcmO; bins=40)
        @test length(spO.series_list) == length(cmO.lines) + 1 + length(ov)
    end

    @testset "image forwarding" begin
        im = BLR.image(rcm, :I; line="Ha")
        @test im isa Plots.Plot
        @test length(im.series_list) == 1
        # `line` has no default, matching the host composite forwarding
        @test_throws UndefKeywordError BLR.image(rcm, :I)
    end

    @testset "plot3d forwarding" begin
        # line=<name> is identical to plotting that line's resident handle directly
        pLine = BLR.plot3d(rcm; line="Ha")
        pDirect = BLR.plot3d(rcm["Ha"])
        @test pLine isa Plots.Plot
        @test length(pLine.series_list) == length(pDirect.series_list)

        # line=nothing (default) overlays every line (+ camera annotation), like the host composite
        pAll = BLR.plot3d(rcm)
        @test pAll isa Plots.Plot
        @test length(pAll.series_list) == length(cm.lines) + 1
        pNoAnnotate = BLR.plot3d(rcm, false)
        @test length(pNoAnnotate.series_list) == length(cm.lines)
    end
end
