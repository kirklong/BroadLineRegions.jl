using RecipesBase
using KernelAbstractions

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
