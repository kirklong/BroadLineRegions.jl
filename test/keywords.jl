using Test, Random, BroadLineRegions

uncheckedIntensity(; r, gain=1.0, kwargs...) = gain / r
registeredIntensity(; r, gain=1.0, kwargs...) = gain / r
BLR.acceptedKeywords(::typeof(registeredIntensity)) = (:gain,)
emptyIntensity(; r, kwargs...) = 1.0 / r
BLR.acceptedKeywords(::typeof(emptyIntensity)) = ()
registeredVelocity(; r, speed=1.0, kwargs...) = speed / sqrt(r)
BLR.acceptedKeywords(::typeof(registeredVelocity)) = (:speed,)
explodingIntensity(; kwargs...) = throw(DomainError(:original, "callback details"))
BLR.acceptedKeywords(::typeof(explodingIntensity)) = ()
contextIntensity(; θₒSystem, kwargs...) = θₒSystem
BLR.acceptedKeywords(::typeof(contextIntensity)) = ()

@testset "ASCII keyword normalization and validation" begin
    for (ascii, canonical) in BLR._ASCII_KEYWORDS
        @test keys(BLR._normalize_keywords((; ascii=>2.0))) == (canonical,)
        @test_throws ArgumentError BLR._normalize_keywords((; ascii=>2.0, canonical=>2.0))
    end
    common = (; F=0.3, i=0.4, rescale=0.8, parallel=false)
    ascii = (; mu=300.0, beta=0.7, theta_o=0.3, gamma=1.2, xi=0.6, r_s=2.0,
        eta_o=0.4, eta_1=0.7, alpha_rm=0.2, tau=0.3)
    unicode = (; μ=300.0, β=0.7, θₒ=0.3, γ=1.2, ξ=0.6, rₛ=2.0,
        ηₒ=0.4, η₁=0.7, αRM=0.2, τ=0.3)
    for rng in (MersenneTwister(13), :philox)
        seed = rng === :philox ? (; seed=77) : (;)
        a = BLR.cloudModel(16; rng=deepcopy(rng), seed..., common..., ascii...)
        u = BLR.cloudModel(16; rng=deepcopy(rng), seed..., common..., unicode...)
        @test isequal(a.params, u.params)
        for field in (:r, :v, :I, :η, :τ, :x, :y, :z)
            @test isequal(BLR.getVariable(a, field; flatten=true), BLR.getVariable(u, field; flatten=true))
        end
        @test isequal(BLR.getProfile(a, :line; bins=8).binSums, BLR.getProfile(u, :line; bins=8).binSums)
    end
    rng = MersenneTwister(19); untouched = copy(rng)
    for kw in ((; muu=1.), (; mu=1., μ=1.), (; θ=0.3), (; θ₀=0.3),
               (; r=1.), (; x=1.), (; constructor=:other), (; nClouds=99))
        @test_throws ArgumentError BLR.cloudModel(4; rng=rng, kw...)
        @test rand(copy(rng)) == rand(copy(untouched))
    end
    # Validation also precedes empty model paths and helper randomness.
    @test_throws ArgumentError BLR.cloudModel(0; typo=1.)
    @test_throws ArgumentError BLR.drawCloud(; rng=rng, typo=1.)
    @test rand(copy(rng)) == rand(copy(untouched))
    @test_throws ArgumentError BLR.vCircularCloud(; r=10., phi_0=0.4, i=0.3, rot=0.1, theta_o=0.2, typo=1)
    @test BLR.response(3.; eta_o=0.2, eta_1=0.7, alpha_rm=0.4) == BLR.response(3.; ηₒ=0.2, η₁=0.7, αRM=0.4)
    Random.seed!(56); a = BLR.DRW(; t=[0.,1.,2.], mu=0.7, tau=2., sigma=0.3)
    Random.seed!(56); u = BLR.DRW(; t=[0.,1.,2.], μ=0.7, τ=2., σ=0.3)
    @test a == u
    @test_logs (:warn, r"Unused.*κ.*σₜ") BLR.cloudModel(4; kappa=0.3, sigma_t=0.2)
    @test_logs BLR.cloudModel(4; rescale=0.3)
end

@testset "Optional callback declarations" begin
    args = (; rng=:philox, seed=14, parallel=false)
    @test BLR.acceptedKeywords(uncheckedIntensity) === nothing
    @test BLR.acceptedKeywords(emptyIntensity) === ()
    m = @test_logs (:warn, r"gain.*not checked.*acceptedKeywords") BLR.cloudModel(4; I=uncheckedIntensity, gain=3., args...)
    @test all(r.I == 3. / r.r for r in m.rings)
    @test m.params.gain == 3.
    @test_logs BLR.cloudModel(4; I=uncheckedIntensity, args...)
    m = @test_logs BLR.cloudModel(4; I=registeredIntensity, gain=4., args...)
    @test all(r.I == 4. / r.r for r in m.rings)
    @test_throws ArgumentError BLR.cloudModel(4; I=emptyIntensity, gain=4., args...)
    @test_throws ArgumentError BLR.cloudModel(4; I=registeredIntensity, gaain=4., args...)
    @test_logs BLR.cloudModel(4; I=registeredIntensity, v=registeredVelocity, gain=4., speed=2., args...)
    m = @test_logs (:warn, r"surprise.*not checked") BLR.cloudModel(4; I=uncheckedIntensity, v=registeredVelocity, surprise=1., args...)
    @test_logs (:warn, r"κ.*not checked") BLR.cloudModel(4; I=uncheckedIntensity, kappa=0.3, args...)
    err = try BLR.drawCloud(; I=explodingIntensity); catch e; e; end
    @test err isa DomainError
    @test occursin("callback details", sprint(showerror, err))
    @test BLR.drawCloud(; I=contextIntensity, theta_o_system=0.35).I == 0.35
    # Direct ring callbacks use the same rules.
    @test_throws ArgumentError BLR.ring(; r=10., i=0.3, v=1., I=emptyIntensity, phi=0.4, delta_a=1., gain=2.)
end

@testset "Forwarding and reconstruction" begin
    c = BLR.drawCloud(; rng=MersenneTwister(1), r_s=4., phi_0=0.6, theta_o=0.2, rot=0.4, i=0.5)
    @test c.v == BLR.vCircularCloud(; r=c.r, rₛ=4., ϕ₀=0.6, θₒ=0.2, rot=0.4, i=0.5, reflect=c.reflect)
    opts = (; nr=5, nϕ=12, f1=1., f2=0.3, f3=0.2, f4=0.1)
    a = BLR.DiskWindModel(50., 2., 1., 0.4; rot=0.3, opts...)
    lo, hi = BLR.get_rMinMaxDiskWind(50., 2., 1.)
    b = BLR.DiskWindModel(lo, hi, 0.4; rot=0.3, alpha=1., opts...)
    for field in (:r, :v, :I, :x, :y, :z)
        @test isequal(BLR.getVariable(a, field; flatten=true), BLR.getVariable(b, field; flatten=true))
    end
    for ring in a.rings
        for j in eachindex(ring.r)
            if isfinite(ring.r[j])
                xyz = BLR.rotate3D_scalar(ring.r[j], ring.ϕ₀[j], 0.4, 0.3, 0.0)
                @test [ring.x[j],ring.y[j],ring.z[j]] ≈ collect(xyz)
            end
        end
    end
    @test BLR.rebuild(a; r_bar=60.).params.r̄ == 60.
    @test BLR.rebuild(a+b; alpha=0.8).params.left.α == 0.8
    cm = BLR.CompositeModel(a; line="A", lineCenter=1.)
    BLR.addLine!(cm; line="B", lineCenter=2., r_bar=60.)
    @test cm["B"].params.r̄ == 60.
    for (fun, extra) in ((BLR.getSpectrum, (; minX=0., maxX=4.)),
                         ((x; kw...)->BLR.getProfile(x,:ratio; lines=("A","B"), kw...), (;)))
        @test_throws ArgumentError fun(cm; typo=1., extra...)
    end
    @test_throws ArgumentError BLR.raytrace!(a; tau_cutoff=0.4, typo=1.)
    @test_throws ArgumentError BLR.rebuild(a; alpha=0.8, α=0.8)
    kw = (; i=0.4, f1=1., f2=0.3, f3=0.2, f4=0.1, alpha=1., phi_min=-1., phi_max=1., rMin=10., rMax=20.)
    @test BLR.IϕDiskWindMask(; r=5., phi=0., kw...) == 0.
    @test BLR.IϕDiskWindMask(; r=[5.,15.,25.], phi=[0.,0.,0.], kw...) == [0., BLR.DiskWindIntensity(; r=15., ϕ=0., i=0.4, f1=1., f2=0.3, f3=0.2, f4=0.1, α=1., rMin=10., rMax=20.), 0.]
    @test_throws ArgumentError BLR.phase(a; U=Float64[], V=Float64[], PA=0., BLRAng=1., typo=1.)
    @test_throws ArgumentError BLR.secondMoment(a; U=Float64[], V=Float64[], PA=0., BLRAng=1., typo=1.)
end

@testset "Resident keyword boundaries" begin
    a = BLR.residentCloudModel(16, 12; mu=200., beta=0.7, theta_o=0.3, r_s=2., eta_o=0.4, tau=0.2)
    b = BLR.residentCloudModel(16, 12; μ=200., β=0.7, θₒ=0.3, rₛ=2., ηₒ=0.4, τ=0.2)
    @test a.ma.r == b.ma.r
    @test a.ma.v == b.ma.v
    @test_throws ArgumentError BLR.residentCloudModel(0,12; typo=1.)
    @test_throws ArgumentError BLR.residentDiskWindModel(50.,2.,1.,0.4; typo=1.)
    @test_throws ArgumentError BLR.gpuCloudModel(0,12; typo=1.)
    @test_throws ArgumentError BLR.gpuDiskWindModel(50.,2.,1.,0.4; typo=1.)
    @test_throws ArgumentError BLR.raytrace!(a; tau_cutoff=1., typo=1.)
    @test_throws ArgumentError BLR.phase(a; U=Float64[],V=Float64[],PA=0.,BLRAng=1.,typo=1.)
end

@testset "Recorded parameters, forwarding skips, and unchanged inputs" begin
    m = BLR.cloudModel(8; rng=:philox, seed=4, I=registeredIntensity, gain=2., parallel=false)
    @test_logs BLR.rebuild(m; gain=3.)
    @test BLR.rebuild(m; gain=3.).params.gain == 3.
    @test_throws ArgumentError BLR.rebuild(merge(m.params, (; typo=1.)))
    @test_throws ArgumentError BLR.rebuild(merge(m.params, (; I=emptyIntensity)))
    u = @test_logs (:warn,) BLR.cloudModel(8; rng=:philox, seed=4, I=uncheckedIntensity, gain=2., parallel=false)
    @test_logs (:warn,) BLR.rebuild(u; gain=3.)
    # A newly registered callback makes old unchecked extras strict.
    @test_throws ArgumentError BLR.rebuild(merge(u.params, (; I=emptyIntensity)))
    saved = deepcopy(m)
    @test_throws ArgumentError BLR.raytrace!(m; misspelled=1.)
    @test isequal(m.params, saved.params)
    @test m.camera.raytraced == saved.camera.raytraced
    @test all(isequal(m.rings[j].I, saved.rings[j].I) for j in eachindex(m.rings))
    cm = BLR.CompositeModel(m; line="A", lineCenter=1.)
    @test_throws ErrorException BLR.addLine!(cm; line="B", lineCenter=2., typo=1.)
    @test cm.lines == ["A"]
    d = BLR.DiskWindModel(300.,900.,0.4; nr=4,n_phi=8,I=BLR.IsotropicIntensity, tau=0.4)
    @test_throws ArgumentError BLR.model(300.,900.,0.4,4,8,BLR.IsotropicIntensity,BLR.vCircularDisk,:linear; r=1.)
    @test_throws ArgumentError BLR.DiskWindModel(300.,900.,0.4; r=1.)
    @test_throws ArgumentError BLR.DiskWindModel(300.,900.,0.4; theta_o=0.4)
    traced = BLR.raytrace!(d + m; tau_cutoff=0.8)
    rebuilt = BLR.rebuild(traced; gain=3.)
    @test rebuilt.params.parent.right.gain == 3.
    rcm = BLR.resident(BLR.CompositeModel(d; line="A", lineCenter=1.))
    BLR.addLine!(cm, d; line="B", lineCenter=2.)
    rcm2 = BLR.resident(cm)
    @test_throws ArgumentError BLR.getSpectrum(rcm; minX=0.,maxX=2.,typo=1.)
    @test_throws ArgumentError BLR.getProfile(rcm2,:ratio; lines=("A","B"),typo=1.)
    @test_throws ArgumentError BLR.getProfile(d,:line; U=Float64[])
    @test_throws ArgumentError BLR.binModel(8; m=d,yVariable=:I,typo=1.)
    @test_throws ArgumentError BLR.binnedVariance([0.,1.],[0.,1.],[1.,1.]; typo=1.)
    @test_throws ArgumentError BLR.binWeightedMean(d,[0.,1.],[1.,1.],[1.,1.],8,nothing; typo=1.)
end

@testset "Custom adapters delegate shared context without repeated warnings" begin
    diskadapter(; r, kwargs...) = BLR.DiskWindIntensity(; r=r, kwargs...) .* (300. ./ r)
    opts = (; f1=1., f2=0., f3=0., f4=0., alpha=1., nr=4, n_phi=8)
    @test_logs (:warn, r"not checked") BLR.DiskWindModel(300.,900.,0.4; I=diskadapter, opts...)
    # The callback context must be restored, including after callback exceptions.
    @test_throws ArgumentError BLR.DiskWindIntensity(; typo=1.)
    @test_throws DomainError BLR.drawCloud(; I=explodingIntensity)
    @test_throws ArgumentError BLR.IsotropicIntensity(; r=1., phi=0., typo=1.)
    @test_logs BLR.cloudModel([0.2], [0.4], [0.3], [0.1], 0.2, 1.; rng=MersenneTwister(1))
end

@testset "Turbulent velocity aliases and explicit inactive resident options" begin
    a = (; sigma_rho_r=0.1, sigma_rho_c=0.2, sigma_theta_r=0.3, sigma_theta_c=0.4,
        theta_e=0.5, sigma_t=0.05, kappa=0.2)
    u = (; σρᵣ=0.1, σρc=0.2, σΘᵣ=0.3, σΘc=0.4, θₑ=0.5, σₜ=0.05, κ=0.2)
    opts = (; I=BLR.cloudIntensity, v=BLR.vCloudTurbulentEllipticalFlow,
        fEllipse=0.7, fFlow=0.6, i=0.4, rng=:philox, seed=90, parallel=false)
    ma = BLR.cloudModel(8; a..., opts...)
    mu = BLR.cloudModel(8; u..., opts...)
    @test ma.params == mu.params
    @test [r.v for r in ma.rings] == [r.v for r in mu.rings]
    @test [r.I for r in ma.rings] == [r.I for r in mu.rings]
    @test_logs (:warn, r"Unused.*κ.*σₜ") BLR.residentCloudModel(4, 7; kappa=0.3, sigma_t=0.2)
    @test_logs BLR.residentCloudModel(4,7)
    @test BLR.vCircularRadialDisk(;r=10.,i=0.4,phi=0.3,v_r_frac=0.2,r_s=2.) ==
        BLR.vCircularRadialDisk(;r=10.,i=0.4,ϕ=0.3,vᵣFrac=0.2,rₛ=2.)
end

@testset "Empty composite forwarding validates before skipping work" begin
    cm = BLR.CompositeModel(BLR.cloudModel(2; rng=:philox,seed=1); line="A",lineCenter=1.)
    empty!(cm.lines)
    @test_throws ArgumentError BLR.raytrace!(cm; typo=1.)
    @test_throws ArgumentError BLR.reset!(cm; typo=1.)
    @test_throws ArgumentError BLR.resident(cm; typo=1.)
    @test_throws ArgumentError BLR.gpu(cm; typo=1.)
end

scaleContextVelocity(; rₛ, kwargs...) = rₛ
BLR.acceptedKeywords(::typeof(scaleContextVelocity)) = ()
lateRegisteredIntensity(; r, gain=1., kwargs...) = gain/r
@testset "Registration transitions and grouped reconstruction warnings" begin
    @test BLR.drawCloud(; r_s=3., v=scaleContextVelocity).v == 3.
    opts = (; I=uncheckedIntensity,gain=2.,rng=:philox,seed=4,parallel=false)
    m = @test_logs (:warn,) BLR.cloudModel(2;opts...)
    @test_logs (:warn,) BLR.rebuild(m+m;gain=3.)
    @test_throws ArgumentError BLR.rebuild(merge((m+m).params, (; typo=1.)))
    @test_logs (:warn,) BLR.cloudModel(2;I=lateRegisteredIntensity,gain=2.,rng=:philox,seed=4)
    @eval BLR.acceptedKeywords(::typeof(lateRegisteredIntensity)) = ()
    @test_throws ArgumentError Base.invokelatest(BLR.cloudModel,2;I=lateRegisteredIntensity,gain=2.,rng=:philox,seed=4)
    @test_logs (:warn, r"Unused.*κ") BLR.IsotropicIntensity(;r=1.,phi=0.,kappa=0.3)
end
