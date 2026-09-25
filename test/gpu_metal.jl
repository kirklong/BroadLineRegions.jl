# Apple-silicon (Metal) GPU correctness tests. These only run when BLR_TEST_METAL=1 AND a functional
# Metal device is present (see the include guard in runtests.jl). Mirrors test/gpu_cuda.jl, but Metal
# has no Float64, so every device run is Float32 and is validated against the CPU backend run on the
# SAME Float32 inputs (tight: identical binning, only accumulation order / device-libm ULPs differ) and,
# where meaningful, against the Float64 host reference at Float32 tolerances.
using Metal
using KernelAbstractions

# Kernels for the Philox known-answer test (defined at top level: @kernel cannot live in a testset).
@kernel function philox_kat_kernel!(out, c, k)
    i = @index(Global)
    o = BLR._philox4x32_bijection(c[4i-3], c[4i-2], c[4i-1], c[4i], k[2i-1], k[2i])
    out[4i-3] = o[1]; out[4i-2] = o[2]; out[4i-1] = o[3]; out[4i] = o[4]
end
@kernel function uniform_kernel!(u, key0, key1)
    i = @index(Global)
    u[i] = BLR._cloud_uniform(Float32, key0, key1, UInt32((i - 1) ÷ 16 + 1), (i - 1) % 16)
end


@testset "Metal GPU raytrace + kernels" begin
    if !Metal.functional()
        @warn "BLR_TEST_METAL=1 but Metal.functional() is false -- skipping Metal GPU tests"
    else
        @info "running Metal GPU tests on $(Metal.device())"
        backend = Metal.MetalBackend()
        cpuB = KernelAbstractions.CPU()

        function approx_eq(a, b; rtol=1e-4, atol=1e-6)
            length(a) == length(b) || return false
            for i in eachindex(a, b)
                if isnan(b[i])
                    isnan(a[i]) || return false
                elseif !isapprox(a[i], b[i]; rtol=rtol, atol=atol)
                    return false
                end
            end
            return true
        end
        # scale-aware Float32 comparison: atol relative to the largest finite reference value
        relclose(a, b; rtol=1e-4) = approx_eq(a, b; rtol=rtol,
            atol=rtol * maximum(abs, filter(isfinite, b); init=0.0) + 1e-30)

        disk(; r1=300.0, r2=900.0, nr=8, nϕ=16, inc=0.4, τ=0.4) =
            BLR.DiskWindModel(r1, r2, inc, nr=nr, nϕ=nϕ, scale=:linear,
                I=BLR.IsotropicIntensity, v=BLR.vCircularDisk, τ=τ, reflect=false)
        clouds(n, seed; μ=600.0, inc=0.4, τ=0.1) =
            BLR.cloudModel(n, μ=μ, β=1.0, F=0.5, θₒ=0.4, i=inc, γ=1.0, ξ=1.0,
                I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=τ, rng=MersenneTwister(seed))

        builders = (
            () -> disk() + clouds(80, 101),
            () -> clouds(80, 101) + disk(),
            () -> disk(r1=250.0, r2=700.0, nr=6, nϕ=12) + disk(r1=500.0, r2=1000.0, nr=6, nϕ=12),
            () -> disk() + clouds(50, 102) + clouds(50, 103, μ=850.0),
        )

        @testset "backend selection + Float64 guard" begin
            if Base.get_extension(BLR, :BroadLineRegionsCUDAExt) === nothing
                @test BLR.defaultGPUBackend() isa Metal.MetalBackend
            end
            m = builders[1]()
            gm = BLR.gpu(m; backend=backend)
            @test gm.backend isa Metal.MetalBackend
            @test gm.ma.I isa Metal.MtlVector{Float32}
            @test gm.ma.reflect isa Metal.MtlVector{Bool}
            @test gm.rt isa BLR.RaytraceMeta
            # no FP64 on Apple GPUs: every entry point rejects T=Float64 with an actionable error
            @test_throws ArgumentError BLR.gpu(m; backend=backend, T=Float64)
            @test_throws ArgumentError BLR.gpu(BLR.flatten(m); backend=backend)
            @test_throws ArgumentError BLR.gpuCloudModel(100, 1; backend=backend, T=Float64)
            @test_throws ArgumentError BLR.gpuDiskWindModel(300.0, 900.0, 0.4; backend=backend, T=Float64,
                f1=1.0, f2=1.0, f3=1.0, f4=1.0, α=1.0)
            @test_throws ArgumentError BLR.raytrace!(builders[1](); backend=backend)   # T defaults to Float64
        end

        @testset "flatten/gpu round-trip" begin
            m = builders[1]()
            ma32 = BLR.flatten(m; T=Float32)
            gm = BLR.gpu(m; backend=backend)
            back = BLR.cpu(gm).ma
            @test isequal(back.I, ma32.I)
            @test isequal(back.x, ma32.x)
            @test back.α == ma32.α
            @test back.reflect == ma32.reflect
        end

        @testset "Philox4x32-10 known-answer test (device == CPU == Random123 KAT)" begin
            # Published Random123 kat_vectors for philox4x32_10: (counter, key) -> output.
            kat = [
                ((0x00000000, 0x00000000, 0x00000000, 0x00000000), (0x00000000, 0x00000000),
                 (0x6627e8d5, 0xe169c58d, 0xbc57ac4c, 0x9b00dbd8)),
                ((0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff), (0xffffffff, 0xffffffff),
                 (0x408f276d, 0x41c83b0e, 0xa20bc7c6, 0x6d5451fd)),
                ((0x243f6a88, 0x85a308d3, 0x13198a2e, 0x03707344), (0xa4093822, 0x299f31d0),
                 (0xd16cfe09, 0x94fdcceb, 0x5001e420, 0x24126ea1)),
            ]
            # KAT vectors + a batch of pseudo-random (counter, key) pairs
            rng = MersenneTwister(77)
            cs = vcat([collect(UInt32, t[1]) for t in kat]..., rand(rng, UInt32, 4 * 256))
            ks = vcat([collect(UInt32, t[2]) for t in kat]..., rand(rng, UInt32, 2 * 256))
            n = length(cs) ÷ 4
            run(b, C, K) = (o = KernelAbstractions.zeros(b, UInt32, 4n);
                philox_kat_kernel!(b)(o, C, K; ndrange=n); KernelAbstractions.synchronize(b); Array(o))
            outCPU = run(cpuB, cs, ks)
            outGPU = run(backend, MtlArray(cs), MtlArray(ks))
            @test outGPU == outCPU                                   # bit-identical, all 259 blocks
            for (j, t) in enumerate(kat)
                @test Tuple(outGPU[4j-3:4j]) == t[3]                  # matches the published KAT
            end

            # The Float32 uniforms the cloud kernel consumes are bit-identical too.
            uC = zeros(Float32, 16 * 512); uG = Metal.zeros(Float32, 16 * 512)
            uniform_kernel!(cpuB)(uC, 0x0000abcd, 0x00001234; ndrange=length(uC))
            uniform_kernel!(backend)(uG, 0x0000abcd, 0x00001234; ndrange=length(uG))
            KernelAbstractions.synchronize(backend)
            @test Array(uG) == uC
            @test all(0 .< uC .< 1)
        end

        @testset "extension sort matches CPU reference" begin
            keys = [2, 1, 2, 1, 0, 2, 0, 2, 1]
            x = Float32[10.0, 4.0, 11.0, 6.0, 3.0, NaN, 9.0, 11.0, 4.0]  # depth ties (two 11.0, two 4.0)
            permCPU = BLR._rt_sortperm_by_key_depth(keys, x)
            permGPU = BLR._rt_sortperm_by_key_depth(MtlArray(keys), MtlArray(x))
            @test permGPU isa Metal.MtlVector{Int}
            permGPU = Array(permGPU)
            @test BLR._rt_sorted_key_depth_pairs(keys, x, permGPU) ==
                  BLR._rt_sorted_key_depth_pairs(keys, x, permCPU)
            @test permGPU == permCPU
            # larger random case with many ties
            rng = MersenneTwister(5)
            k2 = rand(rng, 0:40, 5000); x2 = Float32.(rand(rng, 1:50, 5000)); x2[1:37:end] .= NaN32
            @test Array(BLR._rt_sortperm_by_key_depth(MtlArray(k2), MtlArray(x2))) ==
                  BLR._rt_sortperm_by_key_depth(k2, x2)
        end

        @testset "host raytrace!(backend=Metal, T=Float32) ≈ CPU" begin
            for build in builders
                ref = BLR.raytrace!(build(); τCutOff=1.0)                                       # host Float64
                ref32 = BLR.raytrace!(build(); τCutOff=1.0, backend=cpuB, T=Float32)             # CPU backend, same Float32 inputs
                gpu = BLR.raytrace!(build(); τCutOff=1.0, backend=backend, T=Float32)
                @test gpu.subModelStartInds == ref32.subModelStartInds
                for sym in (:I, :v, :r, :ϕ, :ϕ₀, :τ, :η)
                    @test approx_eq(BLR.getVariable(gpu, sym, flatten=true),
                        BLR.getVariable(ref32, sym, flatten=true); rtol=1e-5, atol=1e-6)
                end
                @test approx_eq(gpu.camera.α, ref32.camera.α; rtol=1e-6, atol=1e-6)
                # and the Float32 device result tracks the Float64 host reference at Float32 precision
                for sym in (:I, :v, :r)
                    @test approx_eq(BLR.getVariable(gpu, sym, flatten=true),
                        BLR.getVariable(ref, sym, flatten=true); rtol=1e-4, atol=1e-6)
                end
                refI = BLR.getVariable(ref, :I, flatten=true); refA = BLR.getVariable(ref, :ΔA, flatten=true)
                gpuI = BLR.getVariable(gpu, :I, flatten=true); gpuA = BLR.getVariable(gpu, :ΔA, flatten=true)
                @test isapprox(sum(gpuI .* gpuA), sum(refI .* refA); rtol=1e-5)
            end
            ref = BLR.raytrace!(builders[1](); IRatios=[1.0, 0.25], τCutOff=1.0)
            gpu = BLR.raytrace!(builders[1](); IRatios=[1.0, 0.25], τCutOff=1.0, backend=backend, T=Float32)
            for sym in (:I, :v, :r)
                @test approx_eq(BLR.getVariable(gpu, sym, flatten=true),
                    BLR.getVariable(ref, sym, flatten=true); rtol=1e-4, atol=1e-6)
            end
        end

        @testset "velocity-dependent τ errors on Metal path" begin
            d = disk()
            for r in d.rings
                r.τ = fill(0.4, length(r.I))
            end
            m = d + clouds(40, 104)
            @test_throws ErrorException BLR.raytrace!(m; backend=backend, T=Float32)
        end

        @testset "profile + variance kernels (Float32 atomics) on MetalBackend" begin
            m = builders[1]()
            ma = BLR.gpu(m; backend=backend).ma
            ma32 = BLR.flatten(m; T=Float32)
            edges = collect(range(-0.08, 0.08, length=41))
            lp = Metal.zeros(Float32, length(edges) - 1)
            BLR._rt_line_profile!(lp, ma, MtlArray(Float32.(edges)); overflow=true, backend=backend)
            lpC = zeros(Float32, length(edges) - 1)
            BLR._rt_line_profile!(lpC, ma32, Float32.(edges); overflow=true, backend=cpuB)
            @test relclose(Array(lp), lpC; rtol=1e-5)
            refLP = BLR.getProfile(m, :line, bins=edges, centered=false, overflow=true).binSums
            @test relclose(Array(lp), refLP; rtol=1e-4)

            θ = ma.α .* 1f-3
            w = ma.I .* ma.ΔA
            nb = length(edges) - 1
            σ² = Metal.zeros(Float32, nb); sumW = similar(σ²); sumWθ = similar(σ²)
            μ = similar(σ²); sumWδ² = similar(σ²)
            BLR._rt_weighted_variance!(σ², sumW, sumWθ, μ, sumWδ², ma.v, θ, w, MtlArray(Float32.(edges));
                overflow=true, backend=backend)
            refσ² = BLR.binnedVariance(Float64.(Array(ma.v)), Float64.(Array(θ)), Float64.(Array(w)),
                bins=edges, overflow=true)[3]
            @test relclose(Array(σ²), refσ²; rtol=1e-3)
        end

        @testset "resident-model observables on Metal ≈ CPU" begin
            m = BLR.raytrace!(builders[1]())
            gm = BLR.gpu(m; backend=backend)
            cm = BLR.resident(m; T=Float32)                 # CPU backend, same Float32 columns
            U = [40.0, -12.0]; V = [5.0, 33.0]; PA = 0.6; BLRAng = 1e-11
            bins = 60
            for sym in (:line, :r, :ϕ, :delay)
                a = BLR.getProfile(gm, sym; bins=bins).binSums
                @test relclose(a, BLR.getProfile(cm, sym; bins=bins).binSums; rtol=1e-5)
                @test relclose(a, BLR.getProfile(m, sym; bins=bins).binSums; rtol=2e-3)
            end
            for sym in (:phase, :moment2)
                a = BLR.getProfile(gm, sym; bins=bins, U=U, V=V, PA=PA, BLRAng=BLRAng).binSums
                @test relclose(a, BLR.getProfile(cm, sym; bins=bins, U=U, V=V, PA=PA, BLRAng=BLRAng).binSums; rtol=1e-4)
            end
            fv = filter(isfinite, BLR.getVariable(m, :v, flatten=true))
            fd = filter(isfinite, BLR.getVariable(m, BLR.t, flatten=true))
            vEdges = collect(range(minimum(fv), maximum(fv), length=16))
            tEdges = collect(range(minimum(fd), maximum(fd), length=21))
            for overflow in (false, true)
                @test relclose(BLR.getΨt(gm, tEdges, overflow), BLR.getΨt(cm, tEdges, overflow); rtol=1e-5)
            end
            @test relclose(vec(BLR.getΨ(gm, vEdges, tEdges)), vec(BLR.getΨ(cm, vEdges, tEdges)); rtol=1e-5)
            @test isapprox(BLR.secondMoment(gm; U=U, V=V, PA=PA, BLRAng=BLRAng, returnAvg=true, bins=bins)[4],
                BLR.secondMoment(m; U=U, V=V, PA=PA, BLRAng=BLRAng, returnAvg=true, bins=bins)[4]; rtol=1e-3)
        end

        @testset "transfer function kernels on MetalBackend" begin
            m = BLR.raytrace!(builders[1]())
            ma = BLR.gpu(m; backend=backend).ma
            delays = BLR._rt_transfer_delays(ma; backend=backend)
            @test approx_eq(Array(delays), Array(ma.η) .* (Array(ma.r) .- Array(ma.x)); rtol=1e-6, atol=1e-6)

            host_delays = Array(delays)
            host_v = Array(ma.v); host_I = Array(ma.I); host_A = Array(ma.ΔA)
            finV = filter(isfinite, host_v); finT = filter(isfinite, host_delays)
            vEdges = Float32.(collect(range(minimum(finV), maximum(finV), length=21)))
            tEdges = Float32.(collect(range(minimum(finT), maximum(finT), length=16)))

            Ψ = Metal.zeros(Float32, length(vEdges)-1, length(tEdges)-1)
            BLR._rt_psi2d!(Ψ, ma.v, delays, ma.I, ma.ΔA, MtlArray(vEdges), MtlArray(tEdges); backend=backend)
            ΨC = zeros(Float32, length(vEdges)-1, length(tEdges)-1)
            BLR._rt_psi2d!(ΨC, host_v, host_delays, host_I, host_A, vEdges, tEdges; backend=cpuB)
            @test relclose(vec(Array(Ψ)), vec(ΨC); rtol=1e-5)

            Ψt = Metal.zeros(Float32, length(tEdges)-1); under = Metal.zeros(Float32, 1); over = Metal.zeros(Float32, 1)
            BLR._rt_psit!(Ψt, under, over, delays, ma.I, ma.ΔA, MtlArray(tEdges); backend=backend)
            ΨtC = zeros(Float32, length(tEdges)-1); uC = zeros(Float32, 1); oC = zeros(Float32, 1)
            BLR._rt_psit!(ΨtC, uC, oC, host_delays, host_I, host_A, tEdges; backend=cpuB)
            @test relclose(Array(Ψt), ΨtC; rtol=1e-5)
            @test isapprox(Array(under)[1], uC[1]; rtol=1e-5, atol=1e-7)
            @test isapprox(Array(over)[1], oC[1]; rtol=1e-5, atol=1e-7)
        end

        @testset "disk deprojection / velocity / intensity kernels on MetalBackend" begin
            inc = 0.4; rMin = 300.0; rMax = 900.0
            m = BLR.DiskWindModel(rMin, rMax, inc, nr=16, nϕ=32, scale=:linear,
                I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.7, f3=0.2, f4=0.9,
                α=1.2, ηₒ=0.4, η₁=0.6, αRM=0.1, rNorm=700.0, reflect=false)
            r3D = BLR.get_r3D(inc, 0.0, 0.0)
            undoTilt = [sin(inc) 0.0 -cos(inc); 0.0 1.0 0.0; cos(inc) 0.0 sin(inc)]
            M = undoTilt * r3D
            args = (inc, 0.0, 0.0, M, r3D, rMin, rMax, 0.4, 0.6, 0.1, 700.0)   # Float64 scalars: the wrapper converts
            αh = Float32.(vec(m.camera.α)); βh = Float32.(vec(m.camera.β))
            cpuOuts = ntuple(_ -> similar(αh), 7)
            BLR._rt_disk_deproject!(cpuOuts..., αh, βh, args...; backend=cpuB)
            α = MtlArray(αh); β = MtlArray(βh)
            gpuOuts = ntuple(_ -> similar(α), 7)
            BLR._rt_disk_deproject!(gpuOuts..., α, β, args...; backend=backend)
            rc = cpuOuts[1]
            isbnd(rv) = isfinite(rv) && (abs(rv - rMin) / rMin < 1e-5 || abs(rv - rMax) / rMax < 1e-5)
            for k in 1:7
                a = Array(gpuOuts[k]); b = cpuOuts[k]
                for j in eachindex(a, b)
                    if isnan(a[j]) != isnan(b[j])
                        @test isbnd(rc[j])
                    elseif isfinite(b[j])
                        @test isapprox(a[j], b[j]; rtol=1e-4, atol=1e-4)
                    end
                end
            end
            v = similar(gpuOuts[1]); vC = similar(cpuOuts[1])
            BLR._rt_v_circular_disk!(v, gpuOuts[1], gpuOuts[2], inc; backend=backend)
            BLR._rt_v_circular_disk!(vC, cpuOuts[1], cpuOuts[2], inc; backend=cpuB)
            I = similar(gpuOuts[1]); IC = similar(cpuOuts[1])
            BLR._rt_disk_wind_i!(I, gpuOuts[1], gpuOuts[2], inc, 1.0, 0.7, 0.2, 0.9, 1.2, rMin, rMax; backend=backend)
            BLR._rt_disk_wind_i!(IC, cpuOuts[1], cpuOuts[2], inc, 1.0, 0.7, 0.2, 0.9, 1.2, rMin, rMax; backend=cpuB)
            fin = isfinite.(Array(I)) .& isfinite.(IC)
            @test count(fin) >= count(isfinite, IC) - 2       # NaN (out-of-range) masks agree up to boundary pixels
            @test maximum(abs.(Array(v)[fin] .- vC[fin])) < 1e-5 * maximum(abs, vC[fin])
            @test maximum(abs.(Array(I)[fin] .- IC[fin])) < 1e-4 * maximum(abs, IC[fin])
        end

        @testset "on-device DiskWind construction on MetalBackend" begin
            fkw = (f1=1.0, f2=0.7, f3=0.2, f4=0.9, α=1.2, ηₒ=0.4, η₁=0.6, αRM=0.1, rNorm=700.0)
            rMin, rMax, inc, nr, nϕ = 311.7, 887.3, 0.4, 32, 64
            cols = (:r, :ϕ, :ϕ₀, :v, :I, :ΔA, :η, :x, :y, :z, :α, :β)
            isbnd(rv) = isfinite(rv) && (abs(rv - rMin) / rMin < 1e-5 || abs(rv - rMax) / rMax < 1e-5)
            for scale in (:linear, :log)
                rmCPU = BLR.residentDiskWindModel(rMin, rMax, inc; nr=nr, nϕ=nϕ, scale=scale,
                    backend=cpuB, T=Float32, fkw...)
                rmGPU = BLR.gpuDiskWindModel(rMin, rMax, inc; nr=nr, nϕ=nϕ, scale=scale,
                    backend=backend, fkw...)
                @test rmGPU.ma.r isa Metal.MtlVector{Float32}
                rc = rmCPU.ma.r
                for c in cols
                    a = Array(getfield(rmGPU.ma, c)); b = getfield(rmCPU.ma, c)
                    ok = true
                    for k in eachindex(a, b)
                        if isnan(a[k]) != isnan(b[k])
                            ok &= isbnd(rc[k])
                        elseif isfinite(b[k])
                            ok &= isapprox(a[k], b[k]; rtol=1e-4, atol=1e-4 * maximum(abs, filter(isfinite, b)))
                        end
                    end
                    @test ok
                end
                # line profile: tight vs the CPU Float32 build, loose vs Float64 (a Float32 build can
                # include/exclude one boundary ring -- the documented B2 masking effect -- which moves a
                # log-grid profile by ~1%). Edges avoid v = 0, where a symmetric disk puts whole
                # columns of points exactly on an edge.
                rm64 = BLR.residentDiskWindModel(rMin, rMax, inc; nr=nr, nϕ=nϕ, scale=scale, backend=cpuB, T=Float64, fkw...)
                e = collect(range(-0.08, 0.08, length=42))
                pG = BLR.getProfile(rmGPU, :line; bins=e).binSums
                @test relclose(pG, BLR.getProfile(rmCPU, :line; bins=e).binSums; rtol=1e-5)
                @test relclose(pG, BLR.getProfile(rm64, :line; bins=e).binSums; rtol=2e-2)
            end
            # r̄/rFac/α form + default backend
            rmD = BLR.gpuDiskWindModel(3000.0, 100.0, 1.0, 75 / 180 * π; nr=64, nϕ=64, f1=1.0, f2=1.0, f3=1.0, f4=1.0,
                backend=backend)
            @test rmD isa BLR.ResidentModel && eltype(rmD.ma.I) == Float32
        end

        @testset "on-device cloud construction on MetalBackend" begin
            cp = (μ=600.0, β=1.0, F=0.5, rₛ=1.0, θₒ=0.5, γ=1.0, ξ=0.8, i=0.4,
                  ηₒ=0.4, η₁=0.6, αRM=0.1, rNorm=700.0)
            N, seed = 20000, 4242
            rmGPU = BLR.gpuCloudModel(N, seed; backend=backend, cp...)
            rmCPU = BLR.residentCloudModel(N, seed; backend=cpuB, T=Float32, cp...)
            @test rmGPU.ma.r isa Metal.MtlVector{Float32}
            # same Philox substreams -> the pure-uniform draws are identical, and the full per-cloud
            # realization agrees with the CPU backend up to device-libm ULPs (a rejection-sampler
            # accept/reject or reflection coin landing exactly on a ULP boundary may flip a tiny fraction)
            @test Array(rmGPU.ma.ϕ₀) == rmCPU.ma.ϕ₀
            @test Array(rmGPU.ma.rot) == rmCPU.ma.rot
            rG = Array(rmGPU.ma.r); rC = rmCPU.ma.r
            @test count(isapprox.(rG, rC; rtol=1e-4)) >= N - 5
            @test count(Array(rmGPU.ma.reflect) .== rmCPU.ma.reflect) >= N - 5
            vG = Array(rmGPU.ma.v); vC = rmCPU.ma.v
            @test count(isapprox.(vG, vC; rtol=1e-3, atol=1e-6)) >= N - 5
            x = Array(rmGPU.ma.x); y = Array(rmGPU.ma.y); z = Array(rmGPU.ma.z)
            @test maximum(abs.(sqrt.(x .^ 2 .+ y .^ 2 .+ z .^ 2) .- rG) ./ rG) < 1e-5
            ks2(a, b) = (A = sort(a); B = sort(b); na = length(A); nb = length(B);
                maximum(abs(searchsortedlast(A, v) / na - searchsortedlast(B, v) / nb) for v in vcat(A, B)))
            mh = BLR.cloudModel(N; μ=cp.μ, β=cp.β, F=cp.F, θₒ=cp.θₒ, γ=cp.γ, ξ=cp.ξ, i=cp.i, rₛ=cp.rₛ,
                I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, rng=:philox, seed=seed)
            @test ks2(Float64.(rG), BLR.getVariable(mh, :r, flatten=true)) < 1.358 * sqrt(2 / N)

            tp = (σρᵣ=0.2, σρc=0.04, σΘᵣ=0.4, σΘc=0.1, θₑ=35 / 180 * π, fEllipse=0.8, fFlow=0.0, σₜ=0.05)
            kw = (μ=cp.μ, β=cp.β, F=cp.F, rₛ=cp.rₛ, θₒ=cp.θₒ, γ=cp.γ, ξ=cp.ξ, i=cp.i, intensity=:cloud, κ=0.4,
                velocity=:turbulent, tp...)
            rmT = BLR.gpuCloudModel(N, seed; backend=backend, kw...)
            rmTC = BLR.residentCloudModel(N, seed; backend=cpuB, T=Float32, kw...)
            vT = Array(rmT.ma.v)
            @test all(isfinite, vT)
            @test count(isapprox.(vT, rmTC.ma.v; rtol=1e-3, atol=1e-6)) >= N - 10
            @test count(isapprox.(Array(rmT.ma.I), rmTC.ma.I; rtol=1e-4, atol=1e-6)) >= N - 10
            mhT = BLR.cloudModel(N; μ=cp.μ, β=cp.β, F=cp.F, rₛ=cp.rₛ, θₒ=cp.θₒ, γ=cp.γ, ξ=cp.ξ, i=cp.i,
                κ=0.4, I=BLR.cloudIntensity, v=BLR.vCloudTurbulentEllipticalFlow, tp..., rng=:philox, seed=seed)
            @test ks2(Float64.(vT), BLR.getVariable(mhT, :v, flatten=true)) < 1.358 * sqrt(2 / N)
        end

        @testset "vCircularRadialDisk on MetalBackend" begin
            fkw = (f1=0.0, f2=0.0, f3=0.0, f4=1.0, α=1.0, ηₒ=0.5, η₁=0.5, αRM=0.0, rNorm=1.0)
            rMin, rMax, inc = 311.7, 887.3, 0.5
            rmCPU = BLR.residentDiskWindModel(rMin, rMax, inc; nr=24, nϕ=48, scale=:log,
                vᵣFrac=0.33, inflow=true, backend=cpuB, T=Float32, fkw...)
            rmGPU = BLR.gpuDiskWindModel(rMin, rMax, inc; nr=24, nϕ=48, scale=:log,
                vᵣFrac=0.33, inflow=true, backend=backend, fkw...)
            vc = rmCPU.ma.v; vg = Array(rmGPU.ma.v)
            fin = isfinite.(vc) .& isfinite.(vg)
            @test maximum(abs.(vc[fin] .- vg[fin])) < 1e-5 * maximum(abs, vc[fin])
        end

        @testset "device-resident raytrace! on MetalBackend" begin
            dk(r1, r2; nr=16, nϕ=32, τ=5.0) = BLR.DiskWindModel(r1, r2, 0.4; nr=nr, nϕ=nϕ, scale=:linear,
                I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.5, f3=0.2, f4=0.3, α=1.0, τ=τ)
            cl(n, seed; μ=600.0, τ=0.1) = BLR.cloudModel(n; μ=μ, β=1.0, F=0.5, θₒ=0.4, i=0.4, γ=1.0,
                ξ=0.8, I=BLR.IsotropicIntensity, v=BLR.vCircularCloud, τ=τ, rng=:philox, seed=seed)
            function check(builder; rfc=false, IR=1.0)
                # CPU backend on the same Float32 columns = the tight reference; host Float64 = sanity
                cref = BLR.raytrace!(BLR.resident(builder(); T=Float32, raytrace=true); IRatios=IR, raytraceFreeClouds=rfc)
                href = BLR.resident(BLR.raytrace!(builder(); IRatios=IR, raytraceFreeClouds=rfc))
                rmg = BLR.gpu(builder(); backend=backend)
                @test rmg.rt isa BLR.RaytraceMeta
                rrt = BLR.raytrace!(rmg; IRatios=IR, raytraceFreeClouds=rfc)
                @test rrt.ma.r isa Metal.MtlVector{Float32}
                @test length(rrt.ma.I) == length(cref.ma.I)
                fC = sum(filter(isfinite, Float64.(cref.ma.I .* cref.ma.ΔA)))
                fD = sum(filter(isfinite, Float64.(Array(rrt.ma.I) .* Array(rrt.ma.ΔA))))
                fH = sum(filter(isfinite, href.ma.I .* href.ma.ΔA))
                @test isapprox(fD, fC; rtol=1e-5)
                @test isapprox(fD, fH; rtol=1e-4)
                e = collect(range(-0.08, 0.08, length=41))
                @test relclose(BLR.getProfile(rrt, :line; bins=e).binSums, BLR.getProfile(cref, :line; bins=e).binSums; rtol=1e-4)
            end
            check(() -> dk(300., 900.) + cl(300, 1))
            check(() -> dk(250., 700.) + dk(500., 1000.))
            check(() -> dk(300., 900.) + cl(300, 4); IR=[1.0, 0.25])
            check(() -> cl(150, 1, μ=300., τ=2.0) + cl(150, 2, μ=320., τ=2.0); rfc=true)
        end

        @testset "on-device-built model raytrace! (Phase 2) on MetalBackend" begin
            rmG = BLR.gpuDiskWindModel(300.0, 900.0, 0.4; nr=16, nϕ=32, scale=:linear,
                f1=1.0, f2=0.5, f3=0.2, f4=0.3, α=1.0, backend=backend) +
                BLR.gpuCloudModel(400, 7; μ=600.0, β=1.0, F=0.5, θₒ=0.4, i=0.4, γ=1.0, ξ=0.8, backend=backend)
            @test rmG.rt isa BLR.RaytraceMeta && rmG.ma.r isa Metal.MtlVector
            resG = BLR.raytrace!(rmG)
            @test resG.ma.r isa Metal.MtlVector
            resC = BLR.raytrace!(BLR.cpu(rmG))
            @test length(resG.ma.I) == length(resC.ma.I)
            @test isapprox(sort(Array(resG.ma.I)), sort(resC.ma.I); rtol=1e-5)
            e = collect(range(-0.08, 0.08, length=41))
            @test relclose(BLR.getProfile(resG, :line; bins=e).binSums, BLR.getProfile(resC, :line; bins=e).binSums; rtol=1e-4)
        end

        @testset "ResidentModel device combine (+) on MetalBackend" begin
            d = BLR.gpuDiskWindModel(311.7, 887.3, 0.4; nr=16, nϕ=32, scale=:log,
                f1=1.0, f2=0.7, f3=0.2, f4=0.9, α=1.2, ηₒ=0.4, η₁=0.6, αRM=0.1, rNorm=700.0, backend=backend)
            c = BLR.gpuCloudModel(4000, 7; μ=600.0, β=1.0, F=0.5, θₒ=0.4, i=0.4, γ=1.0, ξ=0.8, backend=backend)
            s = d + c
            @test s.ma.r isa Metal.MtlVector
            @test s.nSubModels == 2
            @test length(s.ma.I) == length(d.ma.I) + length(c.ma.I)
            @test count(isfinite, BLR.getProfile(s, :line; bins=40).binSums) > 0
            cpuRm = BLR.cpu(c)
            @test cpuRm.ma.r isa Array
            err = try
                d + cpuRm
            catch e
                e
            end
            @test err isa ArgumentError && occursin("different backends", err.msg)
            @test (BLR.cpu(d) + cpuRm) isa BLR.ResidentModel
        end

        @testset "resident composite model on MetalBackend" begin
            dm = disk()
            cl = clouds(80, 105)
            cm = BLR.CompositeModel(dm; line="Ha", lineCenter=6562.8)
            BLR.addLine!(cm, cl; line="Hb", lineCenter=4861.3, fluxRatio=0.35)
            rcm = BLR.gpu(cm; backend=backend)
            rcmC = BLR.resident(cm; T=Float32)
            @test rcm isa BLR.ResidentCompositeModel
            @test rcm["Ha"].ma.I isa Metal.MtlVector{Float32}
            wC = BLR._fluxWeights(cm); wD = BLR._fluxWeights(rcm)
            for line in cm.lines
                @test isapprox(wD[line], wC[line]; rtol=1e-4)
                vminC, vmaxC = BLR._finiteVRange(rcmC[line])
                vminD, vmaxD = BLR._finiteVRange(rcm[line])
                @test vminD == vminC && vmaxD == vmaxC      # min/max reductions are exact
            end
            _, _, fC, tC = BLR.getSpectrum(rcmC; bins=48)
            _, _, fD, tD = BLR.getSpectrum(rcm; bins=48)
            for line in cm.lines
                @test relclose(fD[line], fC[line]; rtol=1e-4)
                @test isapprox(sum(fD[line]), cm.fluxRatios[line]; rtol=1e-4)
            end
            @test relclose(tD, tC; rtol=1e-4)
            vEdges = collect(range(-0.09, 0.09, length=37))
            pC = BLR.getProfile(rcmC, :ratio; lines=("Ha", "Hb"), bins=vEdges)
            pD = BLR.getProfile(rcm, :ratio; lines=("Ha", "Hb"), bins=vEdges)
            @test isnan.(pD.binSums) == isnan.(pC.binSums)
            ok = .!isnan.(pC.binSums)
            @test approx_eq(pD.binSums[ok], pC.binSums[ok]; rtol=1e-4, atol=1e-6)
            @test BLR.lineRatio(rcm, "Ha", "Hb") == BLR.lineRatio(cm, "Ha", "Hb")
        end
    end
end
