# GPU correctness tests. These only run when BLR_TEST_CUDA=1 AND a functional CUDA device is
# present (see the include guard in runtests.jl). Everything here is validated against the CPU
# reference path so the GPU is held to the same numbers, not just "does it run".
using CUDA
using KernelAbstractions

@testset "CUDA GPU raytrace + kernels" begin
    if !CUDA.functional()
        @warn "BLR_TEST_CUDA=1 but CUDA.functional() is false -- skipping GPU tests"
    else
        @info "running GPU tests on $(CUDA.name(CUDA.device()))"
        backend = CUDA.CUDABackend()

        # NaN-aware elementwise approx (raytrace output is compacted so it should be NaN-free,
        # but stay defensive).
        function approx_eq(a, b; rtol=1e-10, atol=1e-12)
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

        @testset "flatten/gpu round-trip" begin
            m = builders[1]()
            ma_cpu = BLR.flatten(m)
            gm = BLR.gpu(m; T=Float64)
            @test gm isa BLR.ResidentModel
            @test gm.ma.I isa CUDA.CuVector{Float64}
            @test gm.ma.reflect isa CUDA.CuVector{Bool}
            back = BLR.cpu(gm).ma
            @test isequal(back.I, ma_cpu.I)
            @test isequal(back.x, ma_cpu.x)
            @test back.α == ma_cpu.α
        end

        @testset "extension sort matches CPU reference" begin
            keys = [2, 1, 2, 1, 0, 2, 0, 2, 1]
            x = [10.0, 4.0, 11.0, 6.0, 3.0, NaN, 9.0, 11.0, 4.0]  # includes a depth tie (two 11.0, two 4.0)
            permCPU = BLR._rt_sortperm_by_key_depth(keys, x)
            permGPU = Array(BLR._rt_sortperm_by_key_depth(CuArray(keys), CuArray(x)))
            # The scan only depends on the (key, depth) sequence; assert that bit-matches, and that
            # the tie-free composite reproduces the CPU permutation exactly.
            @test BLR._rt_sorted_key_depth_pairs(keys, x, permGPU) ==
                  BLR._rt_sorted_key_depth_pairs(keys, x, permCPU)
            @test permGPU == permCPU
        end

        @testset "GPU raytrace! == CPU (Float64)" begin
            for build in builders
                ref = BLR.raytrace!(build(); τCutOff=1.0)
                gpu = BLR.raytrace!(build(); τCutOff=1.0, backend=backend, T=Float64)
                for sym in (:I, :v, :r, :ϕ, :ϕ₀, :τ, :η)
                    @test approx_eq(BLR.getVariable(gpu, sym, flatten=true),
                        BLR.getVariable(ref, sym, flatten=true))
                end
                @test gpu.subModelStartInds == ref.subModelStartInds
                @test approx_eq(gpu.camera.α, ref.camera.α)
                @test approx_eq(gpu.camera.β, ref.camera.β)
                # integrated invariants (robust to any last-bit reduction differences)
                refI = BLR.getVariable(ref, :I, flatten=true); refA = BLR.getVariable(ref, :ΔA, flatten=true)
                gpuI = BLR.getVariable(gpu, :I, flatten=true); gpuA = BLR.getVariable(gpu, :ΔA, flatten=true)
                @test isapprox(sum(gpuI .* gpuA), sum(refI .* refA); rtol=1e-10)
                edges = collect(range(-0.08, 0.08, length=41))
                @test approx_eq(BLR.getProfile(gpu, :line, bins=edges, centered=false).binSums,
                    BLR.getProfile(ref, :line, bins=edges, centered=false).binSums; rtol=1e-9, atol=1e-12)
            end
        end

        @testset "GPU raytrace! ≈ CPU (Float32) + IRatios" begin
            ref = BLR.raytrace!(builders[1](); IRatios=[1.0, 0.25], τCutOff=1.0)
            gpu = BLR.raytrace!(builders[1](); IRatios=[1.0, 0.25], τCutOff=1.0, backend=backend, T=Float32)
            for sym in (:I, :v, :r)
                @test approx_eq(BLR.getVariable(gpu, sym, flatten=true),
                    BLR.getVariable(ref, sym, flatten=true); rtol=1e-4, atol=1e-6)
            end
        end

        @testset "velocity-dependent τ errors on GPU path" begin
            d = disk()
            for r in d.rings
                r.τ = fill(0.4, length(r.I))   # per-point τ vector => velocity-dependent
            end
            m = d + clouds(40, 104)
            @test_throws ErrorException BLR.raytrace!(m; backend=backend, T=Float64)
        end

        @testset "profile + variance kernels on CUDABackend" begin
            m = builders[1]()
            ma = BLR.gpu(m; T=Float64).ma
            edges = collect(range(-0.08, 0.08, length=41))
            dEdges = CuArray(edges)
            lp = CUDA.zeros(Float64, length(edges) - 1)
            BLR._rt_line_profile!(lp, ma, dEdges; overflow=true, backend=backend)
            refLP = BLR.getProfile(m, :line, bins=edges, centered=false, overflow=true).binSums
            @test approx_eq(Array(lp), refLP; rtol=1e-9, atol=1e-12)

            θ = ma.α .* 1e-10
            w = ma.I .* ma.ΔA
            nb = length(edges) - 1
            σ² = CUDA.zeros(Float64, nb); sumW = similar(σ²); sumWθ = similar(σ²)
            μ = similar(σ²); sumWδ² = similar(σ²)
            BLR._rt_weighted_variance!(σ², sumW, sumWθ, μ, sumWδ², ma.v, θ, w, dEdges;
                overflow=true, backend=backend)
            refσ² = BLR.binnedVariance(Array(ma.v), Array(θ), Array(w), bins=edges, overflow=true)[3]
            @test approx_eq(Array(σ²), refσ²; rtol=1e-9, atol=1e-20)
        end

        @testset "resident-model observables on GPU == CPU" begin
            m = BLR.raytrace!(builders[1]())
            gm = BLR.gpu(m; T=Float64)
            @test gm isa BLR.ResidentModel
            U = [40.0, -12.0]; V = [5.0, 33.0]; PA = 0.6; BLRAng = 1e-11
            bins = 60
            for sym in (:line, :r, :ϕ)
                @test approx_eq(BLR.getProfile(gm, sym; bins=bins).binSums,
                    BLR.getProfile(m, sym; bins=bins).binSums; rtol=1e-9, atol=1e-12)
            end
            @test approx_eq(BLR.getProfile(gm, :delay; bins=bins).binSums,
                BLR.getProfile(m, :delay; bins=bins).binSums; rtol=1e-8, atol=1e-8)
            for sym in (:phase, :moment2)
                @test approx_eq(BLR.getProfile(gm, sym; bins=bins, U=U, V=V, PA=PA, BLRAng=BLRAng).binSums,
                    BLR.getProfile(m, sym; bins=bins, U=U, V=V, PA=PA, BLRAng=BLRAng).binSums; rtol=1e-8, atol=1e-10)
            end
            fv = filter(isfinite, BLR.getVariable(m, :v, flatten=true))
            fd = filter(isfinite, BLR.getVariable(m, BLR.t, flatten=true))
            vEdges = collect(range(minimum(fv), maximum(fv), length=16))
            tEdges = collect(range(minimum(fd), maximum(fd), length=21))
            for overflow in (false, true)
                @test approx_eq(BLR.getΨt(gm, tEdges, overflow), BLR.getΨt(m, tEdges, overflow); rtol=1e-8, atol=1e-12)
            end
            @test approx_eq(vec(BLR.getΨ(gm, vEdges, tEdges)), vec(BLR.getΨ(m, vEdges, tEdges)); rtol=1e-8, atol=1e-12)
            @test isapprox(BLR.secondMoment(gm; U=U, V=V, PA=PA, BLRAng=BLRAng, returnAvg=true, bins=bins)[4],
                BLR.secondMoment(m; U=U, V=V, PA=PA, BLRAng=BLRAng, returnAvg=true, bins=bins)[4]; rtol=1e-8)
        end

        @testset "transfer function kernels on CUDABackend" begin
            m = BLR.raytrace!(builders[1]())
            ma = BLR.gpu(m; T=Float64).ma
            delays = BLR._rt_transfer_delays(ma; backend=backend)
            @test approx_eq(Array(delays), Array(ma.η) .* (Array(ma.r) .- Array(ma.x)); rtol=1e-12, atol=1e-12)

            host_delays = Array(delays)
            host_v = Array(ma.v); host_I = Array(ma.I); host_A = Array(ma.ΔA)
            finV = filter(isfinite, host_v); finT = filter(isfinite, host_delays)
            vEdges = collect(range(minimum(finV), maximum(finV), length=21))
            tEdges = collect(range(minimum(finT), maximum(finT), length=16))

            Ψ = CUDA.zeros(Float64, length(vEdges)-1, length(tEdges)-1)
            BLR._rt_psi2d!(Ψ, ma.v, delays, ma.I, ma.ΔA, CuArray(vEdges), CuArray(tEdges); backend=backend)
            Ψref = zeros(length(vEdges)-1, length(tEdges)-1)
            BLR.ΨAccumulate!(Ψref, host_v, host_delays, host_I, host_A, vEdges, tEdges)
            @test approx_eq(vec(Array(Ψ)), vec(Ψref); rtol=1e-9, atol=1e-12)

            Ψt = CUDA.zeros(Float64, length(tEdges)-1); under = CUDA.zeros(Float64, 1); over = CUDA.zeros(Float64, 1)
            BLR._rt_psit!(Ψt, under, over, delays, ma.I, ma.ΔA, CuArray(tEdges); backend=backend)
            Ψtref = zeros(length(tEdges)-1)
            sU, sO = BLR.ΨtAccumulate!(Ψtref, host_delays, host_I, host_A, tEdges)
            @test approx_eq(Array(Ψt), Ψtref; rtol=1e-9, atol=1e-12)
            @test isapprox(Array(under)[1], sU; rtol=1e-9, atol=1e-12)
            @test isapprox(Array(over)[1], sO; rtol=1e-9, atol=1e-12)
        end

        @testset "disk deprojection kernel on CUDABackend" begin
            inc = 0.4; rMin = 300.0; rMax = 900.0
            m = BLR.DiskWindModel(rMin, rMax, inc, nr=16, nϕ=32, scale=:linear,
                I=BLR.DiskWindIntensity, v=BLR.vCircularDisk, f1=1.0, f2=0.7, f3=0.2, f4=0.9,
                α=1.2, ηₒ=0.4, η₁=0.6, αRM=0.1, rNorm=700.0, reflect=false)
            r3D = BLR.get_r3D(inc, 0.0, 0.0)
            undoTilt = [sin(inc) 0.0 -cos(inc); 0.0 1.0 0.0; cos(inc) 0.0 sin(inc)]
            M = undoTilt * r3D
            args = (inc, 0.0, 0.0, M, r3D, rMin, rMax, 0.4, 0.6, 0.1, 700.0)
            # CPU-backend run of this kernel is validated against the model in the CPU testset;
            # here we just confirm the CUDABackend run reproduces it bit-for-bit.
            cpuOuts = ntuple(_ -> similar(m.camera.α), 7)
            BLR._rt_disk_deproject!(cpuOuts..., m.camera.α, m.camera.β, args...;
                backend=KernelAbstractions.CPU())
            α = CuArray(m.camera.α); β = CuArray(m.camera.β)
            gpuOuts = ntuple(_ -> similar(α), 7)
            BLR._rt_disk_deproject!(gpuOuts..., α, β, args...; backend=backend)
            for k in 1:7
                @test approx_eq(vec(Array(gpuOuts[k])), vec(cpuOuts[k]); rtol=1e-10, atol=1e-10)
            end
        end
    end
end
