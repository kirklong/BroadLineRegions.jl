#!/usr/bin/env julia
# W2-B-rt5: CPU Phase A vs CUDA backend raytrace benchmark.
#
# Quick smoke:
#   julia --project=. benchmark/gpu_bench.jl
#
# Full plan-sized run, including separate cloud, disk, and both-scaling regimes:
#   FULL=1 REPEATS=3 julia --project=. benchmark/gpu_bench.jl
#
# CUDA is a weak dependency. Run this in an environment where CUDA.jl is
# available, for example a temporary environment that develops this package and
# adds CUDA.
using BroadLineRegions
using Printf
using Random

try
    @eval using Adapt
    @eval using CUDA
catch err
    @warn "CUDA.jl and Adapt.jl are required for benchmark/gpu_bench.jl" exception=(err, catch_backtrace())
    exit(1)
end

const QUICK_CONFIGS = (
    (regime="quick", nr=16, nphi=32, ncloud=1_000),
)

const CLOUD_SCALING_CONFIGS = (
    (regime="cloud64", nr=64, nphi=128, ncloud=1_000),
    (regime="cloud64", nr=64, nphi=128, ncloud=2_500),
    (regime="cloud64", nr=64, nphi=128, ncloud=5_000),
    (regime="cloud64", nr=64, nphi=128, ncloud=20_000),
    (regime="cloud64", nr=64, nphi=128, ncloud=100_000),
    (regime="cloud64", nr=64, nphi=128, ncloud=1_000_000),
)

const DISK_SCALING_CONFIGS = (
    (regime="disk1k", nr=64, nphi=128, ncloud=1_000),
    (regime="disk1k", nr=128, nphi=256, ncloud=1_000),
    (regime="disk1k", nr=256, nphi=512, ncloud=1_000),
    (regime="disk1k", nr=512, nphi=1024, ncloud=1_000),
    (regime="disk1k", nr=1024, nphi=2048, ncloud=1_000),
    (regime="disk1k", nr=2048, nphi=4096, ncloud=1_000),
)

const BOTH_SCALING_CONFIGS = (
    (regime="both", nr=64, nphi=128, ncloud=1_000),
    (regime="both", nr=128, nphi=256, ncloud=5_000),
    (regime="both", nr=256, nphi=512, ncloud=20_000),
    (regime="both", nr=512, nphi=1024, ncloud=100_000),
    (regime="both", nr=1024, nphi=2048, ncloud=250_000),
    (regime="both", nr=2048, nphi=4096, ncloud=1_000_000),
)

isfull() = get(ENV, "FULL", "0") == "1"
repeats() = parse(Int, get(ENV, "REPEATS", isfull() ? "3" : "1"))
configs() = isfull() ? (CLOUD_SCALING_CONFIGS..., DISK_SCALING_CONFIGS..., BOTH_SCALING_CONFIGS...) : QUICK_CONFIGS

mkdisk(nr, nphi) = BLR.DiskWindModel(500.0, 5.0, 1.0, 45/180*pi, nr=nr, nϕ=nphi,
    scale=:log, f1=1.0, f2=1.0, f3=1.0, f4=1.0, reflect=false, τ=5.0)

mkcloud(n) = BLR.cloudModel(n, μ=500.0, β=1.0, F=0.5, θₒ=30/180*pi, γ=1.0,
    ξ=1.0, i=45/180*pi, I=BLR.IsotropicIntensity, v=BLR.vCircularCloud,
    rescale=1e-5, τ=0.0, rng=MersenneTwister(11))

mkmodel(nr, nphi, ncloud) = mkdisk(nr, nphi) + mkcloud(ncloud)

function best_seconds(f; n=repeats(), sync_cuda=false)
    best = Inf
    for _ in 1:n
        GC.gc()
        sync_cuda && CUDA.synchronize()
        t0 = time_ns()
        f()
        sync_cuda && CUDA.synchronize()
        best = min(best, (time_ns() - t0) / 1e9)
    end
    return best
end

function best_seconds(setup, f; n=repeats(), sync_cuda=false)
    best = Inf
    for _ in 1:n
        GC.gc()
        arg = setup()
        sync_cuda && CUDA.synchronize()
        t0 = time_ns()
        f(arg)
        sync_cuda && CUDA.synchronize()
        best = min(best, (time_ns() - t0) / 1e9)
    end
    return best
end

function time_cpu(nr, nphi, ncloud)
    return best_seconds(() -> mkmodel(nr, nphi, ncloud), BLR.raytrace!)
end

function time_gpu_public(nr, nphi, ncloud, ::Type{T}) where {T<:Real}
    backend = CUDA.CUDABackend()
    return best_seconds(() -> mkmodel(nr, nphi, ncloud),
        m -> BLR.raytrace!(m; backend=backend, T=T); sync_cuda=true)
end

function prepare_gpu_resident(nr, nphi, ncloud, ::Type{T}) where {T<:Real}
    backend = CUDA.CUDABackend()
    m = mkmodel(nr, nphi, ncloud)
    cam_start_inds = BLR.getFlattenedCameraIndices(m)
    grids, pixels = BLR._rt_build_output(m, cam_start_inds)[1:2]
    ma = BLR.gpu(m; T=T)
    grid_arrays = Adapt.adapt(CUDA.CuArray, BLR._rt_grid_arrays(grids; T=T))

    keys = BLR._rt_bin_assign(ma, grid_arrays; backend=backend)
    perm = BLR._rt_sortperm_by_key_depth(ma, keys)
    CUDA.synchronize()

    iratios = fill(1.0, length(m.subModelStartInds))
    scan = Adapt.adapt(CUDA.CuArray,
        BLR._rt_scan_arrays(m, Array(keys), Array(perm), pixels, iratios; T=T))
    scan_out = Adapt.adapt(CUDA.CuArray, BLR._rt_scan_output(length(pixels); T=T))
    BLR._rt_segmented_scan!(scan_out, ma, keys, perm, scan; τCutOff=T(1.0), backend=backend)
    CUDA.synchronize()
    return (; backend, ma, grid_arrays, scan, scan_out)
end

function time_gpu_resident(nr, nphi, ncloud, ::Type{T}) where {T<:Real}
    prep = prepare_gpu_resident(nr, nphi, ncloud, T)
    return best_seconds(sync_cuda=true) do
        keys = BLR._rt_bin_assign(prep.ma, prep.grid_arrays; backend=prep.backend)
        perm = BLR._rt_sortperm_by_key_depth(prep.ma, keys)
        BLR._rt_segmented_scan!(prep.scan_out, prep.ma, keys, perm, prep.scan;
            τCutOff=T(1.0), backend=prep.backend)
    end
end

fmt(x) = isfinite(x) ? @sprintf("%9.3f", x) : "       --"

function main()
    CUDA.functional() || error("CUDA.jl loaded but CUDA.functional() is false")
    println("CUDA device: ", CUDA.name(CUDA.device()))
    println("repeats: ", repeats(), isfull() ? " (FULL=1)" : " (quick smoke; set FULL=1 for plan-sized run)")
    println()
    println("resident = preloaded ModelArrays/grid/scan metadata; times GPU bin-sort-scan only")
    println("public   = BLR.raytrace!(...; backend=CUDA.CUDABackend()) including copy + rebuild")
    println()
    @printf("%8s %10s %10s %12s %12s %11s %11s %11s %11s %11s\n",
        "regime", "disk", "clouds", "points", "work", "CPU public", "GPU64 pub", "GPU64 res", "GPU32 pub", "GPU32 res")
    @printf("%8s %10s %10s %12s %12s %11s %11s %11s %11s %11s\n",
        "-"^8, "-"^10, "-"^10, "-"^12, "-"^12, "-"^11, "-"^11, "-"^11, "-"^11, "-"^11)

    # Compile all code paths before reporting the table.
    BLR.raytrace!(mkmodel(4, 8, 16))
    BLR.raytrace!(mkmodel(4, 8, 16); backend=CUDA.CUDABackend(), T=Float64)
    BLR.raytrace!(mkmodel(4, 8, 16); backend=CUDA.CUDABackend(), T=Float32)
    time_gpu_resident(4, 8, 16, Float64)
    time_gpu_resident(4, 8, 16, Float32)

    for cfg in configs()
        nr = cfg.nr
        nphi = cfg.nphi
        ncloud = cfg.ncloud
        points = nr * nphi + ncloud
        work = nr * nphi * ncloud
        cpu = time_cpu(nr, nphi, ncloud)
        gpu64_public = time_gpu_public(nr, nphi, ncloud, Float64)
        gpu64_resident = time_gpu_resident(nr, nphi, ncloud, Float64)
        gpu32_public = time_gpu_public(nr, nphi, ncloud, Float32)
        gpu32_resident = time_gpu_resident(nr, nphi, ncloud, Float32)
        @printf("%8s %4dx%-5d %10d %12d %12.2e %11s %11s %11s %11s %11s\n",
            cfg.regime, nr, nphi, ncloud, points, work, fmt(cpu), fmt(gpu64_public),
            fmt(gpu64_resident), fmt(gpu32_public), fmt(gpu32_resident))
    end
end

main()
