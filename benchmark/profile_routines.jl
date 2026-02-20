#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using SplitVOF
using StaticArrays
using Profile
using Printf
using Base.Threads: nthreads

const CASES = (
    :splitgrid_alloc,
    :initfgrid,
    :youngs_normals,
    :reconstruct,
    :compute_dt,
    :step,
    :integrate,
    :pipeline,
)

function _parse_args(args)
    cfg = Dict{String,Any}(
        "case" => :all,
        "runs" => 50,
        "profile_runs" => 0,
        "dim" => 2,
        "nx" => 48,
        "ny" => 48,
        "nz" => 24,
        "nc" => 6,
        "cfl" => 0.5,
        "dt_factor" => 0.5,
    )

    # Positional fallback: case runs profile_runs dim nx ny nz nc
    if !isempty(args) && !occursin('=', args[1])
        cfg["case"] = Symbol(args[1])
    end
    if length(args) >= 2 && !occursin('=', args[2])
        cfg["runs"] = parse(Int, args[2])
    end
    if length(args) >= 3 && !occursin('=', args[3])
        cfg["profile_runs"] = parse(Int, args[3])
    end
    if length(args) >= 4 && !occursin('=', args[4])
        cfg["dim"] = parse(Int, args[4])
    end
    if length(args) >= 5 && !occursin('=', args[5])
        cfg["nx"] = parse(Int, args[5])
    end
    if length(args) >= 6 && !occursin('=', args[6])
        cfg["ny"] = parse(Int, args[6])
    end
    if length(args) >= 7 && !occursin('=', args[7])
        cfg["nz"] = parse(Int, args[7])
    end
    if length(args) >= 8 && !occursin('=', args[8])
        cfg["nc"] = parse(Int, args[8])
    end

    for a in args
        occursin('=', a) || continue
        k, v = split(a, '='; limit=2)
        if !haskey(cfg, k)
            @warn "Ignoring unknown option" key=k value=v
            continue
        end

        if k == "case"
            cfg[k] = Symbol(v)
        elseif cfg[k] isa Int
            cfg[k] = parse(Int, v)
        else
            cfg[k] = parse(Float64, v)
        end
    end

    return cfg
end

function _params(cfg, dim::Int)
    nx = Int(cfg["nx"])
    ny = Int(cfg["ny"])
    nz = Int(cfg["nz"])
    cfl = Float64(cfg["cfl"])

    if dim == 2
        return SplitVOFParams(
            nx = nx,
            ny = ny,
            nz = 1,
            xlim = (-1.0, 1.0),
            ylim = (-1.0, 1.0),
            zlim = (0.0, 1.0),
            cfl = cfl,
            reconstruction = YoungsPLIC(),
            advection = StrangSplit(),
        )
    elseif dim == 3
        return SplitVOFParams(
            nx = nx,
            ny = ny,
            nz = nz,
            xlim = (-1.0, 1.0),
            ylim = (-1.0, 1.0),
            zlim = (-1.0, 1.0),
            cfl = cfl,
            reconstruction = YoungsPLIC(),
            advection = StrangSplit(),
        )
    else
        throw(ArgumentError("dim must be 2 or 3"))
    end
end

function _shape(dim::Int)
    if dim == 2
        return circle_levelset((0.0, 0.0), 0.35)
    else
        return sphere_levelset((0.0, 0.0, 0.0), 0.45)
    end
end

function _velocity(dim::Int)
    if dim == 2
        return (x, t) -> SVector(-x[2], x[1])
    else
        return (x, t) -> SVector(0.25, -0.10, 0.15)
    end
end

function _prepare_state(cfg, dim::Int)
    nc = Int(cfg["nc"])
    params = _params(cfg, dim)
    vg = vofgrid(params)
    initfgrid!(vg, _shape(dim); nc=nc)
    return vg, params
end

function _prepare_reconstructed_state(cfg, dim::Int)
    vg, params = _prepare_state(cfg, dim)
    reconstruct!(vg, params)
    return vg, params
end

function build_case(case::Symbol; cfg)
    dim = Int(cfg["dim"])
    nc = Int(cfg["nc"])
    dt_factor = Float64(cfg["dt_factor"])

    if case == :splitgrid_alloc
        params = _params(cfg, dim)
        return let params=params
            () -> begin
                vofgrid(params)
                nothing
            end
        end
    elseif case == :initfgrid
        params = _params(cfg, dim)
        vg = vofgrid(params)
        shape = _shape(dim)
        return let vg=vg, shape=shape, nc=nc
            () -> begin
                initfgrid!(vg, shape; nc=nc)
                nothing
            end
        end
    elseif case == :youngs_normals
        vg, _ = _prepare_state(cfg, dim)
        return let vg=vg
            () -> begin
                youngs_normals!(vg)
                nothing
            end
        end
    elseif case == :reconstruct
        vg, params = _prepare_state(cfg, dim)
        return let vg=vg, params=params
            () -> begin
                reconstruct!(vg, params)
                nothing
            end
        end
    elseif case == :compute_dt
        vg, params = _prepare_reconstructed_state(cfg, dim)
        vel = _velocity(dim)
        return let vg=vg, vel=vel, params=params
            () -> compdt(vg, vel, params)
        end
    elseif case == :step
        vg, params = _prepare_reconstructed_state(cfg, dim)
        vel = _velocity(dim)
        fract0 = copy(fractg(vg))
        dt = dt_factor * compdt(vg, vel, params)
        return let vg=vg, params=params, vel=vel, fract0=fract0, dt=dt
            () -> begin
                fractg(vg) .= fract0
                vg.t = 0.0
                vg.step = 0
                vofadv!(vg, vel, dt, params)
                nothing
            end
        end
    elseif case == :integrate
        vg, params = _prepare_reconstructed_state(cfg, dim)
        vel = _velocity(dim)
        fract0 = copy(fractg(vg))
        dt = dt_factor * compdt(vg, vel, params)
        tf = 2.0 * dt
        return let vg=vg, params=params, vel=vel, fract0=fract0, tf=tf
            () -> begin
                fractg(vg) .= fract0
                vg.t = 0.0
                vg.step = 0
                integrate!(vg, vel, tf, params)
                nothing
            end
        end
    elseif case == :pipeline
        params = _params(cfg, dim)
        shape = _shape(dim)
        vel = _velocity(dim)
        return let params=params, shape=shape, vel=vel, nc=nc, dt_factor=dt_factor
            () -> begin
                vg = vofgrid(params)
                initfgrid!(vg, shape; nc=nc)
                reconstruct!(vg, params)
                dt = dt_factor * compdt(vg, vel, params)
                vofadv!(vg, vel, dt, params)
                nothing
            end
        end
    else
        error("Unknown case: $(case). Valid cases: $(collect(CASES)) or :all")
    end
end

function benchmark_profile(f::F; case::Symbol, cfg::Dict{String,Any}) where {F}
    runs = Int(cfg["runs"])
    profile_runs = Int(cfg["profile_runs"])

    # Warm-up compilation.
    f()

    elapsed = @elapsed begin
        for _ in 1:runs
            f()
        end
    end

    total_alloc_bytes = @allocated begin
        for _ in 1:runs
            f()
        end
    end

    alloc_per_run = total_alloc_bytes / runs
    alloc_free = total_alloc_bytes == 0

    println("Case: $(case) | dim=$(cfg["dim"]) | threads=$(nthreads())")
    @printf("  runs=%d  time_total=%.6f s  time_avg=%.3f us/run\n",
            runs, elapsed, elapsed * 1e6 / runs)
    @printf("  alloc_total=%d B  alloc_avg=%.2f B/run  alloc_free=%s\n",
            total_alloc_bytes, alloc_per_run, alloc_free)

    if profile_runs > 0
        Profile.clear()
        @profile begin
            for _ in 1:profile_runs
                f()
            end
        end
        println("  profile (flat, sorted by count):")
        Profile.print(format=:flat, sortedby=:count, mincount=5)
    end

    return nothing
end

function main(args=ARGS)
    cfg = _parse_args(args)

    if cfg["case"] == :all
        for case in CASES
            f = build_case(case; cfg=cfg)
            benchmark_profile(f; case=case, cfg=cfg)
            println()
        end
    else
        f = build_case(cfg["case"]; cfg=cfg)
        benchmark_profile(f; case=cfg["case"], cfg=cfg)
    end

    println("Tip: allocation-free target means alloc_avg ~ 0 B/run for hot kernels.")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
