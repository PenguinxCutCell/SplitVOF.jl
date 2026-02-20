function compute_dt(vg::SplitVOFGrid, velocity, t::Real=vg.t; cfl::Real=0.5)
    g = vg.grid
    tloc = Float64(t)

    umax = 0.0
    vmax = 0.0

    for j in 1:g.ny
        y = g.ylo + (j - 0.5) * g.dy
        for iface in 1:(g.nx + 1)
            x = g.xlo + (iface - 1) * g.dx
            vel = _velocity_at(velocity, SVector{2,Float64}(x, y), tloc)
            umax = max(umax, abs(vel[1]))
        end
    end

    for jface in 1:(g.ny + 1)
        y = g.ylo + (jface - 1) * g.dy
        for i in 1:g.nx
            x = g.xlo + (i - 0.5) * g.dx
            vel = _velocity_at(velocity, SVector{2,Float64}(x, y), tloc)
            vmax = max(vmax, abs(vel[2]))
        end
    end

    dtx = umax <= eps(Float64) ? Inf : g.dx / umax
    dty = vmax <= eps(Float64) ? Inf : g.dy / vmax

    return Float64(cfl) * min(dtx, dty)
end

function compute_dt(vg::SplitVOFGrid3D, velocity, t::Real=vg.t; cfl::Real=0.5)
    g = vg.grid
    tloc = Float64(t)

    umax = 0.0
    vmax = 0.0
    wmax = 0.0

    for k in 1:g.nz, j in 1:g.ny
        y = g.ylo + (j - 0.5) * g.dy
        z = g.zlo + (k - 0.5) * g.dz
        for iface in 1:(g.nx + 1)
            x = g.xlo + (iface - 1) * g.dx
            vel = _velocity_at(velocity, SVector{3,Float64}(x, y, z), tloc)
            umax = max(umax, abs(vel[1]))
        end
    end

    for k in 1:g.nz, jface in 1:(g.ny + 1)
        y = g.ylo + (jface - 1) * g.dy
        z = g.zlo + (k - 0.5) * g.dz
        for i in 1:g.nx
            x = g.xlo + (i - 0.5) * g.dx
            vel = _velocity_at(velocity, SVector{3,Float64}(x, y, z), tloc)
            vmax = max(vmax, abs(vel[2]))
        end
    end

    for kface in 1:(g.nz + 1), j in 1:g.ny, i in 1:g.nx
        z = g.zlo + (kface - 1) * g.dz
        x = g.xlo + (i - 0.5) * g.dx
        y = g.ylo + (j - 0.5) * g.dy
        vel = _velocity_at(velocity, SVector{3,Float64}(x, y, z), tloc)
        wmax = max(wmax, abs(vel[3]))
    end

    dtx = umax <= eps(Float64) ? Inf : g.dx / umax
    dty = vmax <= eps(Float64) ? Inf : g.dy / vmax
    dtz = wmax <= eps(Float64) ? Inf : g.dz / wmax

    return Float64(cfl) * min(dtx, dty, dtz)
end

function step!(vg::SplitVOFGrid, velocity, dt::Real)
    dtf = Float64(dt)
    dtf >= 0.0 || throw(ArgumentError("dt must be non-negative"))
    dtf == 0.0 && return vg

    t0 = vg.t

    if iseven(vg.step)
        reconstruct!(vg)
        _x_sweep!(vg, velocity, 0.5 * dtf, t0 + 0.25 * dtf)

        reconstruct!(vg)
        _y_sweep!(vg, velocity, dtf, t0 + 0.50 * dtf)

        reconstruct!(vg)
        _x_sweep!(vg, velocity, 0.5 * dtf, t0 + 0.75 * dtf)
    else
        reconstruct!(vg)
        _y_sweep!(vg, velocity, 0.5 * dtf, t0 + 0.25 * dtf)

        reconstruct!(vg)
        _x_sweep!(vg, velocity, dtf, t0 + 0.50 * dtf)

        reconstruct!(vg)
        _y_sweep!(vg, velocity, 0.5 * dtf, t0 + 0.75 * dtf)
    end

    vg.t = t0 + dtf
    vg.step += 1
    return vg
end

function _split_order3d(step::Int)
    r = step % 3
    if r == 0
        return (1, 2, 3)
    elseif r == 1
        return (2, 3, 1)
    else
        return (3, 1, 2)
    end
end

function _apply_sweep!(vg::SplitVOFGrid3D, axis::Int, velocity, dt::Float64, t::Float64)
    if axis == 1
        return _x_sweep!(vg, velocity, dt, t)
    elseif axis == 2
        return _y_sweep!(vg, velocity, dt, t)
    else
        return _z_sweep!(vg, velocity, dt, t)
    end
end

function step!(vg::SplitVOFGrid3D, velocity, dt::Real)
    dtf = Float64(dt)
    dtf >= 0.0 || throw(ArgumentError("dt must be non-negative"))
    dtf == 0.0 && return vg

    tmid = vg.t + 0.5 * dtf
    a, b, c = _split_order3d(vg.step)

    reconstruct!(vg)
    _apply_sweep!(vg, a, velocity, 0.5 * dtf, tmid)

    reconstruct!(vg)
    _apply_sweep!(vg, b, velocity, 0.5 * dtf, tmid)

    reconstruct!(vg)
    _apply_sweep!(vg, c, velocity, dtf, tmid)

    reconstruct!(vg)
    _apply_sweep!(vg, b, velocity, 0.5 * dtf, tmid)

    reconstruct!(vg)
    _apply_sweep!(vg, a, velocity, 0.5 * dtf, tmid)

    vg.t += dtf
    vg.step += 1
    return vg
end

function integrate!(vg::AbstractSplitVOFGrid, velocity, tf::Real;
                    dtmax::Real=Inf, cfl::Real=0.5)
    tf64 = Float64(tf)
    tf64 >= vg.t || throw(ArgumentError("tf must be >= current time"))

    while vg.t <= tf64 - eps(tf64)
        dt_cfl = compute_dt(vg, velocity, vg.t; cfl=cfl)
        dt = min(Float64(dtmax), dt_cfl, tf64 - vg.t)

        if !isfinite(dt)
            dt = tf64 - vg.t
        end

        dt > 0.0 || break
        step!(vg, velocity, dt)
    end

    return vg
end

compdt(vg::AbstractSplitVOFGrid, velocity, t::Real=vg.t; cfl::Real=0.5) =
    compute_dt(vg, velocity, t; cfl=cfl)

compdt(vg::AbstractSplitVOFGrid, velocity, params::SplitVOFParams, t::Real=vg.t) =
    compute_dt(vg, velocity, t; cfl=params.cfl)

vofadv!(vg::AbstractSplitVOFGrid, velocity, dt::Real) = step!(vg, velocity, dt)

function vofadv!(vg::AbstractSplitVOFGrid, velocity, dt::Real, params::SplitVOFParams)
    params.advection isa StrangSplit ||
        throw(ArgumentError("Only StrangSplit() is currently implemented"))
    return step!(vg, velocity, dt)
end

function integrate!(vg::AbstractSplitVOFGrid, velocity, tf::Real, params::SplitVOFParams;
                    dtmax::Real=Inf)
    params.advection isa StrangSplit ||
        throw(ArgumentError("Only StrangSplit() is currently implemented"))
    return integrate!(vg, velocity, tf; dtmax=dtmax, cfl=params.cfl)
end
