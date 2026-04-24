"""
Compute a stable time step from face-velocity CFL constraints.

Supports velocity call signatures `(x,t)`, `(x,y,t)/(x,y,z,t)`,
`(x)`, or coordinate-only variants.
"""
function compute_dt(vg::SplitVOFGrid, velocity::F, t::Real=vg.t; cfl::Real=0.5) where {F<:Function}
    g = vg.grid
    tloc = Float64(t)
    mode = _velocity_mode_2d(velocity, tloc)

    umax = 0.0
    vmax = 0.0

    if mode == 1
        for j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                vel = velocity(SVector{2,Float64}(x, y), tloc)
                umax = max(umax, abs(Float64(vel[1])))
            end
        end

        for jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vel = velocity(SVector{2,Float64}(x, y), tloc)
                vmax = max(vmax, abs(Float64(vel[2])))
            end
        end
    elseif mode == 2
        for j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                vel = velocity(x, y, tloc)
                umax = max(umax, abs(Float64(vel[1])))
            end
        end

        for jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vel = velocity(x, y, tloc)
                vmax = max(vmax, abs(Float64(vel[2])))
            end
        end
    elseif mode == 3
        for j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                vel = velocity(SVector{2,Float64}(x, y))
                umax = max(umax, abs(Float64(vel[1])))
            end
        end

        for jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vel = velocity(SVector{2,Float64}(x, y))
                vmax = max(vmax, abs(Float64(vel[2])))
            end
        end
    else
        for j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                vel = velocity(x, y)
                umax = max(umax, abs(Float64(vel[1])))
            end
        end

        for jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vel = velocity(x, y)
                vmax = max(vmax, abs(Float64(vel[2])))
            end
        end
    end

    dtx = umax <= eps(Float64) ? Inf : g.dx / umax
    dty = vmax <= eps(Float64) ? Inf : g.dy / vmax

    return Float64(cfl) * min(dtx, dty)
end

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

function compute_dt(vg::SplitVOFGrid3D, velocity::F, t::Real=vg.t; cfl::Real=0.5) where {F<:Function}
    g = vg.grid
    tloc = Float64(t)
    mode = _velocity_mode_3d(velocity, tloc)

    umax = 0.0
    vmax = 0.0
    wmax = 0.0

    if mode == 1
        for k in 1:g.nz, j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                vel = velocity(SVector{3,Float64}(x, y, z), tloc)
                umax = max(umax, abs(Float64(vel[1])))
            end
        end

        for k in 1:g.nz, jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vel = velocity(SVector{3,Float64}(x, y, z), tloc)
                vmax = max(vmax, abs(Float64(vel[2])))
            end
        end

        for kface in 1:(g.nz + 1), j in 1:g.ny, i in 1:g.nx
            z = g.zlo + (kface - 1) * g.dz
            x = g.xlo + (i - 0.5) * g.dx
            y = g.ylo + (j - 0.5) * g.dy
            vel = velocity(SVector{3,Float64}(x, y, z), tloc)
            wmax = max(wmax, abs(Float64(vel[3])))
        end
    elseif mode == 2
        for k in 1:g.nz, j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                vel = velocity(x, y, z, tloc)
                umax = max(umax, abs(Float64(vel[1])))
            end
        end

        for k in 1:g.nz, jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vel = velocity(x, y, z, tloc)
                vmax = max(vmax, abs(Float64(vel[2])))
            end
        end

        for kface in 1:(g.nz + 1), j in 1:g.ny, i in 1:g.nx
            z = g.zlo + (kface - 1) * g.dz
            x = g.xlo + (i - 0.5) * g.dx
            y = g.ylo + (j - 0.5) * g.dy
            vel = velocity(x, y, z, tloc)
            wmax = max(wmax, abs(Float64(vel[3])))
        end
    elseif mode == 3
        for k in 1:g.nz, j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                vel = velocity(SVector{3,Float64}(x, y, z))
                umax = max(umax, abs(Float64(vel[1])))
            end
        end

        for k in 1:g.nz, jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vel = velocity(SVector{3,Float64}(x, y, z))
                vmax = max(vmax, abs(Float64(vel[2])))
            end
        end

        for kface in 1:(g.nz + 1), j in 1:g.ny, i in 1:g.nx
            z = g.zlo + (kface - 1) * g.dz
            x = g.xlo + (i - 0.5) * g.dx
            y = g.ylo + (j - 0.5) * g.dy
            vel = velocity(SVector{3,Float64}(x, y, z))
            wmax = max(wmax, abs(Float64(vel[3])))
        end
    else
        for k in 1:g.nz, j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                vel = velocity(x, y, z)
                umax = max(umax, abs(Float64(vel[1])))
            end
        end

        for k in 1:g.nz, jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vel = velocity(x, y, z)
                vmax = max(vmax, abs(Float64(vel[2])))
            end
        end

        for kface in 1:(g.nz + 1), j in 1:g.ny, i in 1:g.nx
            z = g.zlo + (kface - 1) * g.dz
            x = g.xlo + (i - 0.5) * g.dx
            y = g.ylo + (j - 0.5) * g.dy
            vel = velocity(x, y, z)
            wmax = max(wmax, abs(Float64(vel[3])))
        end
    end

    dtx = umax <= eps(Float64) ? Inf : g.dx / umax
    dty = vmax <= eps(Float64) ? Inf : g.dy / vmax
    dtz = wmax <= eps(Float64) ? Inf : g.dz / wmax

    return Float64(cfl) * min(dtx, dty, dtz)
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

"""
Advance one split-advection step in 2D using Strang splitting.
"""
function step!(vg::SplitVOFGrid, velocity, dt::Real)
    return step!(vg, velocity, dt, YoungsPLIC())
end

function step!(vg::SplitVOFGrid, velocity, dt::Real, reconstruction::ReconstructionMethod)
    dtf = Float64(dt)
    dtf >= 0.0 || throw(ArgumentError("dt must be non-negative"))
    dtf == 0.0 && return vg

    t0 = vg.t

    if iseven(vg.step)
        reconstruct!(vg, reconstruction)
        _x_sweep!(vg, velocity, 0.5 * dtf, t0 + 0.25 * dtf)

        reconstruct!(vg, reconstruction)
        _y_sweep!(vg, velocity, dtf, t0 + 0.50 * dtf)

        reconstruct!(vg, reconstruction)
        _x_sweep!(vg, velocity, 0.5 * dtf, t0 + 0.75 * dtf)
    else
        reconstruct!(vg, reconstruction)
        _y_sweep!(vg, velocity, 0.5 * dtf, t0 + 0.25 * dtf)

        reconstruct!(vg, reconstruction)
        _x_sweep!(vg, velocity, dtf, t0 + 0.50 * dtf)

        reconstruct!(vg, reconstruction)
        _y_sweep!(vg, velocity, 0.5 * dtf, t0 + 0.75 * dtf)
    end

    vg.t = t0 + dtf
    vg.step += 1
    return vg
end

"""
Advance one split-advection step in 2D using first-order Lie splitting.
"""
function _step_lie!(vg::SplitVOFGrid, velocity, dt::Real)
    return _step_lie!(vg, velocity, dt, YoungsPLIC())
end

function _step_lie!(vg::SplitVOFGrid, velocity, dt::Real, reconstruction::ReconstructionMethod)
    dtf = Float64(dt)
    dtf >= 0.0 || throw(ArgumentError("dt must be non-negative"))
    dtf == 0.0 && return vg

    tmid = vg.t + 0.5 * dtf

    if iseven(vg.step)
        reconstruct!(vg, reconstruction)
        _x_sweep!(vg, velocity, dtf, tmid)

        reconstruct!(vg, reconstruction)
        _y_sweep!(vg, velocity, dtf, tmid)
    else
        reconstruct!(vg, reconstruction)
        _y_sweep!(vg, velocity, dtf, tmid)

        reconstruct!(vg, reconstruction)
        _x_sweep!(vg, velocity, dtf, tmid)
    end

    vg.t += dtf
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

function _lie_order3d(step::Int)
    r = step % 3
    if r == 0
        return (1, 2, 3)
    elseif r == 1
        return (2, 3, 1)
    else
        return (3, 1, 2)
    end
end

"""
Advance one split-advection step in 3D using symmetric split sweeps.
"""
function step!(vg::SplitVOFGrid3D, velocity, dt::Real)
    return step!(vg, velocity, dt, YoungsPLIC())
end

function step!(vg::SplitVOFGrid3D, velocity, dt::Real, reconstruction::ReconstructionMethod)
    dtf = Float64(dt)
    dtf >= 0.0 || throw(ArgumentError("dt must be non-negative"))
    dtf == 0.0 && return vg

    tmid = vg.t + 0.5 * dtf
    a, b, c = _split_order3d(vg.step)

    reconstruct!(vg, reconstruction)
    _apply_sweep!(vg, a, velocity, 0.5 * dtf, tmid)

    reconstruct!(vg, reconstruction)
    _apply_sweep!(vg, b, velocity, 0.5 * dtf, tmid)

    reconstruct!(vg, reconstruction)
    _apply_sweep!(vg, c, velocity, dtf, tmid)

    reconstruct!(vg, reconstruction)
    _apply_sweep!(vg, b, velocity, 0.5 * dtf, tmid)

    reconstruct!(vg, reconstruction)
    _apply_sweep!(vg, a, velocity, 0.5 * dtf, tmid)

    vg.t += dtf
    vg.step += 1
    return vg
end

"""
Advance one split-advection step in 3D using first-order Lie splitting.
"""
function _step_lie!(vg::SplitVOFGrid3D, velocity, dt::Real)
    return _step_lie!(vg, velocity, dt, YoungsPLIC())
end

function _step_lie!(vg::SplitVOFGrid3D, velocity, dt::Real, reconstruction::ReconstructionMethod)
    dtf = Float64(dt)
    dtf >= 0.0 || throw(ArgumentError("dt must be non-negative"))
    dtf == 0.0 && return vg

    tmid = vg.t + 0.5 * dtf
    a, b, c = _lie_order3d(vg.step)

    reconstruct!(vg, reconstruction)
    _apply_sweep!(vg, a, velocity, dtf, tmid)

    reconstruct!(vg, reconstruction)
    _apply_sweep!(vg, b, velocity, dtf, tmid)

    reconstruct!(vg, reconstruction)
    _apply_sweep!(vg, c, velocity, dtf, tmid)

    vg.t += dtf
    vg.step += 1
    return vg
end

"""
Integrate from current time `vg.t` to final time `tf`.
"""
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

"""
Alias for `compute_dt`.
"""
compdt(vg::AbstractSplitVOFGrid, velocity, t::Real=vg.t; cfl::Real=0.5) =
    compute_dt(vg, velocity, t; cfl=cfl)

compdt(vg::AbstractSplitVOFGrid, velocity, params::SplitVOFParams, t::Real=vg.t) =
    compute_dt(vg, velocity, t; cfl=params.cfl)

"""
Alias for `step!`.
"""
vofadv!(vg::AbstractSplitVOFGrid, velocity, dt::Real) = step!(vg, velocity, dt)

"""
Advance one advection step using algorithm choice in `params`.
"""
function vofadv!(vg::AbstractSplitVOFGrid, velocity, dt::Real, params::SplitVOFParams)
    if params.advection isa StrangSplit
        return step!(vg, velocity, dt, params.reconstruction)
    elseif params.advection isa LieSplit
        return _step_lie!(vg, velocity, dt, params.reconstruction)
    else
        throw(ArgumentError("Unsupported advection method $(typeof(params.advection))"))
    end
end

"""
Integrate to `tf` using algorithm choice in `params`.
"""
function integrate!(vg::AbstractSplitVOFGrid, velocity, tf::Real, params::SplitVOFParams;
                    dtmax::Real=Inf)
    tf64 = Float64(tf)
    tf64 >= vg.t || throw(ArgumentError("tf must be >= current time"))

    if params.advection isa StrangSplit
        while vg.t <= tf64 - eps(tf64)
            dt_cfl = compute_dt(vg, velocity, vg.t; cfl=params.cfl)
            dt = min(Float64(dtmax), dt_cfl, tf64 - vg.t)

            if !isfinite(dt)
                dt = tf64 - vg.t
            end

            dt > 0.0 || break
            step!(vg, velocity, dt, params.reconstruction)
        end
    elseif params.advection isa LieSplit
        while vg.t <= tf64 - eps(tf64)
            dt_cfl = compute_dt(vg, velocity, vg.t; cfl=params.cfl)
            dt = min(Float64(dtmax), dt_cfl, tf64 - vg.t)

            if !isfinite(dt)
                dt = tf64 - vg.t
            end

            dt > 0.0 || break
            _step_lie!(vg, velocity, dt, params.reconstruction)
        end
    else
        throw(ArgumentError("Unsupported advection method $(typeof(params.advection))"))
    end

    return vg
end
