@inline _to_svec2(v::SVector{2,<:Real}) = SVector{2,Float64}(v[1], v[2])
@inline _to_svec2(v::NTuple{2,<:Real}) = SVector{2,Float64}(v[1], v[2])
@inline _to_svec3(v::SVector{3,<:Real}) = SVector{3,Float64}(v[1], v[2], v[3])
@inline _to_svec3(v::NTuple{3,<:Real}) = SVector{3,Float64}(v[1], v[2], v[3])

@inline function _velocity_at(v::SVector{2,<:Real}, x::SVector{2,Float64}, t::Float64)
    return _to_svec2(v)
end

@inline function _velocity_at(v::NTuple{2,<:Real}, x::SVector{2,Float64}, t::Float64)
    return _to_svec2(v)
end

@inline function _velocity_at(v::SVector{3,<:Real}, x::SVector{3,Float64}, t::Float64)
    return _to_svec3(v)
end

@inline function _velocity_at(v::NTuple{3,<:Real}, x::SVector{3,Float64}, t::Float64)
    return _to_svec3(v)
end

function _velocity_at(v::Function, x::SVector{2,Float64}, t::Float64)
    if applicable(v, x, t)
        return _to_svec2(v(x, t))
    elseif applicable(v, x[1], x[2], t)
        return _to_svec2(v(x[1], x[2], t))
    elseif applicable(v, x)
        return _to_svec2(v(x))
    elseif applicable(v, x[1], x[2])
        return _to_svec2(v(x[1], x[2]))
    else
        throw(ArgumentError("velocity must accept (x,t), (x,y,t), (x), or (x,y) for 2D"))
    end
end

function _velocity_at(v::Function, x::SVector{3,Float64}, t::Float64)
    if applicable(v, x, t)
        return _to_svec3(v(x, t))
    elseif applicable(v, x[1], x[2], x[3], t)
        return _to_svec3(v(x[1], x[2], x[3], t))
    elseif applicable(v, x)
        return _to_svec3(v(x))
    elseif applicable(v, x[1], x[2], x[3])
        return _to_svec3(v(x[1], x[2], x[3]))
    else
        throw(ArgumentError("velocity must accept (x,t), (x,y,z,t), (x), or (x,y,z) for 3D"))
    end
end

@inline function _velocity_mode_2d(v::F, t::Float64) where {F<:Function}
    tls = task_local_storage()
    cache_any = get(tls, :_splitvof_velocity_mode_2d, nothing)
    cache = if cache_any === nothing
        c = Dict{DataType,Int}()
        tls[:_splitvof_velocity_mode_2d] = c
        c
    else
        cache_any::Dict{DataType,Int}
    end

    T = typeof(v)
    mode = get(cache, T, 0)
    if mode != 0
        return mode
    end

    x = SVector{2,Float64}(0.0, 0.0)
    if applicable(v, x, t)
        mode = 1
    elseif applicable(v, x[1], x[2], t)
        mode = 2
    elseif applicable(v, x)
        mode = 3
    elseif applicable(v, x[1], x[2])
        mode = 4
    else
        throw(ArgumentError("velocity must accept (x,t), (x,y,t), (x), or (x,y) for 2D"))
    end

    cache[T] = mode
    return mode
end

@inline function _velocity_mode_3d(v::F, t::Float64) where {F<:Function}
    tls = task_local_storage()
    cache_any = get(tls, :_splitvof_velocity_mode_3d, nothing)
    cache = if cache_any === nothing
        c = Dict{DataType,Int}()
        tls[:_splitvof_velocity_mode_3d] = c
        c
    else
        cache_any::Dict{DataType,Int}
    end

    T = typeof(v)
    mode = get(cache, T, 0)
    if mode != 0
        return mode
    end

    x = SVector{3,Float64}(0.0, 0.0, 0.0)
    if applicable(v, x, t)
        mode = 1
    elseif applicable(v, x[1], x[2], x[3], t)
        mode = 2
    elseif applicable(v, x)
        mode = 3
    elseif applicable(v, x[1], x[2], x[3])
        mode = 4
    else
        throw(ArgumentError("velocity must accept (x,t), (x,y,z,t), (x), or (x,y,z) for 3D"))
    end

    cache[T] = mode
    return mode
end

mutable struct _Sweep2DWork
    uface::Matrix{Float64}
    vface::Matrix{Float64}
    fluxx::Matrix{Float64}
    fluxy::Matrix{Float64}
    Anew::Matrix{Float64}
    poly::Polygon2D
end

mutable struct _Sweep3DWork
    uface::Array{Float64,3}
    vface::Array{Float64,3}
    wface::Array{Float64,3}
    fluxx::Array{Float64,3}
    fluxy::Array{Float64,3}
    fluxz::Array{Float64,3}
    Anew::Array{Float64,3}
    poly::Polyhedron3D
end

function _Sweep2DWork(vg::SplitVOFGrid)
    g = vg.grid
    return _Sweep2DWork(
        zeros(Float64, g.nx + 1, g.ny),
        zeros(Float64, g.nx, g.ny + 1),
        zeros(Float64, g.nx + 1, g.ny),
        zeros(Float64, g.nx, g.ny + 1),
        similar(vg.fract),
        _rect_polygon(0.0, 1.0, 0.0, 1.0),
    )
end

function _Sweep3DWork(vg::SplitVOFGrid3D)
    g = vg.grid
    return _Sweep3DWork(
        zeros(Float64, g.nx + 1, g.ny, g.nz),
        zeros(Float64, g.nx, g.ny + 1, g.nz),
        zeros(Float64, g.nx, g.ny, g.nz + 1),
        zeros(Float64, g.nx + 1, g.ny, g.nz),
        zeros(Float64, g.nx, g.ny + 1, g.nz),
        zeros(Float64, g.nx, g.ny, g.nz + 1),
        similar(vg.fract),
        _rect_polyhedron(0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
    )
end

function _get_sweep2d_work(vg::SplitVOFGrid)
    tls = task_local_storage()
    cache_any = get(tls, :_splitvof_sweep2d_work, nothing)
    cache = if cache_any === nothing
        c = Dict{UInt, _Sweep2DWork}()
        tls[:_splitvof_sweep2d_work] = c
        c
    else
        cache_any::Dict{UInt, _Sweep2DWork}
    end

    key = objectid(vg)
    work = get(cache, key, nothing)
    if work === nothing || size(work.Anew) != size(vg.fract)
        work = _Sweep2DWork(vg)
        cache[key] = work
    end
    return work
end

function _get_sweep3d_work(vg::SplitVOFGrid3D)
    tls = task_local_storage()
    cache_any = get(tls, :_splitvof_sweep3d_work, nothing)
    cache = if cache_any === nothing
        c = Dict{UInt, _Sweep3DWork}()
        tls[:_splitvof_sweep3d_work] = c
        c
    else
        cache_any::Dict{UInt, _Sweep3DWork}
    end

    key = objectid(vg)
    work = get(cache, key, nothing)
    if work === nothing || size(work.Anew) != size(vg.fract)
        work = _Sweep3DWork(vg)
        cache[key] = work
    end
    return work
end

function _x_face_velocities!(uface::AbstractArray{Float64,2},
                             vg::SplitVOFGrid, velocity::F, t::Float64) where {F<:Function}
    g = vg.grid
    mode = _velocity_mode_2d(velocity, t)
    if mode == 1
        @inbounds for j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                uface[iface, j] = Float64(velocity(SVector{2,Float64}(x, y), t)[1])
            end
        end
    elseif mode == 2
        @inbounds for j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                uface[iface, j] = Float64(velocity(x, y, t)[1])
            end
        end
    elseif mode == 3
        @inbounds for j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                uface[iface, j] = Float64(velocity(SVector{2,Float64}(x, y))[1])
            end
        end
    else
        @inbounds for j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                uface[iface, j] = Float64(velocity(x, y)[1])
            end
        end
    end
    return uface
end

function _x_face_velocities!(uface::AbstractArray{Float64,2},
                             vg::SplitVOFGrid, velocity, t::Float64)
    g = vg.grid
    @inbounds for j in 1:g.ny
        y = g.ylo + (j - 0.5) * g.dy
        for iface in 1:(g.nx + 1)
            x = g.xlo + (iface - 1) * g.dx
            vel = _velocity_at(velocity, SVector{2,Float64}(x, y), t)
            uface[iface, j] = vel[1]
        end
    end
    return uface
end

function _y_face_velocities!(vface::AbstractArray{Float64,2},
                             vg::SplitVOFGrid, velocity::F, t::Float64) where {F<:Function}
    g = vg.grid
    mode = _velocity_mode_2d(velocity, t)
    if mode == 1
        @inbounds for jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vface[i, jface] = Float64(velocity(SVector{2,Float64}(x, y), t)[2])
            end
        end
    elseif mode == 2
        @inbounds for jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vface[i, jface] = Float64(velocity(x, y, t)[2])
            end
        end
    elseif mode == 3
        @inbounds for jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vface[i, jface] = Float64(velocity(SVector{2,Float64}(x, y))[2])
            end
        end
    else
        @inbounds for jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vface[i, jface] = Float64(velocity(x, y)[2])
            end
        end
    end
    return vface
end

function _y_face_velocities!(vface::AbstractArray{Float64,2},
                             vg::SplitVOFGrid, velocity, t::Float64)
    g = vg.grid
    @inbounds for jface in 1:(g.ny + 1)
        y = g.ylo + (jface - 1) * g.dy
        for i in 1:g.nx
            x = g.xlo + (i - 0.5) * g.dx
            vel = _velocity_at(velocity, SVector{2,Float64}(x, y), t)
            vface[i, jface] = vel[2]
        end
    end
    return vface
end

function _x_face_velocities!(uface::AbstractArray{Float64,3},
                             vg::SplitVOFGrid3D, velocity::F, t::Float64) where {F<:Function}
    g = vg.grid
    mode = _velocity_mode_3d(velocity, t)
    if mode == 1
        @inbounds for k in 1:g.nz, j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                uface[iface, j, k] = Float64(velocity(SVector{3,Float64}(x, y, z), t)[1])
            end
        end
    elseif mode == 2
        @inbounds for k in 1:g.nz, j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                uface[iface, j, k] = Float64(velocity(x, y, z, t)[1])
            end
        end
    elseif mode == 3
        @inbounds for k in 1:g.nz, j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                uface[iface, j, k] = Float64(velocity(SVector{3,Float64}(x, y, z))[1])
            end
        end
    else
        @inbounds for k in 1:g.nz, j in 1:g.ny
            y = g.ylo + (j - 0.5) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for iface in 1:(g.nx + 1)
                x = g.xlo + (iface - 1) * g.dx
                uface[iface, j, k] = Float64(velocity(x, y, z)[1])
            end
        end
    end
    return uface
end

function _x_face_velocities!(uface::AbstractArray{Float64,3},
                             vg::SplitVOFGrid3D, velocity, t::Float64)
    g = vg.grid
    @inbounds for k in 1:g.nz, j in 1:g.ny
        y = g.ylo + (j - 0.5) * g.dy
        z = g.zlo + (k - 0.5) * g.dz
        for iface in 1:(g.nx + 1)
            x = g.xlo + (iface - 1) * g.dx
            vel = _velocity_at(velocity, SVector{3,Float64}(x, y, z), t)
            uface[iface, j, k] = vel[1]
        end
    end
    return uface
end

function _y_face_velocities!(vface::AbstractArray{Float64,3},
                             vg::SplitVOFGrid3D, velocity::F, t::Float64) where {F<:Function}
    g = vg.grid
    mode = _velocity_mode_3d(velocity, t)
    if mode == 1
        @inbounds for k in 1:g.nz, jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vface[i, jface, k] = Float64(velocity(SVector{3,Float64}(x, y, z), t)[2])
            end
        end
    elseif mode == 2
        @inbounds for k in 1:g.nz, jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vface[i, jface, k] = Float64(velocity(x, y, z, t)[2])
            end
        end
    elseif mode == 3
        @inbounds for k in 1:g.nz, jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vface[i, jface, k] = Float64(velocity(SVector{3,Float64}(x, y, z))[2])
            end
        end
    else
        @inbounds for k in 1:g.nz, jface in 1:(g.ny + 1)
            y = g.ylo + (jface - 1) * g.dy
            z = g.zlo + (k - 0.5) * g.dz
            for i in 1:g.nx
                x = g.xlo + (i - 0.5) * g.dx
                vface[i, jface, k] = Float64(velocity(x, y, z)[2])
            end
        end
    end
    return vface
end

function _y_face_velocities!(vface::AbstractArray{Float64,3},
                             vg::SplitVOFGrid3D, velocity, t::Float64)
    g = vg.grid
    @inbounds for k in 1:g.nz, jface in 1:(g.ny + 1)
        y = g.ylo + (jface - 1) * g.dy
        z = g.zlo + (k - 0.5) * g.dz
        for i in 1:g.nx
            x = g.xlo + (i - 0.5) * g.dx
            vel = _velocity_at(velocity, SVector{3,Float64}(x, y, z), t)
            vface[i, jface, k] = vel[2]
        end
    end
    return vface
end

function _z_face_velocities!(wface::AbstractArray{Float64,3},
                             vg::SplitVOFGrid3D, velocity::F, t::Float64) where {F<:Function}
    g = vg.grid
    mode = _velocity_mode_3d(velocity, t)
    if mode == 1
        @inbounds for kface in 1:(g.nz + 1), j in 1:g.ny, i in 1:g.nx
            z = g.zlo + (kface - 1) * g.dz
            x = g.xlo + (i - 0.5) * g.dx
            y = g.ylo + (j - 0.5) * g.dy
            wface[i, j, kface] = Float64(velocity(SVector{3,Float64}(x, y, z), t)[3])
        end
    elseif mode == 2
        @inbounds for kface in 1:(g.nz + 1), j in 1:g.ny, i in 1:g.nx
            z = g.zlo + (kface - 1) * g.dz
            x = g.xlo + (i - 0.5) * g.dx
            y = g.ylo + (j - 0.5) * g.dy
            wface[i, j, kface] = Float64(velocity(x, y, z, t)[3])
        end
    elseif mode == 3
        @inbounds for kface in 1:(g.nz + 1), j in 1:g.ny, i in 1:g.nx
            z = g.zlo + (kface - 1) * g.dz
            x = g.xlo + (i - 0.5) * g.dx
            y = g.ylo + (j - 0.5) * g.dy
            wface[i, j, kface] = Float64(velocity(SVector{3,Float64}(x, y, z))[3])
        end
    else
        @inbounds for kface in 1:(g.nz + 1), j in 1:g.ny, i in 1:g.nx
            z = g.zlo + (kface - 1) * g.dz
            x = g.xlo + (i - 0.5) * g.dx
            y = g.ylo + (j - 0.5) * g.dy
            wface[i, j, kface] = Float64(velocity(x, y, z)[3])
        end
    end
    return wface
end

function _z_face_velocities!(wface::AbstractArray{Float64,3},
                             vg::SplitVOFGrid3D, velocity, t::Float64)
    g = vg.grid
    @inbounds for kface in 1:(g.nz + 1), j in 1:g.ny, i in 1:g.nx
        z = g.zlo + (kface - 1) * g.dz
        x = g.xlo + (i - 0.5) * g.dx
        y = g.ylo + (j - 0.5) * g.dy
        vel = _velocity_at(velocity, SVector{3,Float64}(x, y, z), t)
        wface[i, j, kface] = vel[3]
    end
    return wface
end

function _material_area_x(vg::SplitVOFGrid, i::Int, j::Int,
                          delta::Float64, side::Int,
                          poly::Polygon2D)
    delta <= 0.0 && return 0.0

    g = vg.grid
    if is_empty(vg, i, j)
        return 0.0
    elseif is_full(vg, i, j)
        return delta * g.dy
    end

    xL, xR, yB, yT = cell_bounds(g, i, j)
    _rect_polygon!(poly, xL, xR, yB, yT)

    _, icontp = inte2d!(poly, vg.cplic[i, j], vg.nxint[i, j], vg.nyint[i, j])
    icontp == 0 && return 0.0

    if side > 0
        x0 = xR - delta
        _, icontp_slab = inte2d!(poly, -x0, 1.0, 0.0)
        icontp_slab == 0 && return 0.0
    else
        x1 = xL + delta
        _, icontp_slab = inte2d!(poly, x1, -1.0, 0.0)
        icontp_slab == 0 && return 0.0
    end

    area = toolv2d(poly)
    if !isfinite(area)
        return 0.0
    end
    return clamp(area, 0.0, delta * g.dy)
end

function _material_area_y(vg::SplitVOFGrid, i::Int, j::Int,
                          delta::Float64, side::Int,
                          poly::Polygon2D)
    delta <= 0.0 && return 0.0

    g = vg.grid
    if is_empty(vg, i, j)
        return 0.0
    elseif is_full(vg, i, j)
        return delta * g.dx
    end

    xL, xR, yB, yT = cell_bounds(g, i, j)
    _rect_polygon!(poly, xL, xR, yB, yT)

    _, icontp = inte2d!(poly, vg.cplic[i, j], vg.nxint[i, j], vg.nyint[i, j])
    icontp == 0 && return 0.0

    if side > 0
        y0 = yT - delta
        _, icontp_slab = inte2d!(poly, -y0, 0.0, 1.0)
        icontp_slab == 0 && return 0.0
    else
        y1 = yB + delta
        _, icontp_slab = inte2d!(poly, y1, 0.0, -1.0)
        icontp_slab == 0 && return 0.0
    end

    area = toolv2d(poly)
    if !isfinite(area)
        return 0.0
    end
    return clamp(area, 0.0, delta * g.dx)
end

function _material_volume_x(vg::SplitVOFGrid3D, i::Int, j::Int, k::Int,
                            delta::Float64, side::Int,
                            poly::Polyhedron3D)
    delta <= 0.0 && return 0.0

    g = vg.grid
    if is_empty(vg, i, j, k)
        return 0.0
    elseif is_full(vg, i, j, k)
        return delta * g.dy * g.dz
    end

    xL, xR, yB, yT, zD, zU = cell_bounds(g, i, j, k)
    _rect_polyhedron!(poly, xL, xR, yB, yT, zD, zU)

    _, icontp = inte3d!(poly,
                        vg.cplic[i, j, k],
                        vg.nxint[i, j, k],
                        vg.nyint[i, j, k],
                        vg.nzint[i, j, k])
    icontp == 0 && return 0.0

    if side > 0
        x0 = xR - delta
        _, icontp_slab = inte3d!(poly, -x0, 1.0, 0.0, 0.0)
        icontp_slab == 0 && return 0.0
    else
        x1 = xL + delta
        _, icontp_slab = inte3d!(poly, x1, -1.0, 0.0, 0.0)
        icontp_slab == 0 && return 0.0
    end

    vol = toolv3d(poly)
    if !isfinite(vol)
        return 0.0
    end
    return clamp(vol, 0.0, delta * g.dy * g.dz)
end

function _material_volume_y(vg::SplitVOFGrid3D, i::Int, j::Int, k::Int,
                            delta::Float64, side::Int,
                            poly::Polyhedron3D)
    delta <= 0.0 && return 0.0

    g = vg.grid
    if is_empty(vg, i, j, k)
        return 0.0
    elseif is_full(vg, i, j, k)
        return delta * g.dx * g.dz
    end

    xL, xR, yB, yT, zD, zU = cell_bounds(g, i, j, k)
    _rect_polyhedron!(poly, xL, xR, yB, yT, zD, zU)

    _, icontp = inte3d!(poly,
                        vg.cplic[i, j, k],
                        vg.nxint[i, j, k],
                        vg.nyint[i, j, k],
                        vg.nzint[i, j, k])
    icontp == 0 && return 0.0

    if side > 0
        y0 = yT - delta
        _, icontp_slab = inte3d!(poly, -y0, 0.0, 1.0, 0.0)
        icontp_slab == 0 && return 0.0
    else
        y1 = yB + delta
        _, icontp_slab = inte3d!(poly, y1, 0.0, -1.0, 0.0)
        icontp_slab == 0 && return 0.0
    end

    vol = toolv3d(poly)
    if !isfinite(vol)
        return 0.0
    end
    return clamp(vol, 0.0, delta * g.dx * g.dz)
end

function _material_volume_z(vg::SplitVOFGrid3D, i::Int, j::Int, k::Int,
                            delta::Float64, side::Int,
                            poly::Polyhedron3D)
    delta <= 0.0 && return 0.0

    g = vg.grid
    if is_empty(vg, i, j, k)
        return 0.0
    elseif is_full(vg, i, j, k)
        return delta * g.dx * g.dy
    end

    xL, xR, yB, yT, zD, zU = cell_bounds(g, i, j, k)
    _rect_polyhedron!(poly, xL, xR, yB, yT, zD, zU)

    _, icontp = inte3d!(poly,
                        vg.cplic[i, j, k],
                        vg.nxint[i, j, k],
                        vg.nyint[i, j, k],
                        vg.nzint[i, j, k])
    icontp == 0 && return 0.0

    if side > 0
        z0 = zU - delta
        _, icontp_slab = inte3d!(poly, -z0, 0.0, 0.0, 1.0)
        icontp_slab == 0 && return 0.0
    else
        z1 = zD + delta
        _, icontp_slab = inte3d!(poly, z1, 0.0, 0.0, -1.0)
        icontp_slab == 0 && return 0.0
    end

    vol = toolv3d(poly)
    if !isfinite(vol)
        return 0.0
    end
    return clamp(vol, 0.0, delta * g.dx * g.dy)
end

function _soft_clip!(A::AbstractArray{<:Real}, tol::Float64)
    @inbounds for idx in eachindex(A)
        a = Float64(A[idx])
        if a < -tol || a > 1.0 + tol
            throw(DomainError(a, "VOF fraction left [0,1] by more than tolerance"))
        elseif a < 0.0
            A[idx] = 0.0
        elseif a > 1.0
            A[idx] = 1.0
        end
    end
    return A
end

function _x_sweep!(vg::SplitVOFGrid, velocity, dt::Float64, t::Float64)
    dt <= 0.0 && return vg
    g = vg.grid
    work = _get_sweep2d_work(vg)

    uface = work.uface
    flux = work.fluxx
    Anew = work.Anew
    poly = work.poly

    _x_face_velocities!(uface, vg, velocity, t)
    fill!(flux, 0.0)

    @inbounds for j in 1:g.ny, iface in 1:(g.nx + 1)
        u = uface[iface, j]
        if abs(u) <= eps(Float64)
            continue
        end

        iL = iface == 1 ? g.nx : (iface - 1)
        iR = iface == (g.nx + 1) ? 1 : iface

        donor_i, side, sgn = if u > 0.0
            (iL, 1, 1.0)
        else
            (iR, -1, -1.0)
        end

        delta = min(abs(u) * dt, g.dx)
        areaf = _material_area_x(vg, donor_i, j, delta, side, poly)
        flux[iface, j] = sgn * areaf
    end

    Aold = vg.fract
    invvol = 1.0 / (g.dx * g.dy)
    @inbounds for j in 1:g.ny, i in 1:g.nx
        Anew[i, j] = Aold[i, j] - (flux[i + 1, j] - flux[i, j]) * invvol
    end

    _soft_clip!(Anew, vg.tol)
    copyto!(vg.fract, Anew)
    return vg
end

function _y_sweep!(vg::SplitVOFGrid, velocity, dt::Float64, t::Float64)
    dt <= 0.0 && return vg
    g = vg.grid
    work = _get_sweep2d_work(vg)

    vface = work.vface
    flux = work.fluxy
    Anew = work.Anew
    poly = work.poly

    _y_face_velocities!(vface, vg, velocity, t)
    fill!(flux, 0.0)

    @inbounds for jface in 1:(g.ny + 1), i in 1:g.nx
        v = vface[i, jface]
        if abs(v) <= eps(Float64)
            continue
        end

        jB = jface == 1 ? g.ny : (jface - 1)
        jT = jface == (g.ny + 1) ? 1 : jface

        donor_j, side, sgn = if v > 0.0
            (jB, 1, 1.0)
        else
            (jT, -1, -1.0)
        end

        delta = min(abs(v) * dt, g.dy)
        areaf = _material_area_y(vg, i, donor_j, delta, side, poly)
        flux[i, jface] = sgn * areaf
    end

    Aold = vg.fract
    invvol = 1.0 / (g.dx * g.dy)
    @inbounds for j in 1:g.ny, i in 1:g.nx
        Anew[i, j] = Aold[i, j] - (flux[i, j + 1] - flux[i, j]) * invvol
    end

    _soft_clip!(Anew, vg.tol)
    copyto!(vg.fract, Anew)
    return vg
end

function _x_sweep!(vg::SplitVOFGrid3D, velocity, dt::Float64, t::Float64)
    dt <= 0.0 && return vg
    g = vg.grid
    work = _get_sweep3d_work(vg)

    uface = work.uface
    flux = work.fluxx
    Anew = work.Anew
    poly = work.poly

    _x_face_velocities!(uface, vg, velocity, t)
    fill!(flux, 0.0)

    @inbounds for k in 1:g.nz, j in 1:g.ny, iface in 1:(g.nx + 1)
        u = uface[iface, j, k]
        if abs(u) <= eps(Float64)
            continue
        end

        iL = iface == 1 ? g.nx : (iface - 1)
        iR = iface == (g.nx + 1) ? 1 : iface

        donor_i, side, sgn = if u > 0.0
            (iL, 1, 1.0)
        else
            (iR, -1, -1.0)
        end

        delta = min(abs(u) * dt, g.dx)
        volf = _material_volume_x(vg, donor_i, j, k, delta, side, poly)
        flux[iface, j, k] = sgn * volf
    end

    Aold = vg.fract
    invvol = 1.0 / (g.dx * g.dy * g.dz)
    @inbounds for k in 1:g.nz, j in 1:g.ny, i in 1:g.nx
        Anew[i, j, k] = Aold[i, j, k] - (flux[i + 1, j, k] - flux[i, j, k]) * invvol
    end

    _soft_clip!(Anew, vg.tol)
    copyto!(vg.fract, Anew)
    return vg
end

function _y_sweep!(vg::SplitVOFGrid3D, velocity, dt::Float64, t::Float64)
    dt <= 0.0 && return vg
    g = vg.grid
    work = _get_sweep3d_work(vg)

    vface = work.vface
    flux = work.fluxy
    Anew = work.Anew
    poly = work.poly

    _y_face_velocities!(vface, vg, velocity, t)
    fill!(flux, 0.0)

    @inbounds for k in 1:g.nz, jface in 1:(g.ny + 1), i in 1:g.nx
        v = vface[i, jface, k]
        if abs(v) <= eps(Float64)
            continue
        end

        jB = jface == 1 ? g.ny : (jface - 1)
        jT = jface == (g.ny + 1) ? 1 : jface

        donor_j, side, sgn = if v > 0.0
            (jB, 1, 1.0)
        else
            (jT, -1, -1.0)
        end

        delta = min(abs(v) * dt, g.dy)
        volf = _material_volume_y(vg, i, donor_j, k, delta, side, poly)
        flux[i, jface, k] = sgn * volf
    end

    Aold = vg.fract
    invvol = 1.0 / (g.dx * g.dy * g.dz)
    @inbounds for k in 1:g.nz, j in 1:g.ny, i in 1:g.nx
        Anew[i, j, k] = Aold[i, j, k] - (flux[i, j + 1, k] - flux[i, j, k]) * invvol
    end

    _soft_clip!(Anew, vg.tol)
    copyto!(vg.fract, Anew)
    return vg
end

function _z_sweep!(vg::SplitVOFGrid3D, velocity, dt::Float64, t::Float64)
    dt <= 0.0 && return vg
    g = vg.grid
    work = _get_sweep3d_work(vg)

    wface = work.wface
    flux = work.fluxz
    Anew = work.Anew
    poly = work.poly

    _z_face_velocities!(wface, vg, velocity, t)
    fill!(flux, 0.0)

    @inbounds for kface in 1:(g.nz + 1), j in 1:g.ny, i in 1:g.nx
        w = wface[i, j, kface]
        if abs(w) <= eps(Float64)
            continue
        end

        kD = kface == 1 ? g.nz : (kface - 1)
        kU = kface == (g.nz + 1) ? 1 : kface

        donor_k, side, sgn = if w > 0.0
            (kD, 1, 1.0)
        else
            (kU, -1, -1.0)
        end

        delta = min(abs(w) * dt, g.dz)
        volf = _material_volume_z(vg, i, j, donor_k, delta, side, poly)
        flux[i, j, kface] = sgn * volf
    end

    Aold = vg.fract
    invvol = 1.0 / (g.dx * g.dy * g.dz)
    @inbounds for k in 1:g.nz, j in 1:g.ny, i in 1:g.nx
        Anew[i, j, k] = Aold[i, j, k] - (flux[i, j, k + 1] - flux[i, j, k]) * invvol
    end

    _soft_clip!(Anew, vg.tol)
    copyto!(vg.fract, Anew)
    return vg
end
