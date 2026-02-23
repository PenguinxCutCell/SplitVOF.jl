function _rect_polygon(xL::Float64, xR::Float64, yB::Float64, yT::Float64)
    nvmax = 16
    vertp = zeros(Float64, nvmax, 2)
    ipv = zeros(Int, nvmax)

    poly = Polygon2D(vertp, ipv, 4, 4)
    _rect_polygon!(poly, xL, xR, yB, yT)
    return poly
end

function _rect_polygon!(poly::Polygon2D, xL::Float64, xR::Float64, yB::Float64, yT::Float64)
    poly.ntp = 4
    poly.ntv = 4
    ipv = poly.ipv
    vertp = poly.vertp

    ipv[1] = 1; ipv[2] = 2; ipv[3] = 3; ipv[4] = 4
    vertp[1, 1] = xL; vertp[1, 2] = yB
    vertp[2, 1] = xR; vertp[2, 2] = yB
    vertp[3, 1] = xR; vertp[3, 2] = yT
    vertp[4, 1] = xL; vertp[4, 2] = yT
    return poly
end

function _rect_polyhedron(xL::Float64, xR::Float64,
                          yB::Float64, yT::Float64,
                          zD::Float64, zU::Float64)
    nsmax = 64
    nvmax = 96

    ipv = zeros(Int, nsmax, nvmax)
    nipv = zeros(Int, nsmax)
    vertp = zeros(Float64, nvmax, 3)
    xns = zeros(Float64, nsmax)
    yns = zeros(Float64, nsmax)
    zns = zeros(Float64, nsmax)

    # Face normals match VOFTools.cubicmesh orientation.
    xns[1] =  1.0; yns[1] =  0.0; zns[1] =  0.0
    xns[2] =  0.0; yns[2] = -1.0; zns[2] =  0.0
    xns[3] =  0.0; yns[3] =  0.0; zns[3] = -1.0
    xns[4] =  0.0; yns[4] =  1.0; zns[4] =  0.0
    xns[5] =  0.0; yns[5] =  0.0; zns[5] =  1.0
    xns[6] = -1.0; yns[6] =  0.0; zns[6] =  0.0

    for is in 1:6
        nipv[is] = 4
    end

    ipv[1, 1] = 1; ipv[1, 2] = 2; ipv[1, 3] = 3; ipv[1, 4] = 4
    ipv[2, 1] = 2; ipv[2, 2] = 1; ipv[2, 3] = 5; ipv[2, 4] = 6
    ipv[3, 1] = 3; ipv[3, 2] = 2; ipv[3, 3] = 6; ipv[3, 4] = 7
    ipv[4, 1] = 4; ipv[4, 2] = 3; ipv[4, 3] = 7; ipv[4, 4] = 8
    ipv[5, 1] = 1; ipv[5, 2] = 4; ipv[5, 3] = 8; ipv[5, 4] = 5
    ipv[6, 1] = 6; ipv[6, 2] = 5; ipv[6, 3] = 8; ipv[6, 4] = 7

    poly = Polyhedron3D(vertp, ipv, nipv, xns, yns, zns, 6, 8, 8)
    _rect_polyhedron!(poly, xL, xR, yB, yT, zD, zU)
    return poly
end

function _rect_polyhedron!(poly::Polyhedron3D,
                           xL::Float64, xR::Float64,
                           yB::Float64, yT::Float64,
                           zD::Float64, zU::Float64)
    poly.nts = 6
    poly.ntp = 8
    poly.ntv = 8

    nipv = poly.nipv
    ipv = poly.ipv
    xns = poly.xns
    yns = poly.yns
    zns = poly.zns
    vertp = poly.vertp

    xns[1] =  1.0; yns[1] =  0.0; zns[1] =  0.0
    xns[2] =  0.0; yns[2] = -1.0; zns[2] =  0.0
    xns[3] =  0.0; yns[3] =  0.0; zns[3] = -1.0
    xns[4] =  0.0; yns[4] =  1.0; zns[4] =  0.0
    xns[5] =  0.0; yns[5] =  0.0; zns[5] =  1.0
    xns[6] = -1.0; yns[6] =  0.0; zns[6] =  0.0

    for is in 1:6
        nipv[is] = 4
    end

    ipv[1, 1] = 1; ipv[1, 2] = 2; ipv[1, 3] = 3; ipv[1, 4] = 4
    ipv[2, 1] = 2; ipv[2, 2] = 1; ipv[2, 3] = 5; ipv[2, 4] = 6
    ipv[3, 1] = 3; ipv[3, 2] = 2; ipv[3, 3] = 6; ipv[3, 4] = 7
    ipv[4, 1] = 4; ipv[4, 2] = 3; ipv[4, 3] = 7; ipv[4, 4] = 8
    ipv[5, 1] = 1; ipv[5, 2] = 4; ipv[5, 3] = 8; ipv[5, 4] = 5
    ipv[6, 1] = 6; ipv[6, 2] = 5; ipv[6, 3] = 8; ipv[6, 4] = 7

    vertp[1, 1] = xR; vertp[1, 2] = yB; vertp[1, 3] = zU
    vertp[2, 1] = xR; vertp[2, 2] = yB; vertp[2, 3] = zD
    vertp[3, 1] = xR; vertp[3, 2] = yT; vertp[3, 3] = zD
    vertp[4, 1] = xR; vertp[4, 2] = yT; vertp[4, 3] = zU
    vertp[5, 1] = xL; vertp[5, 2] = yB; vertp[5, 3] = zU
    vertp[6, 1] = xL; vertp[6, 2] = yB; vertp[6, 3] = zD
    vertp[7, 1] = xL; vertp[7, 2] = yT; vertp[7, 3] = zD
    vertp[8, 1] = xL; vertp[8, 2] = yT; vertp[8, 3] = zU
    return poly
end

struct _Phi2DXY{F}
    f::F
end

@inline (p::_Phi2DXY)(x::Float64, y::Float64) = Float64(p.f(x, y))

struct _Phi2DSVec{F}
    f::F
end

@inline (p::_Phi2DSVec)(x::Float64, y::Float64) = Float64(p.f(SVector{2,Float64}(x, y)))

struct _Phi3DXYZ{F}
    f::F
end

@inline (p::_Phi3DXYZ)(x::Float64, y::Float64, z::Float64) = Float64(p.f(x, y, z))

struct _Phi3DSVec{F}
    f::F
end

@inline (p::_Phi3DSVec)(x::Float64, y::Float64, z::Float64) =
    Float64(p.f(SVector{3,Float64}(x, y, z)))

function _as_phi2d(f::Function)
    if applicable(f, 0.0, 0.0)
        return _Phi2DXY(f)
    elseif applicable(f, SVector{2,Float64}(0.0, 0.0))
        return _Phi2DSVec(f)
    else
        throw(ArgumentError("level-set function must accept (x, y) or SVector{2}"))
    end
end

function _as_phi3d(f::Function)
    if applicable(f, 0.0, 0.0, 0.0)
        return _Phi3DXYZ(f)
    elseif applicable(f, SVector{3,Float64}(0.0, 0.0, 0.0))
        return _Phi3DSVec(f)
    else
        throw(ArgumentError("level-set function must accept (x, y, z) or SVector{3}"))
    end
end

function _get_init_poly2d(vg::SplitVOFGrid)
    tls = task_local_storage()
    cache_any = get(tls, :_splitvof_init_poly2d, nothing)
    cache = if cache_any === nothing
        c = Dict{UInt, Polygon2D}()
        tls[:_splitvof_init_poly2d] = c
        c
    else
        cache_any::Dict{UInt, Polygon2D}
    end

    key = objectid(vg)
    poly = get(cache, key, nothing)
    if poly === nothing
        poly = _rect_polygon(0.0, 1.0, 0.0, 1.0)
        cache[key] = poly
    end
    return poly
end

function _get_init_poly3d(vg::SplitVOFGrid3D)
    tls = task_local_storage()
    cache_any = get(tls, :_splitvof_init_poly3d, nothing)
    cache = if cache_any === nothing
        c = Dict{UInt, Polyhedron3D}()
        tls[:_splitvof_init_poly3d] = c
        c
    else
        cache_any::Dict{UInt, Polyhedron3D}
    end

    key = objectid(vg)
    poly = get(cache, key, nothing)
    if poly === nothing
        poly = _rect_polyhedron(0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
        cache[key] = poly
    end
    return poly
end

function _initfgrid_phi!(vg::SplitVOFGrid, phi, nc::Int, tol::Float64)
    g = vg.grid
    poly = _get_init_poly2d(vg)

    for j in 1:g.ny, i in 1:g.nx
        xL = g.xlo + (i - 1) * g.dx
        xR = xL + g.dx
        yB = g.ylo + (j - 1) * g.dy
        yT = yB + g.dy
        _rect_polygon!(poly, xL, xR, yB, yT)
        α = initf2d(phi, poly, nc, tol)
        vg.fract[i, j] = clamp(α, 0.0, 1.0)
    end

    vg.t = 0.0
    vg.step = 0
    return vg
end

function _initfgrid_phi!(vg::SplitVOFGrid3D, phi, nc::Int, tol::Float64)
    g = vg.grid
    poly = _get_init_poly3d(vg)

    for k in 1:g.nz, j in 1:g.ny, i in 1:g.nx
        xL = g.xlo + (i - 1) * g.dx
        xR = xL + g.dx
        yB = g.ylo + (j - 1) * g.dy
        yT = yB + g.dy
        zD = g.zlo + (k - 1) * g.dz
        zU = zD + g.dz
        _rect_polyhedron!(poly, xL, xR, yB, yT, zD, zU)
        α = initf3d(phi, poly, nc, tol)
        vg.fract[i, j, k] = clamp(α, 0.0, 1.0)
    end

    vg.t = 0.0
    vg.step = 0
    return vg
end

"""
Initialize 2D volume fractions from an implicit level-set function.

`f` can accept either `(x, y)` or `SVector{2}` and is interpreted as:
- `f > 0`: inside liquid
- `f < 0`: outside liquid
"""
function initfgrid!(vg::SplitVOFGrid, f::Function; nc::Int=8, tol::Float64=10.0)
    return _initfgrid_phi!(vg, _as_phi2d(f), nc, tol)
end

"""
Initialize 3D volume fractions from an implicit level-set function.

`f` can accept either `(x, y, z)` or `SVector{3}` and is interpreted as:
- `f > 0`: inside liquid
- `f < 0`: outside liquid
"""
function initfgrid!(vg::SplitVOFGrid3D, f::Function; nc::Int=6, tol::Float64=10.0)
    return _initfgrid_phi!(vg, _as_phi3d(f), nc, tol)
end

"""
Build a circular implicit function for 2D initialization.
"""
function circle_levelset(center::NTuple{2,<:Real}=(0.5, 0.5), radius::Real=0.25)
    cx, cy = Float64(center[1]), Float64(center[2])
    r2 = Float64(radius)^2
    return (x, y) -> r2 - ((x - cx)^2 + (y - cy)^2)
end

"""
Build a spherical implicit function for 3D initialization.
"""
function sphere_levelset(center::NTuple{3,<:Real}=(0.5, 0.5, 0.5), radius::Real=0.25)
    cx, cy, cz = Float64(center[1]), Float64(center[2]), Float64(center[3])
    r2 = Float64(radius)^2
    return (x, y, z) -> r2 - ((x - cx)^2 + (y - cy)^2 + (z - cz)^2)
end

"""
Build the standard Zalesak slotted-disk implicit function.
"""
function zalesak_levelset(; center=(0.5, 0.5),
                           radius=0.3,
                           slot_width=0.1,
                           slot_height=0.25)
    cx, cy = Float64(center[1]), Float64(center[2])
    r2 = Float64(radius)^2

    x0 = cx - Float64(slot_width) / 2
    x1 = cx + Float64(slot_width) / 2
    y0 = cy - Float64(radius)
    y1 = y0 + Float64(slot_height)

    return function (x, y)
        ϕdisk = r2 - ((x - cx)^2 + (y - cy)^2)
        ϕslot = min(x - x0, x1 - x, y - y0, y1 - y)
        return min(ϕdisk, -ϕslot)
    end
end
