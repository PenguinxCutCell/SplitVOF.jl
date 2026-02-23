"""
Abstract supertype for interface-reconstruction algorithms.
"""
abstract type ReconstructionMethod end

"""
Youngs gradient-based normal reconstruction with PLIC volume enforcement.
"""
struct YoungsPLIC <: ReconstructionMethod end

"""
Abstract supertype for advection/integration algorithms.
"""
abstract type AdvectionMethod end

"""
Directional Strang splitting advection scheme.
"""
struct StrangSplit <: AdvectionMethod end

"""
Simulation parameters for split geometric VOF on a Cartesian grid.

Key fields:
- `nx, ny, nz`: number of cells
- `xlim, ylim, zlim`: domain bounds
- `cfl`: target CFL for `compute_dt`/`integrate!`
- `tol`: interface tolerance for empty/full detection
- `reconstruction`, `advection`: algorithm selectors
"""
Base.@kwdef struct SplitVOFParams
    nx::Int = 64
    ny::Int = 64
    nz::Int = 1
    xlim::Tuple{Float64,Float64} = (0.0, 1.0)
    ylim::Tuple{Float64,Float64} = (0.0, 1.0)
    zlim::Tuple{Float64,Float64} = (0.0, 1.0)
    cfl::Float64 = 0.5
    tol::Float64 = 1.0e-12
    reconstruction::ReconstructionMethod = YoungsPLIC()
    advection::AdvectionMethod = StrangSplit()
end

"""
Uniform 2D Cartesian grid descriptor.
"""
struct CartesianGrid2D{T<:AbstractFloat}
    xlo::T
    xhi::T
    ylo::T
    yhi::T
    nx::Int
    ny::Int
    dx::T
    dy::T
end

"""
Build a 2D Cartesian grid from domain limits and cell counts.
"""
function CartesianGrid2D(xlim::NTuple{2,<:Real},
                         ylim::NTuple{2,<:Real},
                         n::NTuple{2,Int})
    nx, ny = n
    nx > 0 || throw(ArgumentError("nx must be > 0"))
    ny > 0 || throw(ArgumentError("ny must be > 0"))

    xlo, xhi = Float64(xlim[1]), Float64(xlim[2])
    ylo, yhi = Float64(ylim[1]), Float64(ylim[2])
    xhi > xlo || throw(ArgumentError("xhi must be > xlo"))
    yhi > ylo || throw(ArgumentError("yhi must be > ylo"))

    dx = (xhi - xlo) / nx
    dy = (yhi - ylo) / ny
    return CartesianGrid2D(xlo, xhi, ylo, yhi, nx, ny, dx, dy)
end

"""
Uniform 3D Cartesian grid descriptor.
"""
struct CartesianGrid3D{T<:AbstractFloat}
    xlo::T
    xhi::T
    ylo::T
    yhi::T
    zlo::T
    zhi::T
    nx::Int
    ny::Int
    nz::Int
    dx::T
    dy::T
    dz::T
end

"""
Build a 3D Cartesian grid from domain limits and cell counts.
"""
function CartesianGrid3D(xlim::NTuple{2,<:Real},
                         ylim::NTuple{2,<:Real},
                         zlim::NTuple{2,<:Real},
                         n::NTuple{3,Int})
    nx, ny, nz = n
    nx > 0 || throw(ArgumentError("nx must be > 0"))
    ny > 0 || throw(ArgumentError("ny must be > 0"))
    nz > 0 || throw(ArgumentError("nz must be > 0"))

    xlo, xhi = Float64(xlim[1]), Float64(xlim[2])
    ylo, yhi = Float64(ylim[1]), Float64(ylim[2])
    zlo, zhi = Float64(zlim[1]), Float64(zlim[2])

    xhi > xlo || throw(ArgumentError("xhi must be > xlo"))
    yhi > ylo || throw(ArgumentError("yhi must be > ylo"))
    zhi > zlo || throw(ArgumentError("zhi must be > zlo"))

    dx = (xhi - xlo) / nx
    dy = (yhi - ylo) / ny
    dz = (zhi - zlo) / nz
    return CartesianGrid3D(xlo, xhi, ylo, yhi, zlo, zhi, nx, ny, nz, dx, dy, dz)
end

"""
2D split-VOF state container.
"""
mutable struct SplitVOFGrid{T<:AbstractFloat}
    grid::CartesianGrid2D{T}
    fract::Matrix{T}
    nxint::Matrix{T}
    nyint::Matrix{T}
    cplic::Matrix{T}
    t::T
    step::Int
    tol::T
end

"""
3D split-VOF state container.
"""
mutable struct SplitVOFGrid3D{T<:AbstractFloat}
    grid::CartesianGrid3D{T}
    fract::Array{T,3}
    nxint::Array{T,3}
    nyint::Array{T,3}
    nzint::Array{T,3}
    cplic::Array{T,3}
    t::T
    step::Int
    tol::T
end

const AbstractSplitVOFGrid = Union{SplitVOFGrid, SplitVOFGrid3D}

"""
Return spatial dimension (`2` or `3`) of a split-VOF state.
"""
dimension(::SplitVOFGrid) = 2
dimension(::SplitVOFGrid3D) = 3

"""
Allocate and initialize a split-VOF state from `SplitVOFParams`.

Returns `SplitVOFGrid` when `nz <= 1`, otherwise `SplitVOFGrid3D`.
"""
function splitgrid(params::SplitVOFParams=SplitVOFParams())
    if params.nz <= 1
        grid = CartesianGrid2D(params.xlim, params.ylim, (params.nx, params.ny))
        fract = zeros(Float64, params.nx, params.ny)
        nxint = zeros(Float64, params.nx, params.ny)
        nyint = zeros(Float64, params.nx, params.ny)
        cplic = zeros(Float64, params.nx, params.ny)
        return SplitVOFGrid(grid, fract, nxint, nyint, cplic, 0.0, 0, params.tol)
    else
        grid = CartesianGrid3D(params.xlim, params.ylim, params.zlim,
                               (params.nx, params.ny, params.nz))
        fract = zeros(Float64, params.nx, params.ny, params.nz)
        nxint = zeros(Float64, params.nx, params.ny, params.nz)
        nyint = zeros(Float64, params.nx, params.ny, params.nz)
        nzint = zeros(Float64, params.nx, params.ny, params.nz)
        cplic = zeros(Float64, params.nx, params.ny, params.nz)
        return SplitVOFGrid3D(grid, fract, nxint, nyint, nzint, cplic, 0.0, 0, params.tol)
    end
end

"""
Allocate a 2D split-VOF state from an existing Cartesian grid.
"""
function splitgrid(grid::CartesianGrid2D; tol::Real=1.0e-12)
    fract = zeros(Float64, grid.nx, grid.ny)
    nxint = zeros(Float64, grid.nx, grid.ny)
    nyint = zeros(Float64, grid.nx, grid.ny)
    cplic = zeros(Float64, grid.nx, grid.ny)
    return SplitVOFGrid(grid, fract, nxint, nyint, cplic, 0.0, 0, Float64(tol))
end

"""
Allocate a 3D split-VOF state from an existing Cartesian grid.
"""
function splitgrid(grid::CartesianGrid3D; tol::Real=1.0e-12)
    fract = zeros(Float64, grid.nx, grid.ny, grid.nz)
    nxint = zeros(Float64, grid.nx, grid.ny, grid.nz)
    nyint = zeros(Float64, grid.nx, grid.ny, grid.nz)
    nzint = zeros(Float64, grid.nx, grid.ny, grid.nz)
    cplic = zeros(Float64, grid.nx, grid.ny, grid.nz)
    return SplitVOFGrid3D(grid, fract, nxint, nyint, nzint, cplic, 0.0, 0, Float64(tol))
end

"""
Alias for `splitgrid(params)`.
"""
vofgrid(params::SplitVOFParams=SplitVOFParams()) = splitgrid(params)

"""
Return the volume-fraction array stored in the split-VOF state.
"""
fractg(vg::AbstractSplitVOFGrid) = vg.fract

Base.size(vg::SplitVOFGrid) = size(vg.fract)
Base.size(vg::SplitVOFGrid3D) = size(vg.fract)

Base.getindex(vg::SplitVOFGrid, i::Int, j::Int) = vg.fract[i, j]
Base.getindex(vg::SplitVOFGrid3D, i::Int, j::Int, k::Int) = vg.fract[i, j, k]

function Base.setindex!(vg::SplitVOFGrid, v, i::Int, j::Int)
    vg.fract[i, j] = v
    return vg
end

function Base.setindex!(vg::SplitVOFGrid3D, v, i::Int, j::Int, k::Int)
    vg.fract[i, j, k] = v
    return vg
end

"""
Total material volume (area in 2D) represented by the volume-fraction field.
"""
volume(vg::SplitVOFGrid) = sum(vg.fract) * vg.grid.dx * vg.grid.dy
volume(vg::SplitVOFGrid3D) = sum(vg.fract) * vg.grid.dx * vg.grid.dy * vg.grid.dz

"""
L1 error of the volume-fraction field against a reference field.
"""
function l1_error(vg::SplitVOFGrid, ref::AbstractMatrix{<:Real})
    size(ref) == size(vg.fract) || throw(ArgumentError("reference has wrong size"))
    return sum(abs.(vg.fract .- ref)) * vg.grid.dx * vg.grid.dy
end

function l1_error(vg::SplitVOFGrid3D, ref::AbstractArray{<:Real,3})
    size(ref) == size(vg.fract) || throw(ArgumentError("reference has wrong size"))
    return sum(abs.(vg.fract .- ref)) * vg.grid.dx * vg.grid.dy * vg.grid.dz
end

"""
L2 error of the volume-fraction field against a reference field.
"""
function l2_error(vg::SplitVOFGrid, ref::AbstractMatrix{<:Real})
    size(ref) == size(vg.fract) || throw(ArgumentError("reference has wrong size"))
    s = 0.0
    @inbounds for j in 1:vg.grid.ny, i in 1:vg.grid.nx
        d = vg.fract[i, j] - ref[i, j]
        s += d * d
    end
    return sqrt(s * vg.grid.dx * vg.grid.dy)
end

function l2_error(vg::SplitVOFGrid3D, ref::AbstractArray{<:Real,3})
    size(ref) == size(vg.fract) || throw(ArgumentError("reference has wrong size"))
    s = 0.0
    @inbounds for k in 1:vg.grid.nz, j in 1:vg.grid.ny, i in 1:vg.grid.nx
        d = vg.fract[i, j, k] - ref[i, j, k]
        s += d * d
    end
    return sqrt(s * vg.grid.dx * vg.grid.dy * vg.grid.dz)
end

"""
L∞ error of the volume-fraction field against a reference field.
"""
function linf_error(vg::SplitVOFGrid, ref::AbstractMatrix{<:Real})
    size(ref) == size(vg.fract) || throw(ArgumentError("reference has wrong size"))
    e = 0.0
    @inbounds for j in 1:vg.grid.ny, i in 1:vg.grid.nx
        e = max(e, abs(vg.fract[i, j] - ref[i, j]))
    end
    return e
end

function linf_error(vg::SplitVOFGrid3D, ref::AbstractArray{<:Real,3})
    size(ref) == size(vg.fract) || throw(ArgumentError("reference has wrong size"))
    e = 0.0
    @inbounds for k in 1:vg.grid.nz, j in 1:vg.grid.ny, i in 1:vg.grid.nx
        e = max(e, abs(vg.fract[i, j, k] - ref[i, j, k]))
    end
    return e
end

"""
Compute bulk shape descriptors from the volume-fraction field.

Returns a named tuple:
- `volume`
- `centroid`
- `second_moments` (central second-moment tensor)
"""
function shape_metrics(vg::SplitVOFGrid)
    g = vg.grid
    dA = g.dx * g.dy
    vol = 0.0
    sx = 0.0
    sy = 0.0
    @inbounds for j in 1:g.ny, i in 1:g.nx
        α = vg.fract[i, j]
        x = g.xlo + (i - 0.5) * g.dx
        y = g.ylo + (j - 0.5) * g.dy
        w = α * dA
        vol += w
        sx += w * x
        sy += w * y
    end

    if vol <= eps(Float64)
        return (volume=0.0,
                centroid=SVector{2,Float64}(0.0, 0.0),
                second_moments=SMatrix{2,2,Float64,4}(0.0, 0.0, 0.0, 0.0))
    end

    cx = sx / vol
    cy = sy / vol
    mxx = 0.0
    mxy = 0.0
    myy = 0.0
    @inbounds for j in 1:g.ny, i in 1:g.nx
        α = vg.fract[i, j]
        x = g.xlo + (i - 0.5) * g.dx
        y = g.ylo + (j - 0.5) * g.dy
        dx = x - cx
        dy = y - cy
        w = α * dA
        mxx += w * dx * dx
        mxy += w * dx * dy
        myy += w * dy * dy
    end

    M2 = SMatrix{2,2,Float64,4}(mxx, mxy, mxy, myy)
    return (volume=vol, centroid=SVector{2,Float64}(cx, cy), second_moments=M2)
end

function shape_metrics(vg::SplitVOFGrid3D)
    g = vg.grid
    dV = g.dx * g.dy * g.dz
    vol = 0.0
    sx = 0.0
    sy = 0.0
    sz = 0.0
    @inbounds for k in 1:g.nz, j in 1:g.ny, i in 1:g.nx
        α = vg.fract[i, j, k]
        x = g.xlo + (i - 0.5) * g.dx
        y = g.ylo + (j - 0.5) * g.dy
        z = g.zlo + (k - 0.5) * g.dz
        w = α * dV
        vol += w
        sx += w * x
        sy += w * y
        sz += w * z
    end

    if vol <= eps(Float64)
        return (volume=0.0,
                centroid=SVector{3,Float64}(0.0, 0.0, 0.0),
                second_moments=SMatrix{3,3,Float64,9}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    end

    cx = sx / vol
    cy = sy / vol
    cz = sz / vol
    mxx = 0.0; mxy = 0.0; mxz = 0.0
    myy = 0.0; myz = 0.0; mzz = 0.0
    @inbounds for k in 1:g.nz, j in 1:g.ny, i in 1:g.nx
        α = vg.fract[i, j, k]
        x = g.xlo + (i - 0.5) * g.dx
        y = g.ylo + (j - 0.5) * g.dy
        z = g.zlo + (k - 0.5) * g.dz
        dx = x - cx
        dy = y - cy
        dz = z - cz
        w = α * dV
        mxx += w * dx * dx
        mxy += w * dx * dy
        mxz += w * dx * dz
        myy += w * dy * dy
        myz += w * dy * dz
        mzz += w * dz * dz
    end

    M2 = SMatrix{3,3,Float64,9}(mxx, mxy, mxz,
                                mxy, myy, myz,
                                mxz, myz, mzz)
    return (volume=vol, centroid=SVector{3,Float64}(cx, cy, cz), second_moments=M2)
end
