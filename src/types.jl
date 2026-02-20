abstract type ReconstructionMethod end
struct YoungsPLIC <: ReconstructionMethod end

abstract type AdvectionMethod end
struct StrangSplit <: AdvectionMethod end

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

dimension(::SplitVOFGrid) = 2
dimension(::SplitVOFGrid3D) = 3

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

function splitgrid(grid::CartesianGrid2D; tol::Real=1.0e-12)
    fract = zeros(Float64, grid.nx, grid.ny)
    nxint = zeros(Float64, grid.nx, grid.ny)
    nyint = zeros(Float64, grid.nx, grid.ny)
    cplic = zeros(Float64, grid.nx, grid.ny)
    return SplitVOFGrid(grid, fract, nxint, nyint, cplic, 0.0, 0, Float64(tol))
end

function splitgrid(grid::CartesianGrid3D; tol::Real=1.0e-12)
    fract = zeros(Float64, grid.nx, grid.ny, grid.nz)
    nxint = zeros(Float64, grid.nx, grid.ny, grid.nz)
    nyint = zeros(Float64, grid.nx, grid.ny, grid.nz)
    nzint = zeros(Float64, grid.nx, grid.ny, grid.nz)
    cplic = zeros(Float64, grid.nx, grid.ny, grid.nz)
    return SplitVOFGrid3D(grid, fract, nxint, nyint, nzint, cplic, 0.0, 0, Float64(tol))
end

vofgrid(params::SplitVOFParams=SplitVOFParams()) = splitgrid(params)
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

volume(vg::SplitVOFGrid) = sum(vg.fract) * vg.grid.dx * vg.grid.dy
volume(vg::SplitVOFGrid3D) = sum(vg.fract) * vg.grid.dx * vg.grid.dy * vg.grid.dz

function l1_error(vg::SplitVOFGrid, ref::AbstractMatrix{<:Real})
    size(ref) == size(vg.fract) || throw(ArgumentError("reference has wrong size"))
    return sum(abs.(vg.fract .- ref)) * vg.grid.dx * vg.grid.dy
end

function l1_error(vg::SplitVOFGrid3D, ref::AbstractArray{<:Real,3})
    size(ref) == size(vg.fract) || throw(ArgumentError("reference has wrong size"))
    return sum(abs.(vg.fract .- ref)) * vg.grid.dx * vg.grid.dy * vg.grid.dz
end
