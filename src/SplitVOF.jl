module SplitVOF

using LinearAlgebra
using StaticArrays
using VOFTools

include("types.jl")
include("grid.jl")
include("init.jl")
include("reconstruction.jl")
include("flux.jl")
include("advection.jl")

export SplitVOFParams,
    ReconstructionMethod,
    AdvectionMethod,
    YoungsPLIC,
    StrangSplit,
    CartesianGrid2D,
    CartesianGrid3D,
    SplitVOFGrid,
    SplitVOFGrid3D,
    splitgrid,
    vofgrid,
    fractg,
    dimension,
    cell_center,
    cell_bounds,
    volume,
    initfgrid!,
    circle_levelset,
    sphere_levelset,
    zalesak_levelset,
    reconstruct!,
    youngs_normals!,
    is_empty,
    is_full,
    is_mixed,
    compute_dt,
    compdt,
    step!,
    vofadv!,
    integrate!,
    l1_error,
    l2_error,
    linf_error,
    shape_metrics

end # module SplitVOF
