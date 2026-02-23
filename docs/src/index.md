```@meta
CurrentModule = SplitVOF
```

# SplitVOF.jl

`SplitVOF.jl` implements a split geometric VOF (PLIC) advection method on periodic Cartesian grids in 2D/3D.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/PenguinxCutCell/SplitVOF.jl")
```

## Quick Start

```julia
using SplitVOF
using StaticArrays

params = SplitVOFParams(
    nx=32, ny=32, nz=32,
    xlim=(-1.0, 1.0), ylim=(-1.0, 1.0), zlim=(-1.0, 1.0),
    cfl=0.5,
    reconstruction=YoungsPLIC(),
    advection=StrangSplit(),
)

vg = vofgrid(params)
initfgrid!(vg, sphere_levelset((0.0, 0.0, 0.0), 0.45); nc=6)

u = (x, t) -> SVector(0.25, -0.10, 0.15)
integrate!(vg, u, 0.2, params)

println("volume = ", volume(vg))
println("shape = ", shape_metrics(vg))
```

## Package Structure

- [Grid And Params](@ref grid-and-params)
- [Initialization](@ref initialization)
- [Reconstruction And Advection](@ref reconstruction-and-advection)
- [Validation And Metrics](@ref validation-and-metrics)
- [Reference](@ref reference)
