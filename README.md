# SplitVOF.jl

A geometric split-advection VOF (PLIC) implementation on periodic Cartesian grids.

## Implemented baseline

- PLIC reconstruction in mixed cells
- Youngs normals
- Geometric swept-slab fluxes from upwind donor cells
- Strang directional splitting with sweep-order alternation
- Conservative finite-volume update of volume fractions

## Current dimensions

- 2D (`nz=1`)
- 3D (`nz>1`)

## Algorithm structure (gVOF-style)

`SplitVOFParams` carries algorithm choices as type tokens:

- `reconstruction = YoungsPLIC()`
- `advection = StrangSplit()`

The current implementation supports these baseline methods only, but the API is structured so additional split/unsplit and reconstruction/advection methods can be added later.

## Geometry backend

This package reuses `VOFTools.jl` geometric routines:

- `initf2d` / `initf3d` for initialization
- `enforv2dsz` / `enforv3dsz` for PLIC volume enforcement
- `inte2d!` / `inte3d!` + `toolv2d` / `toolv3d` for geometric clipping/volume

## Quick start (3D)

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
initfgrid!(vg, sphere_levelset((0.0, 0.0, 0.0), 0.45); nc=5)

u = (x, t) -> SVector(0.2, -0.1, 0.05)
integrate!(vg, u, 0.2, params)

println("volume = ", volume(vg))
```

## Notes

- v0.1 is a robust baseline implementation intended for verification before higher-order or unsplit variants.
- Boundary conditions are currently periodic in all directions.
