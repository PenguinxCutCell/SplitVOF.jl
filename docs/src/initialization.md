# [Initialization](@id initialization)

## API

- `initfgrid!`
- `circle_levelset`
- `sphere_levelset`
- `zalesak_levelset`

Initialization uses geometric integration from `VOFTools.jl` (`initf2d`/`initf3d`) and is designed to be allocation-free in hot loops.
