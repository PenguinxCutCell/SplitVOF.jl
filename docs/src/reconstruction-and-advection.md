# [Reconstruction And Advection](@id reconstruction-and-advection)

## API

- `reconstruct!`, `youngs_normals!`
- `is_empty`, `is_full`, `is_mixed`
- `compute_dt`, `compdt`
- `step!`, `vofadv!`, `integrate!`

Current baseline:

- Youngs normal estimation
- PLIC reconstruction
- Directionally split Strang advection
- Geometric swept-slab fluxes
