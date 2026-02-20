using SplitVOF
using StaticArrays

params = SplitVOFParams(
    nx=128,
    ny=128,
    xlim=(-0.5, 0.5),
    ylim=(-0.5, 0.5),
    cfl=0.5,
)

vg = splitgrid(params)
initfgrid!(vg, zalesak_levelset(center=(0.0, 0.0), radius=0.3, slot_width=0.1, slot_height=0.25); nc=10)

α0 = copy(vg.fract)
V0 = volume(vg)

ω = 2π
u = (x, t) -> SVector(-ω * x[2], ω * x[1])

# One full revolution around the origin
integrate!(vg, u, 1.0; cfl=params.cfl)

V1 = volume(vg)
rel_mass = abs(V1 - V0) / V0

# Normalize L1 shape error by initial slotted-disk area
area0 = sum(α0) * vg.grid.dx * vg.grid.dy
shape_err = l1_error(vg, α0) / area0

println("Relative mass error = ", rel_mass)
println("Normalized L1 shape error = ", shape_err)
