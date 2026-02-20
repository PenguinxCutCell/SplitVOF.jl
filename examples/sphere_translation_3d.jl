using SplitVOF
using StaticArrays

params = SplitVOFParams(
    nx=32, ny=32, nz=32,
    xlim=(-1.0, 1.0), ylim=(-1.0, 1.0), zlim=(-1.0, 1.0),
    cfl=0.5,
)

vg = vofgrid(params)
initfgrid!(vg, sphere_levelset((0.0, 0.0, 0.0), 0.45); nc=5)

v0 = volume(vg)

vel = (x, t) -> SVector(0.25, -0.05, 0.10)
integrate!(vg, vel, 0.1, params)

v1 = volume(vg)
println("relative mass error = ", abs(v1 - v0) / v0)
