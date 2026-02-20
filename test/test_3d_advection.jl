@testset "3D split advection one-step" begin
    params = SplitVOFParams(nx=12, ny=12, nz=12,
                            xlim=(-1.0, 1.0), ylim=(-1.0, 1.0), zlim=(-1.0, 1.0),
                            cfl=0.45,
                            reconstruction=YoungsPLIC(),
                            advection=StrangSplit())
    vg = splitgrid(params)

    initfgrid!(vg, sphere_levelset((0.0, 0.0, 0.0), 0.45); nc=4)
    reconstruct!(vg, params)

    v0 = volume(vg)
    vel = (x, t) -> SVector(0.35, -0.10, 0.20)

    dt = 0.4 * compdt(vg, vel, params)
    vofadv!(vg, vel, dt, params)

    v1 = volume(vg)
    @test abs(v1 - v0) / v0 < 1e-8
    @test all((-1e-12 .<= vg.fract) .& (vg.fract .<= 1 + 1e-12))
    @test vg.step == 1
end
