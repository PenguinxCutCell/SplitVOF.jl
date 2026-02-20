@testset "3D initialization" begin
    params = SplitVOFParams(nx=20, ny=20, nz=20,
                            xlim=(-1.0, 1.0), ylim=(-1.0, 1.0), zlim=(-1.0, 1.0))
    vg = splitgrid(params)

    R = 0.45
    initfgrid!(vg, sphere_levelset((0.0, 0.0, 0.0), R); nc=5)

    @test dimension(vg) == 3
    @test all(0.0 .<= vg.fract .<= 1.0)

    vol_exact = 4.0 / 3.0 * π * R^3
    @test abs(volume(vg) - vol_exact) < 5.0e-2
end
