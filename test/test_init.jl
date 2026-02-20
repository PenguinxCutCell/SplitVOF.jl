@testset "Initialization" begin
    params = SplitVOFParams(nx=48, ny=48, xlim=(-1.0, 1.0), ylim=(-1.0, 1.0))
    vg = splitgrid(params)
    R = 0.35
    initfgrid!(vg, circle_levelset((0.0, 0.0), R); nc=8)

    @test all(0.0 .<= vg.fract .<= 1.0)

    area_exact = π * R^2
    @test abs(volume(vg) - area_exact) < 1.5e-2
end
