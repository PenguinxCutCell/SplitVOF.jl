@testset "Split advection: translation" begin
    params = SplitVOFParams(nx=48, ny=48, xlim=(-1.0, 1.0), ylim=(-1.0, 1.0), cfl=0.5)
    vg = splitgrid(params)
    R = 0.30
    initfgrid!(vg, circle_levelset((0.0, 0.0), R); nc=8)

    alpha0 = copy(vg.fract)
    v0 = volume(vg)

    u = (x, t) -> SVector(0.5, 0.0)
    tf = 4.0  # one full periodic lap in x over domain length 2.0
    integrate!(vg, u, tf; cfl=params.cfl)

    v1 = volume(vg)
    @test abs(v1 - v0) / v0 < 1e-12
    @test all((-1e-12 .<= vg.fract) .& (vg.fract .<= 1 + 1e-12))

    rel_shape_err = l1_error(vg, alpha0) / (π * R^2)
    @test rel_shape_err < 0.08
end
