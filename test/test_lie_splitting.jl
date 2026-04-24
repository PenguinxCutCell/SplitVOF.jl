@testset "Lie split advection: 2D periodic translation" begin
    params = SplitVOFParams(nx=48, ny=48,
                            xlim=(-1.0, 1.0), ylim=(-1.0, 1.0),
                            cfl=0.5,
                            advection=LieSplit())
    vg = splitgrid(params)

    R = 0.30
    initfgrid!(vg, circle_levelset((0.0, 0.0), R); nc=8)

    alpha0 = copy(vg.fract)
    v0 = volume(vg)
    s0 = shape_metrics(vg)

    vel = (x, t) -> SVector(0.5, 0.0)
    tf = 4.0
    integrate!(vg, vel, tf, params)

    v1 = volume(vg)
    s1 = shape_metrics(vg)
    rel_shape_err = l1_error(vg, alpha0) / (π * R^2)

    @test abs(v1 - v0) / v0 < 1e-12
    @test all((-1e-12 .<= vg.fract) .& (vg.fract .<= 1 + 1e-12))
    @test rel_shape_err < 0.11
    @test norm(s1.centroid - s0.centroid) < 2.5e-2
    @test norm(s1.second_moments - s0.second_moments) / norm(s0.second_moments) < 8.0e-2
    @test vg.step > 0
end

@testset "Lie split advection: 3D one-step translation" begin
    params = SplitVOFParams(nx=12, ny=12, nz=12,
                            xlim=(-1.0, 1.0), ylim=(-1.0, 1.0), zlim=(-1.0, 1.0),
                            cfl=0.45,
                            reconstruction=YoungsPLIC(),
                            advection=LieSplit())
    vg = splitgrid(params)

    initfgrid!(vg, sphere_levelset((0.0, 0.0, 0.0), 0.45); nc=4)
    reconstruct!(vg, params)

    v0 = volume(vg)
    s0 = shape_metrics(vg)
    vel = (x, t) -> SVector(0.35, -0.10, 0.20)

    dt = 0.4 * compdt(vg, vel, params)
    vofadv!(vg, vel, dt, params)

    v1 = volume(vg)
    s1 = shape_metrics(vg)

    @test abs(v1 - v0) / v0 < 1e-8
    @test all((-1e-12 .<= vg.fract) .& (vg.fract .<= 1 + 1e-12))
    @test vg.step == 1
    @test norm((s1.centroid - s0.centroid) - vel(SVector(0.0, 0.0, 0.0), 0.0) * dt) < 4.0e-2
    @test norm(s1.second_moments - s0.second_moments) / norm(s0.second_moments) < 1.0e-1
end
