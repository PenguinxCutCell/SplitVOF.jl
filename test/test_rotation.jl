@testset "Split advection: solid body rotation" begin
    params = SplitVOFParams(nx=40, ny=40, xlim=(-1.0, 1.0), ylim=(-1.0, 1.0), cfl=0.45)
    vg = splitgrid(params)

    initfgrid!(vg, zalesak_levelset(center=(0.0, 0.0), radius=0.35, slot_width=0.12, slot_height=0.30); nc=8)
    v0 = volume(vg)
    s0 = shape_metrics(vg)
    tr0 = tr(s0.second_moments)
    det0 = det(s0.second_moments)

    ω = 2π
    u = (x, t) -> SVector(-ω * x[2], ω * x[1])

    integrate!(vg, u, 0.2; cfl=params.cfl)

    v1 = volume(vg)
    s1 = shape_metrics(vg)
    tr1 = tr(s1.second_moments)
    det1 = det(s1.second_moments)
    @test abs(v1 - v0) / v0 < 1e-10
    @test all((-1e-12 .<= vg.fract) .& (vg.fract .<= 1 + 1e-12))
    @test vg.step > 0
    @test norm(s1.centroid - s0.centroid) < 3.0e-2
    @test abs(tr1 - tr0) / tr0 < 6.0e-2
    @test abs(det1 - det0) / max(det0, eps(Float64)) < 1.2e-1
end
