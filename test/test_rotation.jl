@testset "Split advection: solid body rotation" begin
    params = SplitVOFParams(nx=40, ny=40, xlim=(-1.0, 1.0), ylim=(-1.0, 1.0), cfl=0.45)
    vg = splitgrid(params)

    initfgrid!(vg, zalesak_levelset(center=(0.0, 0.0), radius=0.35, slot_width=0.12, slot_height=0.30); nc=8)
    v0 = volume(vg)

    ω = 2π
    u = (x, t) -> SVector(-ω * x[2], ω * x[1])

    integrate!(vg, u, 0.2; cfl=params.cfl)

    v1 = volume(vg)
    @test abs(v1 - v0) / v0 < 1e-10
    @test all((-1e-12 .<= vg.fract) .& (vg.fract .<= 1 + 1e-12))
    @test vg.step > 0
end
