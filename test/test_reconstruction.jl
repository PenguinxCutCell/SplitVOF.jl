@testset "Youngs + PLIC reconstruction" begin
    params = SplitVOFParams(nx=40, ny=40, xlim=(-1.0, 1.0), ylim=(-1.0, 1.0))
    vg = splitgrid(params)
    initfgrid!(vg, circle_levelset((0.0, 0.0), 0.4); nc=8)

    reconstruct!(vg)

    mixed = [(i, j) for j in 1:vg.grid.ny for i in 1:vg.grid.nx if is_mixed(vg, i, j)]
    @test !isempty(mixed)

    for (i, j) in mixed
        nmag = hypot(vg.nxint[i, j], vg.nyint[i, j])
        @test abs(nmag - 1.0) < 1e-10
        @test isfinite(vg.cplic[i, j])
    end

    # Check local volume enforcement for a sample of mixed cells
    nsamp = min(length(mixed), 20)
    dxdy = vg.grid.dx * vg.grid.dy
    for k in 1:nsamp
        i, j = mixed[k]
        xL, xR, yB, yT = cell_bounds(vg.grid, i, j)
        poly = SplitVOF._rect_polygon(xL, xR, yB, yT)
        _, icontp = inte2d!(poly, vg.cplic[i, j], vg.nxint[i, j], vg.nyint[i, j])
        area = icontp == 0 ? 0.0 : toolv2d(poly)
        @test abs(area / dxdy - vg.fract[i, j]) < 1e-7
    end
end
