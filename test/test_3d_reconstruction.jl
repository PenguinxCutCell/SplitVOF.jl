@testset "3D Youngs + PLIC reconstruction" begin
    params = SplitVOFParams(nx=14, ny=14, nz=14,
                            xlim=(-1.0, 1.0), ylim=(-1.0, 1.0), zlim=(-1.0, 1.0))
    vg = splitgrid(params)
    initfgrid!(vg, sphere_levelset((0.0, 0.0, 0.0), 0.5); nc=4)

    reconstruct!(vg)

    mixed = [(i, j, k)
             for k in 1:vg.grid.nz, j in 1:vg.grid.ny, i in 1:vg.grid.nx
             if is_mixed(vg, i, j, k)]
    @test !isempty(mixed)

    for (i, j, k) in mixed
        nmag = sqrt(vg.nxint[i, j, k]^2 + vg.nyint[i, j, k]^2 + vg.nzint[i, j, k]^2)
        @test abs(nmag - 1.0) < 1e-10
        @test isfinite(vg.cplic[i, j, k])
    end

    nsamp = min(length(mixed), 10)
    dv = vg.grid.dx * vg.grid.dy * vg.grid.dz
    for q in 1:nsamp
        i, j, k = mixed[q]
        xL, xR, yB, yT, zD, zU = cell_bounds(vg.grid, i, j, k)
        poly = SplitVOF._rect_polyhedron(xL, xR, yB, yT, zD, zU)
        _, icontp = inte3d!(poly,
                            vg.cplic[i, j, k],
                            vg.nxint[i, j, k],
                            vg.nyint[i, j, k],
                            vg.nzint[i, j, k])
        v = icontp == 0 ? 0.0 : toolv3d(poly)
        @test abs(v / dv - vg.fract[i, j, k]) < 1e-6
    end
end
