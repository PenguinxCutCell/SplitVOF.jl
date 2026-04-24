"""
Return `true` if cell is empty (`α <= tol`).
"""
@inline is_empty(vg::SplitVOFGrid, i::Int, j::Int) = vg.fract[i, j] <= vg.tol

"""
Return `true` if cell is full (`α >= 1 - tol`).
"""
@inline is_full(vg::SplitVOFGrid, i::Int, j::Int) = vg.fract[i, j] >= 1.0 - vg.tol

"""
Return `true` if cell is mixed (`tol < α < 1 - tol`).
"""
@inline is_mixed(vg::SplitVOFGrid, i::Int, j::Int) = !is_empty(vg, i, j) && !is_full(vg, i, j)

@inline is_empty(vg::SplitVOFGrid3D, i::Int, j::Int, k::Int) = vg.fract[i, j, k] <= vg.tol
@inline is_full(vg::SplitVOFGrid3D, i::Int, j::Int, k::Int) = vg.fract[i, j, k] >= 1.0 - vg.tol
@inline is_mixed(vg::SplitVOFGrid3D, i::Int, j::Int, k::Int) = !is_empty(vg, i, j, k) && !is_full(vg, i, j, k)

function _reconstruct2d_verts()
    tls = task_local_storage()
    verts_any = get(tls, :_splitvof_reconstruct2d_verts, nothing)
    if verts_any === nothing
        verts = Matrix{Float64}(undef, 4, 2)
        tls[:_splitvof_reconstruct2d_verts] = verts
        return verts
    end
    return verts_any::Matrix{Float64}
end

function _reconstruct3d_verts()
    tls = task_local_storage()
    verts_any = get(tls, :_splitvof_reconstruct3d_verts, nothing)
    if verts_any === nothing
        verts = Matrix{Float64}(undef, 8, 3)
        tls[:_splitvof_reconstruct3d_verts] = verts
        return verts
    end
    return verts_any::Matrix{Float64}
end

function _youngs_gradient(vg::SplitVOFGrid, i::Int, j::Int)
    g = vg.grid
    C = vg.fract

    im = _wrap(i - 1, g.nx)
    ip = _wrap(i + 1, g.nx)
    jm = _wrap(j - 1, g.ny)
    jp = _wrap(j + 1, g.ny)

    dCdx = (C[ip, jm] + 2.0 * C[ip, j] + C[ip, jp] -
            C[im, jm] - 2.0 * C[im, j] - C[im, jp]) / (8.0 * g.dx)

    dCdy = (C[im, jp] + 2.0 * C[i, jp] + C[ip, jp] -
            C[im, jm] - 2.0 * C[i, jm] - C[ip, jm]) / (8.0 * g.dy)

    return dCdx, dCdy
end

function _youngs_gradient(vg::SplitVOFGrid3D, i::Int, j::Int, k::Int)
    g = vg.grid
    C = vg.fract

    im = _wrap(i - 1, g.nx)
    ip = _wrap(i + 1, g.nx)
    jm = _wrap(j - 1, g.ny)
    jp = _wrap(j + 1, g.ny)
    km = _wrap(k - 1, g.nz)
    kp = _wrap(k + 1, g.nz)

    sumx = 0.0
    sumy = 0.0
    sumz = 0.0

    for dj in -1:1
        wj = dj == 0 ? 2.0 : 1.0
        jj = _wrap(j + dj, g.ny)
        for dk in -1:1
            wk = dk == 0 ? 2.0 : 1.0
            kk = _wrap(k + dk, g.nz)
            w = wj * wk
            sumx += w * (C[ip, jj, kk] - C[im, jj, kk])
        end
    end

    for di in -1:1
        wi = di == 0 ? 2.0 : 1.0
        ii = _wrap(i + di, g.nx)
        for dk in -1:1
            wk = dk == 0 ? 2.0 : 1.0
            kk = _wrap(k + dk, g.nz)
            w = wi * wk
            sumy += w * (C[ii, jp, kk] - C[ii, jm, kk])
        end
    end

    for di in -1:1
        wi = di == 0 ? 2.0 : 1.0
        ii = _wrap(i + di, g.nx)
        for dj in -1:1
            wj = dj == 0 ? 2.0 : 1.0
            jj = _wrap(j + dj, g.ny)
            w = wi * wj
            sumz += w * (C[ii, jj, kp] - C[ii, jj, km])
        end
    end

    dCdx = sumx / (32.0 * g.dx)
    dCdy = sumy / (32.0 * g.dy)
    dCdz = sumz / (32.0 * g.dz)
    return dCdx, dCdy, dCdz
end

"""
Compute Youngs interface normals in 2D and store in `nxint`, `nyint`.
"""
function youngs_normals!(vg::SplitVOFGrid)
    g = vg.grid
    C = vg.fract

    for j in 1:g.ny, i in 1:g.nx
        if is_mixed(vg, i, j)
            gx, gy = _youngs_gradient(vg, i, j)
            m = hypot(gx, gy)

            if m <= 1.0e-14
                im = _wrap(i - 1, g.nx)
                ip = _wrap(i + 1, g.nx)
                jm = _wrap(j - 1, g.ny)
                jp = _wrap(j + 1, g.ny)
                gx = (C[ip, j] - C[im, j]) / (2.0 * g.dx)
                gy = (C[i, jp] - C[i, jm]) / (2.0 * g.dy)
                m = hypot(gx, gy)
            end

            if m <= 1.0e-14
                vg.nxint[i, j] = 1.0
                vg.nyint[i, j] = 0.0
            else
                vg.nxint[i, j] = gx / m
                vg.nyint[i, j] = gy / m
            end
        else
            vg.nxint[i, j] = 1.0
            vg.nyint[i, j] = 0.0
        end
    end

    return vg
end

"""
Compute Youngs interface normals in 3D and store in `nxint`, `nyint`, `nzint`.
"""
function youngs_normals!(vg::SplitVOFGrid3D)
    g = vg.grid
    C = vg.fract

    for k in 1:g.nz, j in 1:g.ny, i in 1:g.nx
        if is_mixed(vg, i, j, k)
            gx, gy, gz = _youngs_gradient(vg, i, j, k)
            m = sqrt(gx^2 + gy^2 + gz^2)

            if m <= 1.0e-14
                im = _wrap(i - 1, g.nx)
                ip = _wrap(i + 1, g.nx)
                jm = _wrap(j - 1, g.ny)
                jp = _wrap(j + 1, g.ny)
                km = _wrap(k - 1, g.nz)
                kp = _wrap(k + 1, g.nz)
                gx = (C[ip, j, k] - C[im, j, k]) / (2.0 * g.dx)
                gy = (C[i, jp, k] - C[i, jm, k]) / (2.0 * g.dy)
                gz = (C[i, j, kp] - C[i, j, km]) / (2.0 * g.dz)
                m = sqrt(gx^2 + gy^2 + gz^2)
            end

            if m <= 1.0e-14
                vg.nxint[i, j, k] = 1.0
                vg.nyint[i, j, k] = 0.0
                vg.nzint[i, j, k] = 0.0
            else
                vg.nxint[i, j, k] = gx / m
                vg.nyint[i, j, k] = gy / m
                vg.nzint[i, j, k] = gz / m
            end
        else
            vg.nxint[i, j, k] = 1.0
            vg.nyint[i, j, k] = 0.0
            vg.nzint[i, j, k] = 0.0
        end
    end

    return vg
end

"""
Compute height-function interface normals in 2D and store in `nxint`, `nyint`.
"""
function height_function_normals!(vg::SplitVOFGrid)
    g = vg.grid
    C = vg.fract

    @inbounds for j in 1:g.ny, i in 1:g.nx
        if is_mixed(vg, i, j)
            im = _wrap(i - 1, g.nx)
            ip = _wrap(i + 1, g.nx)
            hm = sum(C[im, _wrap(j + s, g.ny)] for s in -3:3)
            hp = sum(C[ip, _wrap(j + s, g.ny)] for s in -3:3)
            sx = (hp - hm) * g.dy / (2.0 * g.dx)

            jb = _wrap(j - 1, g.ny)
            jt = _wrap(j + 1, g.ny)
            vm = sum(C[_wrap(i + s, g.nx), jb] for s in -3:3)
            vp = sum(C[_wrap(i + s, g.nx), jt] for s in -3:3)
            sy = (vp - vm) * g.dx / (2.0 * g.dy)

            gx, gy = _youngs_gradient(vg, i, j)
            if abs(gx) >= abs(gy)
                nx = -sx
                ny = 1.0
            else
                nx = 1.0
                ny = -sy
            end

            m = hypot(nx, ny)
            if m <= 1.0e-14
                youngs_normals!(vg)
                return vg
            end
            nx /= m
            ny /= m

            # Keep the same orientation convention as Youngs normals.
            if nx * gx + ny * gy < 0.0
                nx = -nx
                ny = -ny
            end
            vg.nxint[i, j] = nx
            vg.nyint[i, j] = ny
        else
            vg.nxint[i, j] = 1.0
            vg.nyint[i, j] = 0.0
        end
    end

    return vg
end

function height_function_normals!(::SplitVOFGrid3D)
    throw(ArgumentError("HeightFunctionNormal is currently implemented only for 2D SplitVOF grids"))
end

"""
Compute ELVIRA-style interface normals in 2D and store in `nxint`, `nyint`.
"""
function elvira_normals!(vg::SplitVOFGrid)
    g = vg.grid
    verts = _reconstruct2d_verts()
    poly = _rect_polygon(0.0, 1.0, 0.0, 1.0)
    candidates = SVector{2, Float64}[]

    @inbounds for j in 1:g.ny, i in 1:g.nx
        if is_mixed(vg, i, j)
            empty!(candidates)
            _elvira_candidates!(candidates, vg, i, j)

            best_n = SVector(1.0, 0.0)
            best_score = Inf
            for n in candidates
                score = _elvira_score(vg, i, j, n, verts, poly)
                if score < best_score
                    best_score = score
                    best_n = n
                end
            end

            vg.nxint[i, j] = best_n[1]
            vg.nyint[i, j] = best_n[2]
        else
            vg.nxint[i, j] = 1.0
            vg.nyint[i, j] = 0.0
        end
    end

    return vg
end

function elvira_normals!(::SplitVOFGrid3D)
    throw(ArgumentError("ELVIRANormal is currently implemented only for 2D SplitVOF grids"))
end

function _elvira_candidates!(candidates, vg::SplitVOFGrid, i::Int, j::Int)
    g = vg.grid

    for dj in -1:1, di in -1:1
        ii = _wrap(i + di, g.nx)
        jj = _wrap(j + dj, g.ny)
        gx, gy = _youngs_gradient(vg, ii, jj)
        _push_normal_candidate!(candidates, gx, gy)

        im = _wrap(ii - 1, g.nx)
        ip = _wrap(ii + 1, g.nx)
        jm = _wrap(jj - 1, g.ny)
        jp = _wrap(jj + 1, g.ny)
        gx = (vg.fract[ip, jj] - vg.fract[im, jj]) / (2.0 * g.dx)
        gy = (vg.fract[ii, jp] - vg.fract[ii, jm]) / (2.0 * g.dy)
        _push_normal_candidate!(candidates, gx, gy)
    end

    for (nx, ny) in (
            (1.0, 0.0),
            (0.0, 1.0),
            (1.0, 1.0),
            (1.0, -1.0),
        )
        _push_normal_candidate!(candidates, nx, ny)
    end

    isempty(candidates) && push!(candidates, SVector(1.0, 0.0))
    return candidates
end

function _push_normal_candidate!(candidates, nx::Real, ny::Real)
    m = hypot(Float64(nx), Float64(ny))
    m <= 1.0e-14 && return candidates

    n = SVector(Float64(nx) / m, Float64(ny) / m)
    _push_unique_normal!(candidates, n)
    _push_unique_normal!(candidates, -n)
    return candidates
end

function _push_unique_normal!(candidates, n)
    any(c -> abs(dot(c, n)) > 1.0 - 1.0e-10 && dot(c, n) > 0.0, candidates) ||
        push!(candidates, n)
    return candidates
end

function _elvira_score(vg::SplitVOFGrid, i::Int, j::Int, n, verts, poly)
    g = vg.grid
    c = _plic_constant_2d(vg, i, j, n[1], n[2], verts)
    isfinite(c) || return Inf

    score = 0.0
    @inbounds for dj in -1:1, di in -1:1
        ii = _wrap(i + di, g.nx)
        jj = _wrap(j + dj, g.ny)
        pred = _plane_fraction_2d(vg, ii, jj, c, n[1], n[2], poly)
        err = pred - clamp(vg.fract[ii, jj], 0.0, 1.0)
        score += err * err
    end
    return score
end

function _plic_constant_2d(vg::SplitVOFGrid, i::Int, j::Int, nx::Real, ny::Real, verts)
    g = vg.grid
    xL, xR, yB, yT = cell_bounds(g, i, j)
    verts[1, 1] = xL; verts[1, 2] = yB
    verts[2, 1] = xR; verts[2, 2] = yB
    verts[3, 1] = xR; verts[3, 2] = yT
    verts[4, 1] = xL; verts[4, 2] = yT
    vliq = clamp(vg.fract[i, j], 0.0, 1.0) * g.dx * g.dy
    return enforv2dsz(g.dx, g.dy, vliq, verts, Float64(nx), Float64(ny))
end

function _plane_fraction_2d(vg::SplitVOFGrid, i::Int, j::Int, c, nx, ny, poly)
    g = vg.grid
    xL, xR, yB, yT = cell_bounds(g, i, j)
    _rect_polygon!(poly, xL, xR, yB, yT)
    _, icontp = inte2d!(poly, Float64(c), Float64(nx), Float64(ny))
    area = icontp == 0 ? 0.0 : toolv2d(poly)
    return clamp(area / (g.dx * g.dy), 0.0, 1.0)
end

"""
Reconstruct 2D PLIC planes (`cplic`) from the current volume fractions.
"""
function reconstruct!(vg::SplitVOFGrid)
    return reconstruct!(vg, YoungsPLIC())
end

function reconstruct!(vg::SplitVOFGrid, ::YoungsPLIC)
    return reconstruct!(vg, Reconstruction(YoungsNormal(), PLIC()))
end

function reconstruct!(vg::SplitVOFGrid, ::SLIC)
    return reconstruct!(vg, Reconstruction(YoungsNormal(), SLIC()))
end

function reconstruct!(vg::SplitVOFGrid, reconstruction::Reconstruction)
    _compute_normals!(vg, reconstruction.normal)
    _apply_interface!(vg, reconstruction.interface)
    _enforce_plic_volume!(vg)
    return vg
end

function _enforce_plic_volume!(vg::SplitVOFGrid)
    g = vg.grid
    cell_area = g.dx * g.dy
    verts = _reconstruct2d_verts()

    for j in 1:g.ny, i in 1:g.nx
        α = clamp(vg.fract[i, j], 0.0, 1.0)

        if is_mixed(vg, i, j)
            xL, xR, yB, yT = cell_bounds(g, i, j)
            verts[1, 1] = xL; verts[1, 2] = yB
            verts[2, 1] = xR; verts[2, 2] = yB
            verts[3, 1] = xR; verts[3, 2] = yT
            verts[4, 1] = xL; verts[4, 2] = yT

            vliq = α * cell_area
            vg.cplic[i, j] = enforv2dsz(g.dx, g.dy, vliq,
                                        verts,
                                        vg.nxint[i, j],
                                        vg.nyint[i, j])
        elseif is_full(vg, i, j)
            vg.cplic[i, j] = 1.0e6
        else
            vg.cplic[i, j] = -1.0e6
        end
    end

    return vg
end

"""
Reconstruct 3D PLIC planes (`cplic`) from the current volume fractions.
"""
function reconstruct!(vg::SplitVOFGrid3D)
    return reconstruct!(vg, YoungsPLIC())
end

function reconstruct!(vg::SplitVOFGrid3D, ::YoungsPLIC)
    return reconstruct!(vg, Reconstruction(YoungsNormal(), PLIC()))
end

function reconstruct!(vg::SplitVOFGrid3D, ::SLIC)
    return reconstruct!(vg, Reconstruction(YoungsNormal(), SLIC()))
end

function reconstruct!(vg::SplitVOFGrid3D, reconstruction::Reconstruction)
    _compute_normals!(vg, reconstruction.normal)
    _apply_interface!(vg, reconstruction.interface)
    _enforce_plic_volume!(vg)
    return vg
end

_compute_normals!(vg::AbstractSplitVOFGrid, ::YoungsNormal) = youngs_normals!(vg)
_compute_normals!(vg::AbstractSplitVOFGrid, ::HeightFunctionNormal) = height_function_normals!(vg)
_compute_normals!(vg::AbstractSplitVOFGrid, ::ELVIRANormal) = elvira_normals!(vg)

_apply_interface!(vg::AbstractSplitVOFGrid, ::PLIC) = vg
_apply_interface!(vg::AbstractSplitVOFGrid, ::SLIC) = _quantize_slic_normals!(vg)

function _enforce_plic_volume!(vg::SplitVOFGrid3D)
    g = vg.grid
    cell_vol = g.dx * g.dy * g.dz
    verts = _reconstruct3d_verts()

    for k in 1:g.nz, j in 1:g.ny, i in 1:g.nx
        α = clamp(vg.fract[i, j, k], 0.0, 1.0)

        if is_mixed(vg, i, j, k)
            xL, xR, yB, yT, zD, zU = cell_bounds(g, i, j, k)
            verts[1, 1] = xR; verts[1, 2] = yB; verts[1, 3] = zU
            verts[2, 1] = xR; verts[2, 2] = yB; verts[2, 3] = zD
            verts[3, 1] = xR; verts[3, 2] = yT; verts[3, 3] = zD
            verts[4, 1] = xR; verts[4, 2] = yT; verts[4, 3] = zU
            verts[5, 1] = xL; verts[5, 2] = yB; verts[5, 3] = zU
            verts[6, 1] = xL; verts[6, 2] = yB; verts[6, 3] = zD
            verts[7, 1] = xL; verts[7, 2] = yT; verts[7, 3] = zD
            verts[8, 1] = xL; verts[8, 2] = yT; verts[8, 3] = zU

            vliq = α * cell_vol
            vg.cplic[i, j, k] = enforv3dsz(g.dx, g.dy, g.dz, vliq,
                                           verts,
                                           vg.nxint[i, j, k],
                                           vg.nyint[i, j, k],
                                           vg.nzint[i, j, k])
        elseif is_full(vg, i, j, k)
            vg.cplic[i, j, k] = 1.0e6
        else
            vg.cplic[i, j, k] = -1.0e6
        end
    end

    return vg
end

function _quantize_slic_normals!(vg::SplitVOFGrid)
    g = vg.grid
    @inbounds for j in 1:g.ny, i in 1:g.nx
        if abs(vg.nxint[i, j]) >= abs(vg.nyint[i, j])
            vg.nxint[i, j] = signbit(vg.nxint[i, j]) ? -1.0 : 1.0
            vg.nyint[i, j] = 0.0
        else
            vg.nxint[i, j] = 0.0
            vg.nyint[i, j] = signbit(vg.nyint[i, j]) ? -1.0 : 1.0
        end
    end
    return vg
end

function _quantize_slic_normals!(vg::SplitVOFGrid3D)
    g = vg.grid
    @inbounds for k in 1:g.nz, j in 1:g.ny, i in 1:g.nx
        ax = abs(vg.nxint[i, j, k])
        ay = abs(vg.nyint[i, j, k])
        az = abs(vg.nzint[i, j, k])
        if ax >= ay && ax >= az
            vg.nxint[i, j, k] = signbit(vg.nxint[i, j, k]) ? -1.0 : 1.0
            vg.nyint[i, j, k] = 0.0
            vg.nzint[i, j, k] = 0.0
        elseif ay >= az
            vg.nxint[i, j, k] = 0.0
            vg.nyint[i, j, k] = signbit(vg.nyint[i, j, k]) ? -1.0 : 1.0
            vg.nzint[i, j, k] = 0.0
        else
            vg.nxint[i, j, k] = 0.0
            vg.nyint[i, j, k] = 0.0
            vg.nzint[i, j, k] = signbit(vg.nzint[i, j, k]) ? -1.0 : 1.0
        end
    end
    return vg
end

"""
Dispatch reconstruction based on `params.reconstruction`.
"""
function reconstruct!(vg::AbstractSplitVOFGrid, params::SplitVOFParams)
    return reconstruct!(vg, params.reconstruction)
end
