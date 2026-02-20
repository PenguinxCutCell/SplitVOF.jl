@inline function _wrap(i::Int, n::Int)
    if i < 1
        return i + n
    elseif i > n
        return i - n
    end
    return i
end

@inline function cell_bounds(g::CartesianGrid2D, i::Int, j::Int)
    xL = g.xlo + (i - 1) * g.dx
    xR = xL + g.dx
    yB = g.ylo + (j - 1) * g.dy
    yT = yB + g.dy
    return xL, xR, yB, yT
end

@inline function cell_bounds(g::CartesianGrid3D, i::Int, j::Int, k::Int)
    xL = g.xlo + (i - 1) * g.dx
    xR = xL + g.dx
    yB = g.ylo + (j - 1) * g.dy
    yT = yB + g.dy
    zD = g.zlo + (k - 1) * g.dz
    zU = zD + g.dz
    return xL, xR, yB, yT, zD, zU
end

@inline function cell_center(g::CartesianGrid2D, i::Int, j::Int)
    x = g.xlo + (i - 0.5) * g.dx
    y = g.ylo + (j - 0.5) * g.dy
    return SVector{2,Float64}(x, y)
end

@inline function cell_center(g::CartesianGrid3D, i::Int, j::Int, k::Int)
    x = g.xlo + (i - 0.5) * g.dx
    y = g.ylo + (j - 0.5) * g.dy
    z = g.zlo + (k - 0.5) * g.dz
    return SVector{3,Float64}(x, y, z)
end

@inline cell_center(vg::SplitVOFGrid, i::Int, j::Int) = cell_center(vg.grid, i, j)
@inline cell_center(vg::SplitVOFGrid3D, i::Int, j::Int, k::Int) = cell_center(vg.grid, i, j, k)
