export GridHelper, InterpolationHelper, ElementHelper

struct InterpolationHelper{S, T}
    points::S
    quadrature::S
    interpolation::S 
    cartesianindex::T
end

function InterpolationHelper(g::DiscontinuousSpectralElementGrid)
    porders = polynomialorders(g)
    if length(porders) == 3
        npx, npy, npz = porders
        rx, wx = GaussQuadrature.legendre(npx+1, both)
        ωx = baryweights(rx)
        ry, wy = GaussQuadrature.legendre(npy+1, both)
        ωy = baryweights(ry)
        rz, wz = GaussQuadrature.legendre(npz+1, both)
        ωz = baryweights(rz)
        linlocal = reshape(collect(1:(npx+1)*(npy+1)*(npz+1)), (npx+1, npy+1, npz+1))
        return InterpolationHelper((rx, ry, rz), (wx, wy, wz), (ωx, ωy, ωz), linlocal)
    elseif length(porders) == 2
        npx, npy = porders
        rx, wx = GaussQuadrature.legendre(npx+1, both)
        ωx = baryweights(rx)
        ry, wy = GaussQuadrature.legendre(npy+1, both)
        ωy = baryweights(ry)
        linlocal = reshape(collect(1:(npx+1)*(npy+1)), (npx+1, npy+1))
        return InterpolationHelper((rx, ry), (wx, wy), (ωx, ωy), linlocal)
    else
        println("Not supported")  
        return nothing      
    end
    return nothing
end

struct ElementHelper{S, T, U, V, W}
    cellcenters::S
    coordinates::T
    cartesiansizes::U
    permutation::V
    cartesianindex::W
end

addup(xC, tol) = sum(abs.(xC[1] .- xC) .≤ tol)

function ElementHelper(g::DiscontinuousSpectralElementGrid)
    porders = polynomialorders(g)
    x, y, z = coordinates(g)
    xC, yC, zC = cellcenters(g)
    ne = size(x)[2]
    ex = round(Int64, ne / addup(xC, 10^4 * eps(maximum(abs.(x)))))
    ey = round(Int64, ne / addup(yC, 10^4 * eps(maximum(abs.(y)))))
    ez = round(Int64, ne / addup(zC, 10^4 * eps(maximum(abs.(z)))))

    check = ne == ex * ey * ez
    check ? true : error("improper counting")
    p = getperm(xC, yC, zC, ex, ey, ez)

    if length(porders) == 3
        npx, npy, npz = porders
        lin = reshape(collect(1:length(xC)), (ex, ey, ez))
        return ElementHelper((xC, yC, zC), (x, y, z), (ex, ey, ez), p, lin)
    elseif length(porders) == 2
        npx, npy = porders
        lin = reshape(collect(1:length(xC)), (ex, ey))
        check = ne == ex * ey 
        check ? true : error("improper counting")
        return ElementHelper((xC, yC), (x, y), (ex, ey), p, lin)
    else
        println("no constructor for polynomial order = ", porders)
        return nothing
    end
    return nothing
end

struct GridHelper{S,T,V}
    interpolation::S
    element::T
    grid::V
end

function GridHelper(g::DiscontinuousSpectralElementGrid)
    return GridHelper(InterpolationHelper(g), ElementHelper(g), g)
end

function getvalue(f, location, gridhelper::GridHelper)
    ih = gridhelper.interpolation
    eh = gridhelper.element
    porders = polynomialorders(gridhelper.grid)
    if length(porders) == 3
        npx, npy, npz = polynomialorders(gridhelper.grid)
        fl = reshape(f, (npx+1, npy+1, npz+1, prod(eh.cartesiansizes)))
        ip = getvalue(fl, eh.cellcenters..., location, 
                eh.permutation, eh.cartesianindex, ih.cartesianindex, 
                eh.coordinates..., ih.points..., ih.interpolation...)
        return ip
    elseif length(porders) == 2
        npx, npy = polynomialorders(gridhelper.grid)
        fl = reshape(f, (npx+1, npy+1, prod(eh.cartesiansizes)))
        ip = getvalue(fl, eh.cellcenters..., location, 
                eh.permutation, eh.cartesianindex, ih.cartesianindex, 
                eh.coordinates..., ih.points..., ih.interpolation...)
        return ip
    end
    return nothing
end