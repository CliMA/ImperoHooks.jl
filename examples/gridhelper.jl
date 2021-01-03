using ImperoHooks
using ClimateMachine
using ClimateMachine.Mesh.Grids
using Impero, MPI
using Makie, GLMakie, AbstractPlotting
using ImageTransformations, Colors
using AbstractPlotting.MakieLayout
using GaussQuadrature, Random, Test
using ClimateMachine.Mesh.Elements
import ClimateMachine.Mesh.Elements: baryweights
import Impero: Circle # Impero circle is more important?
include(pwd() * "/examples/find_elements.jl")
include(pwd() * "/examples/lagrange_interpolant.jl")

struct InterpolationHelper{S, T}
    points::S
    quadrature::S
    interpolation::S 
    cartesianindex::T
end

function InterpolationHelper(g::DiscontinuousSpectralElementGrid)
    porders = polynomialorders(grid)
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
    porders = polynomialorders(grid)
    x, y, z = coordinates(grid)
    xC, yC, zC = cellcenters(grid)
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
    return GridHelper(InterpolationHelper(grid), ElementHelper(grid), grid)
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

## 
# front load costs
# quadrature points
Ω = Circle(-1,1) × Circle(-1,1) × Circle(-1,1)
dims = ndims(Ω)
ex, ey, ez = (3,4,5) .* 10
(npx, npy, npz) = (1,2,3) # needs to be at least these polynomial orders due to last test
ClimateMachine.gpu_allowscalar(true)
# Define Grid: might want to loop over element sizes and polynomial orders
grid = DiscontinuousSpectralElementGrid(Ω, elements = (ex, ey, ez), polynomialorder = (npx, npy, npz), array = ArrayType)
gridhelper = GridHelper(grid)

npx, npy, npz = polynomialorders(grid)
x, y, z = coordinates(grid)
xC, yC, zC = cellcenters(grid)

ne = size(x)[2]
ex = round(Int64, ne / addup(xC, 10^4 * eps(maximum(abs.(x)))))
ey = round(Int64, ne / addup(yC, 10^4 * eps(maximum(abs.(y)))))
ez = round(Int64, ne / addup(zC, 10^4 * eps(maximum(abs.(z)))))
check = ne == ex * ey * ez
check ? true : error("improper counting")
p = getperm(xC, yC, zC, ex, ey, ez)
lin = reshape(collect(1:length(xC)), (ex, ey, ez))


ne = size(x)[2]
rx, wx = GaussQuadrature.legendre(npx+1, both)
ωx = baryweights(rx)
ry, wy = GaussQuadrature.legendre(npy+1, both)
ωy = baryweights(ry)
rz, wz = GaussQuadrature.legendre(npz+1, both)
ωz = baryweights(rz)
linlocal = reshape(collect(1:(npx+1)*(npy+1)*(npz+1)), (npx+1, npy+1, npz+1))

# computation


# now test
@testset "Check Interpolation" begin
    f = copy(z)
    fl = reshape(f, (npx+1, npy+1, npz+1, ex*ey*ez))
    location = (0.0, 0.0, -0.3)
    @test location[3] ≈ getvalue(fl, xC, yC, zC, location, p, lin, linlocal, x, y, z, rx, ry, rz, ωx, ωy, ωz)
    @. f = x
    fl = reshape(f, (npx+1, npy+1, npz+1, ex*ey*ez))
    location = (0.1, 0.0, -0.3)
    @test location[1] ≈ getvalue(fl, xC, yC, zC, location, p, lin, linlocal, x, y, z, rx, ry, rz, ωx, ωy, ωz)
    @. f = y
    fl = reshape(f, (npx+1, npy+1, npz+1, ex*ey*ez))
    location = (0.1, -0.9, -0.3)
    @test location[2] ≈ getvalue(fl, xC, yC, zC, location, p, lin, linlocal, x, y, z, rx, ry, rz, ωx, ωy, ωz)
end
##
@testset "Check Interpolation" begin
    f = copy(z)
    location = (0.0, 0.0, -0.3)
    @test location[3] ≈ getvalue(f, location, gridhelper)
    @. f = x
    location = (0.1, 0.0, -0.3)
    @test location[1] ≈ getvalue(f, location, gridhelper)
    @. f = y
    location = (0.1, -0.9, -0.3)
    @test location[2] ≈ getvalue(f, location, gridhelper)
    @. f = x + y^2 + z^3
    location = (0.1, -0.9, -0.3)
    value = location[1] + location[2]^2 + location[3]^3
    @test value ≈ getvalue(f, location, gridhelper)
end


