using ImperoHooks
using ClimateMachine
using ClimateMachine.Mesh.Grids
using Impero, MPI
using Makie, GLMakie, AbstractPlotting
using ImageTransformations, Colors
using AbstractPlotting.MakieLayout
using GaussQuadrature, Random, Test
include(pwd() * "/examples/permutation.jl")
include(pwd() * "/examples/lagrange_interpolant.jl")
##
ClimateMachine.init()
const ArrayType = ClimateMachine.array_type()
const mpicomm = MPI.COMM_WORLD
const FT = Float64
Ω = Circle(-1,1) × Circle(-1,1) × Circle(-1,1)
dims = ndims(Ω)
ex, ey, ez = (3,4,5) .* 10
(npx, npy, npz) = (4,4,4)
ClimateMachine.gpu_allowscalar(true)
# Define Grid: might want to loop over element sizes and polynomial orders
grid = DiscontinuousSpectralElementGrid(Ω, elements = (ex, ey, ez), polynomialorder = (npx, npy, npz), array = ArrayType)

x, y, z = coordinates(grid)
xC, yC, zC = cellcenters(grid)

##
# Check
# (npx, npy, npz) = polynomialorders(grid)

#=
@code_warntype checkgl(newx, rx)
@code_warntype checkgl_2(newx, rx)
@btime icheck = checkgl(newx, rx)
@btime lagrange_eval(f, newx, newy, newz, rx, ry, rz, ωx, ωy, ωz)
@btime lagrange_eval_nocheck(f, newx, newy, newz, rx, ry, rz, ωx, ωy, ωz)
=#

## 
# front load costs
# quadrature points
npx, npy, npz = polynomialorders(grid)
x, y, z = coordinates(grid)
ne = size(x)[2]
rx, wx = GaussQuadrature.legendre(npx+1, both)
ωx = baryweights(rx)
ry, wy = GaussQuadrature.legendre(npy+1, both)
ωy = baryweights(ry)
rz, wz = GaussQuadrature.legendre(npz+1, both)
ωz = baryweights(rz)
# element getter
xC, yC, zC = cellcenters(grid)
ex = round(Int64, ne / sum(xC[1] .≈ xC))
ey = round(Int64, ne / sum(yC[1] .≈ yC))
ez = round(Int64, ne / sum(zC[1] .≈ zC))
check = ne == ex * ey * ez
check ? true : error("improper counting")
p = getperm(xC, yC, zC, ex, ey, ez)
lin = reshape(collect(1:length(xC)), (ex, ey, ez))
linlocal = reshape(collect(1:(npx+1)*(npy+1)*(npz+1)), (npx+1, npy+1, npz+1))


# computation


# now test
f = copy(z)
fl = reshape(f, (npx+1, npy+1, npz+1, ex*ey*ez))
location = (0.0, 0.0, -0.3)
getvalue(fl, xC, yC, zC, location, p, lin, linlocal, x, y, z, rx, ry, rz, ωx, ωy, ωz)
#=
@btime getvalue(fl, xC, yC, zC, location, p, lin, linlocal, x, y, z, rx, ry, rz, ωx, ωy, ωz)
=#