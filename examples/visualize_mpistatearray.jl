using ImperoHooks
using ClimateMachine
using ClimateMachine.Mesh.Grids
using Impero, MPI
using Makie, GLMakie, AbstractPlotting
using ImageTransformations, Colors
using AbstractPlotting.MakieLayout
using GaussQuadrature, Random, Test
include(pwd() * "/examples/permulation.jl")

ClimateMachine.init()
const ArrayType = ClimateMachine.array_type()
const mpicomm = MPI.COMM_WORLD
const FT = Float64
Ω = Circle(-1,1) × Circle(-1,1) × Circle(-1,1)
dims = ndims(Ω)
ex, ey, ez = (2,3,4) 
(npx, npy, npz) = (2,2,2)
ClimateMachine.gpu_allowscalar(true)
# Define Grid: might want to loop over element sizes and polynomial orders
grid = DiscontinuousSpectralElementGrid(Ω, elements = (ex, ey, ez), polynomialorder = (npx, npy, npz), array = ArrayType)

x, y, z = coordinates(grid)
xC, yC, zC = cellcenters(grid)

# (npx, npy, npz) = polynomialorders(grid)
rx, wx = GaussQuadrature.legendre(npx+1, both)
ωx = baryweights(rx)
ry, wy = GaussQuadrature.legendre(npx+1, both)
ωy = baryweights(ry)
rz, wz = GaussQuadrature.legendre(npx+1, both)
ωz = baryweights(rz)

newx, newy, newz = (0,0,0)
ωx ./ (newx .- rx)
numerator = 0.0
denominator = 0.0

for k in 1:(npz+1)
    Δz = (newz .- rz[k])
    polez = ωz[j] ./ Δz
    for j in 1:(npy+1)
        Δy = (newy .- ry[j])
        polej = ωy[j] ./ Δy
        for i in 1:(npx+1)
            Δx = newx - rx[i]
            polex = ωx[i] ./ Δx
            pole = polex * poley * polez
            numerator += f[i,j,k] * pole
            denominator += pole
        end
    end
end
##
