using ImperoHooks
using ClimateMachine
using ClimateMachine.Mesh.Grids
using Impero, MPI
using Makie, GLMakie, AbstractPlotting
using ImageTransformations, Colors
using AbstractPlotting.MakieLayout
using GaussQuadrature, Random, Test
include(pwd() * "/examples/permutation.jl")
##
ClimateMachine.init()
const ArrayType = ClimateMachine.array_type()
const mpicomm = MPI.COMM_WORLD
const FT = Float64
Ω = Circle(-1,1) × Circle(-1,1) × Circle(-1,1)
dims = ndims(Ω)
ex, ey, ez = (1,1,1) 
(npx, npy, npz) = (2,2,2)
ClimateMachine.gpu_allowscalar(true)
# Define Grid: might want to loop over element sizes and polynomial orders
grid = DiscontinuousSpectralElementGrid(Ω, elements = (ex, ey, ez), polynomialorder = (npx, npy, npz), array = ArrayType)

x, y, z = coordinates(grid)
xC, yC, zC = cellcenters(grid)

# Functions
function checkgl(x, rx)
    for i in eachindex(rx)
        if abs(x - rx[i]) ≤ eps(rx[i])  return i end
    end
    return true
end

function loopchecker(icheck, rx, newx)
    if kcheck
        Δz = (newz .- rz[k])
        polez = ωz[k] ./ Δz
        kk = k
    else
        Δz = 0.0
        polez = 1.0
        k = npz + 1
        kk = kcheck
    end
    return 
end

function lagrange_eval(f, newx, newy, newz, rx, ry, rz, ωx, ωy, ωz)
    icheck = checkgl(newx, rx)
    jcheck = checkgl(newy, ry)
    kcheck = checkgl(newz, rz)
    numerator = zeros(1)
    denominator = zeros(1)
    for k in eachindex(rz)
        if kcheck
            Δz = (newz .- rz[k])
            polez = ωz[k] ./ Δz
            kk = k
        else
            polez = 1.0
            k = eachindex(rz)[end]
            kk = kcheck
        end
        for j in eachindex(ry)
            if jcheck
                Δy = (newy .- ry[j])
                poley = ωy[j] ./ Δy
                jj = j
            else
                poley = 1.0
                j = eachindex(ry)[end]
                jj = jcheck
            end
            for i in eachindex(rx)
                if icheck
                    Δx = (newx .- rx[i])
                    polex = ωx[i] ./ Δx
                    ii = i
                else
                    polex = 1.0
                    i = eachindex(rx)[end]
                    ii = icheck
                end
                numerator[1] += f[ii,jj,kk] * polex * poley * polez
                denominator[1] += polex * poley * polez
            end
        end
    end
    return numerator[1] / denominator[1]
end
##
# Check
# (npx, npy, npz) = polynomialorders(grid)
rx, wx = GaussQuadrature.legendre(npx+1, both)
ωx = baryweights(rx)
ry, wy = GaussQuadrature.legendre(npy+1, both)
ωy = baryweights(ry)
rz, wz = GaussQuadrature.legendre(npz+1, both)
ωz = baryweights(rz)

f = copy(reshape(x, (npx+1, npy+1, npz+1)))
newx, newy, newz = (0.3, 0.1, 0.1)
lagrange_eval(f, newx, newy, newz, rx, ry, rz, ωx, ωy, ωz)

@btime lagrange_eval(f, newx, newy, newz, rx, ry, rz, ωx, ωy, ωz)


##
gx = LagrangeInterpolant(copy(rx), rx, wx, ωx)
gy = LagrangeInterpolant(copy(ry), ry, wy, ωy)
gz = LagrangeInterpolant(copy(rz), rz, wz, ωz)

for k in 1:(npz+1)
    for j in 1:(npy+1)
        gx.values .= view(f, :, j, k)
        fjk = gx(newx)
    end
end

