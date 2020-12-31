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
ex, ey, ez = (3,4,5) 
(npx, npy, npz) = (4,4,4)
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

function checkgl_2(x, rx)
    for i in eachindex(rx)
        if abs(x - rx[i]) ≤ eps(rx[i])  return i end
    end
    return 0
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
    icheck = checkgl_2(newx, rx)
    jcheck = checkgl_2(newy, ry)
    kcheck = checkgl_2(newz, rz)
    numerator = zeros(1)
    denominator = zeros(1)
    for k in eachindex(rz)
        if kcheck ==0
            Δz = (newz .- rz[k])
            polez = ωz[k] ./ Δz
            kk = k
        else
            polez = 1.0
            k = eachindex(rz)[end]
            kk = kcheck
        end
        for j in eachindex(ry)
            if jcheck ==0
                Δy = (newy .- ry[j])
                poley = ωy[j] ./ Δy
                jj = j
            else
                poley = 1.0
                j = eachindex(ry)[end]
                jj = jcheck
            end
            for i in eachindex(rx)
                if icheck == 0
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

function lagrange_eval_nocheck(f, newx, newy, newz, rx, ry, rz, ωx, ωy, ωz)
    numerator = zeros(1)
    denominator = zeros(1)
    for k in eachindex(rz)
            Δz = (newz .- rz[k])
            polez = ωz[k] ./ Δz
            kk = k
        for j in eachindex(ry)
                Δy = (newy .- ry[j])
                poley = ωy[j] ./ Δy
                jj = j
            for i in eachindex(rx)
                    Δx = (newx .- rx[i])
                    polex = ωx[i] ./ Δx
                    ii = i
                    numerator[1] += f[ii,jj,kk] * polex * poley * polez
                    denominator[1] += polex * poley * polez
            end
        end
    end
    return numerator[1] / denominator[1]
end

function lagrange_eval_nocheck(f, newx, newy, rx, ry, ωx, ωy)
    numerator = zeros(1)
    denominator = zeros(1)
    for j in eachindex(ry)
            Δy = (newy .- ry[j])
            poley = ωy[j] ./ Δy
        for i in eachindex(rx)
                Δx = (newx .- rx[i])
                polex = ωx[i] ./ Δx
                numerator[1] += f[i,j] * polex * poley
                denominator[1] += polex * poley
        end
    end
    return numerator[1] / denominator[1]
end

function lagrange_eval_nocheck(f, newx, rx, ωx)
    numerator = zeros(1)
    denominator = zeros(1)
    for i in eachindex(rx)
            Δx = (newx .- rx[i])
            polex = ωx[i] ./ Δx
            numerator[1] += f[i] * polex * poley
            denominator[1] += polex * poley
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

#=
f = copy(x)
fl = reshape(f, (npx+1, npy+1, npz+1, ex*ey*ez))
newx, newy, newz = (0.3, 0.1, 0.1)
lagrange_eval(fl, newx, newy, newz, rx, ry, rz, ωx, ωy, ωz)
=#

#=
@code_warntype checkgl(newx, rx)
@code_warntype checkgl_2(newx, rx)
@btime icheck = checkgl(newx, rx)
@btime lagrange_eval(f, newx, newy, newz, rx, ry, rz, ωx, ωy, ωz)
@btime lagrange_eval_nocheck(f, newx, newy, newz, rx, ry, rz, ωx, ωy, ωz)
=#

## 
xC, yC, zC = cellcenters(grid)
p = getperm(xC, yC, zC, ex, ey, ez)
lin = reshape(collect(1:length(xC)), (ex, ey, ez))
linlocal = reshape(collect(1:(npx+1)*(npy+1)*(npz+1)), (npx+1, npy+1, npz+1))

f = copy(x)
fl = reshape(f, (npx+1, npy+1, npz+1, ex*ey*ez))

location = (0.23, -0.3, 0.5)

# computation
function getvalue(fl, xC, yC, zC, location, p, lin, linlocal, x, y, z, rx, ry, rz, ωx, ωy, ωz)
    e = findelement(xC, yC, zC, location, p, lin)
    xmax = x[linlocal[npx+1,1,1], e]
    xmin = x[linlocal[1,1,1], e]
    ymax = y[linlocal[1,npy+1,1], e]
    ymin = y[linlocal[1,1,1], e]
    zmax = z[linlocal[1,1,npz+1], e]
    zmin = z[linlocal[1,1,1], e]

    newx = 2 * (location[1] - xmin) / (xmax - xmin) - 1
    newy = 2 * (location[2] - ymin) / (ymax - ymin) - 1
    newz = 2 * (location[3] - zmin) / (zmax - zmin) - 1

    return lagrange_eval(view(fl,:,:,:,e), newx, newy, newz, rx, ry, rz, ωx, ωy, ωz)
end

@btime getvalue(fl, xC, yC, zC, location, p, lin, linlocal, x, y, z, rx, ry, rz, ωx, ωy, ωz)