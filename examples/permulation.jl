using ImperoHooks
using ClimateMachine
using ClimateMachine.Mesh.Grids
using Impero, MPI
using Makie, GLMakie, AbstractPlotting
using ImageTransformations, Colors
using AbstractPlotting.MakieLayout
using GaussQuadrature, Random, Test

ClimateMachine.init()
const ArrayType = ClimateMachine.array_type()
const mpicomm = MPI.COMM_WORLD
const FT = Float64
Ω = Circle(-1,1) × Circle(-1,1) × Circle(-1,1)
dims = ndims(Ω)
ex, ey, ez = (2,3,4)
ClimateMachine.gpu_allowscalar(true)
# Define Grid: might want to loop over element sizes and polynomial orders
grid = DiscontinuousSpectralElementGrid(Ω, elements = (ex, ey, ez), polynomialorder = (4,4,4), array = ArrayType)

x, y, z = coordinates(grid)
##
xC, yC, zC = cellcenters(grid)
function getperm(xC, yC, zC, ex, ey, ez)
    pz = sortperm(zC)
    tmpZ = reshape(zC[pz], (ex, ey, ez))
    tmpY = reshape(yC[pz], (ex*ey, ez))
    tmp_py = [sortperm(tmpY[:,i]) for i in 1:ez]
    py = zeros(Int64,length(pz))
    for i in eachindex(tmp_py)
        n = length(tmp_py[i]) 
        ii = (i-1) * n + 1
        py[ii : ii + n-1] .= tmp_py[i] .+ ii .- 1
    end
    newY = reshape(yC[pz][py], (ex, ey, ez))
    newZ = reshape(zC[pz][py], (ex, ey, ez))
    tmp_px = [sortperm(reshape(xC[pz][py], (ex, ey*ez))[:,i]) for i in 1:ey*ez]
    px = zeros(Int64,length(pz))
    for i in eachindex(tmp_px)
        n = length(tmp_px[i]) 
        ii = (i-1) * n + 1
        px[ii : ii + n-1] .= tmp_px[i] .+ ii .- 1
    end
    p = [pz[py[px[i]]] for i in eachindex(px)]
    return p
end
p = getperm(xC, yC, zC, ex, ey, ez)
newX = reshape(xC[p], (ex, ey, ez))
newY = reshape(yC[p], (ex, ey, ez))
newZ = reshape(zC[p], (ex, ey, ez))
lin = reshape(collect(1:length(xC)), (ex, ey, ez))

@testset "Unique Permutation" begin
    @test prod( sort(p) .== collect(1:length(p)))
end
@testset "Check lin" begin
    for i in 1:10
        rng = MersenneTwister(i);
        ind = lin[rand(rng, 1:ex), rand(rng, 1:ey), rand(rng, 1:ez)]
        @test newX[ind] == xC[p[ind]]
        @test newY[ind] == yC[p[ind]]
        @test newZ[ind] == zC[p[ind]]
    end
end

@testset "Check Permutation" begin
    tol = eps(2.0)
    for i in 1:ex
        @test prod( abs.(newX[i, :, :] .- newX[i, 1, 1]) .< tol)
    end
    for i in 1:ey
        @test prod( abs.(newY[:, i, :] .- newY[1, i, 1]) .< tol)
    end
    for i in 1:ez
        @test prod( abs.(newZ[:, :, i] .- newZ[1, 1, i]) .< tol)
    end
end