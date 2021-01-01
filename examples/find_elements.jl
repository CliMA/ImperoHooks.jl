using ImperoHooks
using ClimateMachine
using ClimateMachine.Mesh.Grids
using Impero, MPI
using Makie, GLMakie, AbstractPlotting
using ImageTransformations, Colors
using AbstractPlotting.MakieLayout
using GaussQuadrature, Random, Test
include(pwd() * "/examples/permutation.jl")

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
xC, yC, zC = cellcenters(grid)

location = (0.0,0.0,0.0)
function distance(x,y,z,location)
    return (x .- location[1]) .^ 2 + (y .- location[2]) .^ 2 + (z .- location[3]) .^ 2
end

function findlocation(x, y, z, location)
    currentmin = ones(1)
    minind = ones(eltype(eachindex(x)), 1)
    currentmin[1] = (x[1] .- location[1]) .^ 2 + (y[1] .- location[2]) .^ 2 + (z[1] .- location[3]) .^ 2
    for i in eachindex(x)
        current = (x[i] .- location[1]) .^ 2 + (y[i] .- location[2]) .^ 2 + (z[i] .- location[3]) .^ 2
        if current < currentmin[1]
            currentmin[1] = current
            minind[1] = i
        end
    end
    return currentmin, minind
end

function findlocation(x, y, z, xC, yC, zC, location)
    e = findlocation(xC, yC, zC, location)[2]
    minval, ijk = findlocation(view(x,:,e), view(y,:,e), view(z,:,e), location)
    return [minval, (ijk, e)]
end

function findlocation(x, y, z, location)
    e = findlocation(xC, yC, zC, location)[2]
    minval, ijk = findlocation(view(x,:,e), view(y,:,e), view(z,:,e), location)
    return [minval, (ijk, e)]
end

function findelement(x,y,z, location)
    for e in 1:size(x)[2]
        xmin, xmax = extrema(x[:,e])
        if  xmin ≤ location[1] ≤ xmax
            ymin, ymax = extrema(y[:,e])
            if ymin ≤ location[2] ≤ ymax
                zmin, zmax = extrema(z[:,e])
                if zmin ≤ location[3] ≤ zmax
                    return e
                end
            end
        end
    end
    error("location was not found in rectangle")
    return nothing
end

function findelement(xC, yC, zC, location, p, lin)
    ex, ey, ez = size(lin)
    # i 
    currentmin = ones(1)
    minind = ones(Int64, 1)
    currentmin[1] = abs.(xC[p[lin[1,1,1]]] .- location[1])
    for i in 2:ex
        current = abs.(xC[p[lin[i,1,1]]] .- location[1])
        if current < currentmin[1]
            currentmin[1] = current
            minind[1] = i
        end
    end
    i = minind[1]
    # j 
    currentmin[1] = abs.(yC[p[lin[1,1,1]]] .- location[2])
    for i in 2:ey
        current = abs.(yC[p[lin[1,i,1]]] .- location[2])
        if current < currentmin[1]
            currentmin[1] = current
            minind[1] = i
        end
    end
    j = minind[1]
    # k 
    currentmin[1] = abs.(zC[p[lin[1,1,1]]] .- location[3])
    for i in 2:ez
        current = abs.(zC[p[lin[1,1, i]]] .- location[3])
        if current < currentmin[1]
            currentmin[1] = current
            minind[1] = i
        end
    end
    k = minind[1]
    return p[lin[i,j,k]]
end


xC, yC, zC = cellcenters(grid)
p = getperm(xC, yC, zC, ex, ey, ez)
lin = reshape(collect(1:length(xC)), (ex, ey, ez))
# best case = low index, worst case = large index
ind = length(xC)
location = (xC[ind], yC[ind], zC[ind])
## Check Times
#=
@btime argmin(distance(x, y, z, location))
@btime findlocation(x, y, z,location)
@btime findlocation(xC, yC, zC, location)
@btime cellcenters(grid)
@btime findlocation(x, y, z, xC, yC, zC, location)
@btime findelement(x, y, z, location)
@btime getperm(xC, yC, zC, ex, ey, ez)
@btime findelement(xC, yC, zC, location, p, lin)
=#
location = (xC[2], yC[2], zC[2])
e = findelement(xC, yC, zC, location, p, lin)
minimum(distance(x[:,e],y[:,e],z[:,e],location))
e = findelement(x, y, z, location)
minimum(distance(x[:,e],y[:,e],z[:,e],location))