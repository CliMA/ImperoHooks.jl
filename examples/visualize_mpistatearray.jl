using ImperoHooks
using ClimateMachine
using ClimateMachine.Mesh.Grids
using Impero, MPI
using Makie, GLMakie, AbstractPlotting
using ImageTransformations, Colors
using AbstractPlotting.MakieLayout
using GaussQuadrature, Random, Test
using Statistics, Printf
using ClimateMachine.MPIStateArrays
include(pwd() * "/examples/gridhelper.jl")
##
ClimateMachine.init()
const ArrayType = ClimateMachine.array_type()
const mpicomm = MPI.COMM_WORLD
const FT = Float64
# Domain
Ω = Circle(-1,1) × Circle(-1,1) × Circle(-1,1)
# Discretization
ex, ey, ez = (3,4,5) .* 10
(npx, npy, npz) = (4,4,4)
ClimateMachine.gpu_allowscalar(true)
grid = DiscontinuousSpectralElementGrid(Ω, elements = (ex, ey, ez), polynomialorder = (npx, npy, npz), array = ArrayType)
gridhelper = GridHelper(grid)
x, y, z = coordinates(grid)
xC, yC, zC = cellcenters(grid)
##
f = MPIStateArray{FT}(mpicomm, ArrayType, size(x)..., 1)
f .=  1 .- (x .^2 + y .^2 + z .^2)

newsize = (10, 10, 10)
nx, ny, nz = newsize
newx = range(-1,1, length = nx)
newy = range(-1,1, length = ny)
newz = range(-1,1, length = nz)
newf = zeros((nx,ny,nz))
for i in 1:nx
    for j in 1:ny
        for k in 1:nz
            location = (newx[i], newy[j], newz[k])
            newf[i,j,k] = getvalue(f, location, gridhelper)
        end
    end
end
##
scene, layout = layoutscene(resolution = (1920, 860))
lscene = layout[2:4,2:4] = LScene(scene)
layout[1, 2:4] = LText(scene, "MPI State Array", textsize = 50)

upperclim_slider = LSlider(scene, range = range(0, 1, length = 101), startvalue = 0.99)
upperclim_node = upperclim_slider.value
lowerclim_slider = LSlider(scene, range = range(0, 1, length = 101), startvalue = 0.01)
lowerclim_node = lowerclim_slider.value

clims = @lift((quantile(f[:], $lowerclim_node) , quantile(f[:], $upperclim_node)))
cmap_rgb = to_colormap(:balance);
volume!(lscene, -1..1, -1..1, -1..1, newf, colormap = cmap_rgb, colorrange = clims)

    lowerclim_string = @lift("lower clim quantile = " *  @sprintf("%0.2f", $lowerclim_node) * ", value = " * @sprintf("%0.1e", $clims[1]))
    upperclim_string = @lift("upper clim quantile = " *  @sprintf("%0.2f", $upperclim_node) * ", value = " * @sprintf("%0.1e", $clims[2]))
layout[1, 1] = vgrid!(
    LText(scene, lowerclim_string, width = nothing),
    lowerclim_slider,
    LText(scene, upperclim_string, width = nothing),
    upperclim_slider,
)
display(scene)

