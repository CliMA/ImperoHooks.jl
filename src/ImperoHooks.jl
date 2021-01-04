module ImperoHooks

using Impero
using ClimateMachine

using ClimateMachine.Mesh.Grids
using ClimateMachine.Mesh.Elements
import ClimateMachine.Mesh.Elements: baryweights
using GaussQuadrature
using Base.Threads
using Makie, GLMakie, AbstractPlotting
using ImageTransformations, Colors
using AbstractPlotting.MakieLayout, LinearAlgebra

include("kernels.jl")
include("grid.jl")
include("calculus.jl")
include("utils.jl")
include("permutations.jl")
include("find_element.jl")
include("lagrange_interpolation.jl")
include("fields.jl")
include("gridhelper.jl")

end # module
