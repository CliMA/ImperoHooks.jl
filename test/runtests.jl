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


include(pwd() * "/test/gradient_test.jl")
include(pwd() * "/test/divergence_test.jl")
include(pwd() * "/test/gradient_divergence_test.jl")

include(pwd() * "/test/impero_algebra_test.jl")
include(pwd() * "/test/impero_curl_test.jl")

include(pwd() * "/test/gridhelper_test.jl")
include(pwd() * "/test/permutation_test.jl")