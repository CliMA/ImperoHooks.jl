import Base: ndims, getindex
abstract type TensorField{S,T} end

struct ScalarField{S,T} <: TensorField{0,0}
    data::S
    grid::T
end

function (ϕ::ScalarField)(x::Tuple{Number, Number, Number})
    return getvalue(ϕ.data, x, ϕ.grid)
end

function (ϕ::ScalarField)(x::Number,y::Number,z::Number)
    return getvalue(ϕ.data, (x,y,z), ϕ.grid)
end

getindex(ϕ::ScalarField, i::Int) = ϕ.data[i]

ndims(ϕ::ScalarField) = ndims(ϕ.data)

function broadcast!(identity, ϕ::ScalarField, x::AbstractArray) 
    ϕ.data .= x 
    return nothing
end

function broadcast!(identity, ϕ::ScalarField, φ::ScalarField) 
    ϕ.data .= φ.data
    return nothing
end

## Test
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
using Base.Threads, LinearAlgebra
include(pwd() * "/examples/gridhelper.jl")
ex, ey, ez = (3,5,7)
(npx, npy, npz) = (1,1,1)
ClimateMachine.gpu_allowscalar(true)
grid = DiscontinuousSpectralElementGrid(Ω, elements = (ex, ey, ez), polynomialorder = (npx, npy, npz), array = ArrayType)
gridhelper = GridHelper(grid)
x, y, z = coordinates(grid)
ϕ =  ScalarField(copy(x), gridhelper)

ϕ(0.1, 0.2, 0.3)
ϕ.data .= y
ϕ(0.1, 0.2, 0.3)
