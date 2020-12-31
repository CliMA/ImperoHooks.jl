using GaussQuadrature, Plots
using ClimateMachine.Mesh.Elements
using ClimateMachine


import Base: getindex
struct LagrangeInterpolant{S}
    values::S
    points::S
    weights::S
    baryweights::S
end

getindex(f::LagrangeInterpolant, i) = getindex(f.values, i)
# From Barycentric Lagrange Interpolation
function (f::LagrangeInterpolant)(x::Number)
    numerator = 0.0
    denominator = 0.0
    for i in eachindex(f.points)
        Δx = x - f.points[i]
        if abs(Δx) ≤ eps(f.points[i]) return f[i] end
        pole = f.baryweights[i] / Δx
        numerator   += pole * f[i]
        denominator += pole
    end
    return numerator / denominator
end

## Example
rx, wx = GaussQuadrature.legendre(npx+1, both)
ωx = baryweights(rx)
ry, wy = GaussQuadrature.legendre(npx+1, both)
ωy = baryweights(ry)
rz, wz = GaussQuadrature.legendre(npx+1, both)
ωz = baryweights(rz)

g = LagrangeInterpolant(copy(rx) .^ 2, rx, wx, ωx)

#=
@btime g.(-1:0.01:1)
@btime [g(x) for x in -1:0.01:1]
=#
Plots.gr(size = (300,300))
h(x) = sin(π*x)
g.values .= @. h(g.points)
newx = collect(-1:0.01:1)
Plots.plot(newx, g.(newx), label = "Interpolated Values")
Plots.plot!(newx, h.(newx), label = "True Values", linestyle = :dash)
Plots.scatter!(g.points, g.values, legend = :bottomright, label = "Interpolation Points")