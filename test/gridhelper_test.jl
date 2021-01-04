Ω = Circle(-1,1) × Circle(-1,1) × Circle(-1,1)
dims = ndims(Ω)
ex, ey, ez = (3,4,5) .* 10
(npx, npy, npz) = (1,2,3) # needs to be at least these polynomial orders due to last test
ClimateMachine.gpu_allowscalar(true)
# Define Grid: might want to loop over element sizes and polynomial orders
grid = DiscontinuousSpectralElementGrid(Ω, elements = (ex, ey, ez), polynomialorder = (npx, npy, npz), array = ArrayType)
gridhelper = GridHelper(grid)
addup(xC, tol) = sum(abs.(xC[1] .- xC) .≤ tol)
npx, npy, npz = polynomialorders(grid)
x, y, z = coordinates(grid)
xC, yC, zC = cellcenters(grid)

ne = size(x)[2]
ex = round(Int64, ne / addup(xC, 10^4 * eps(maximum(abs.(x)))))
ey = round(Int64, ne / addup(yC, 10^4 * eps(maximum(abs.(y)))))
ez = round(Int64, ne / addup(zC, 10^4 * eps(maximum(abs.(z)))))
check = ne == ex * ey * ez
check ? true : error("improper counting")
p = getperm(xC, yC, zC, ex, ey, ez)
lin = reshape(collect(1:length(xC)), (ex, ey, ez))


ne = size(x)[2]
rx, wx = GaussQuadrature.legendre(npx+1, both)
ωx = baryweights(rx)
ry, wy = GaussQuadrature.legendre(npy+1, both)
ωy = baryweights(ry)
rz, wz = GaussQuadrature.legendre(npz+1, both)
ωz = baryweights(rz)
linlocal = reshape(collect(1:(npx+1)*(npy+1)*(npz+1)), (npx+1, npy+1, npz+1))

# computation


# now test
@testset "Check getvalue Interpolation" begin
    f = copy(z)
    fl = reshape(f, (npx+1, npy+1, npz+1, ex*ey*ez))
    location = (0.0, 0.0, -0.3)
    @test location[3] ≈ getvalue(fl, xC, yC, zC, location, p, lin, linlocal, x, y, z, rx, ry, rz, ωx, ωy, ωz)
    @. f = x
    fl = reshape(f, (npx+1, npy+1, npz+1, ex*ey*ez))
    location = (0.1, 0.0, -0.3)
    @test location[1] ≈ getvalue(fl, xC, yC, zC, location, p, lin, linlocal, x, y, z, rx, ry, rz, ωx, ωy, ωz)
    @. f = y
    fl = reshape(f, (npx+1, npy+1, npz+1, ex*ey*ez))
    location = (0.1, -0.9, -0.3)
    @test location[2] ≈ getvalue(fl, xC, yC, zC, location, p, lin, linlocal, x, y, z, rx, ry, rz, ωx, ωy, ωz)
end
##
@testset "Check GridHelper Interpolation" begin
    f = copy(z)
    location = (0.0, 0.0, -0.3)
    @test location[3] ≈ getvalue(f, location, gridhelper)
    @. f = x
    location = (0.1, 0.0, -0.3)
    @test location[1] ≈ getvalue(f, location, gridhelper)
    @. f = y
    location = (0.1, -0.9, -0.3)
    @test location[2] ≈ getvalue(f, location, gridhelper)
    @. f = x + y^2 + z^3
    location = (0.1, -0.9, -0.3)
    value = location[1] + location[2]^2 + location[3]^3
    @test value ≈ getvalue(f, location, gridhelper)
end