Ω = Circle(-1,1) × Circle(-1,1) × Circle(-1,1)
ex, ey, ez = (3,5,7)[1:ndims(Ω)]
(npx, npy, npz) = (1,1,1)[1:ndims(Ω)]
ClimateMachine.gpu_allowscalar(true)
grid = DiscontinuousSpectralElementGrid(Ω, elements = (ex, ey, ez), polynomialorder = (npx, npy, npz), array = ArrayType)
gridhelper = GridHelper(grid)
x, y, z = coordinates(grid)
ϕ =  ScalarField(copy(x), gridhelper)
@testset "Check Interpolation 3D" begin
    ϕ .= x
    @test ϕ(0.1, 0.2, 0.3) ≈ 0.1
    ϕ .= y
    @test ϕ(0.1, 0.2, 0.3) ≈ 0.2
    ϕ .= z
    @test ϕ(0.1, 0.2, 0.3) ≈ 0.3
end
##
Ω = Circle(0,2) × Circle(-4,2)
elements =  (3,5,7)[1:ndims(Ω)]
polynomialorder = (1,1,1)[1:ndims(Ω)]
ClimateMachine.gpu_allowscalar(true)
grid = DiscontinuousSpectralElementGrid(Ω, elements = elements, polynomialorder = polynomialorder, array = ArrayType)
gridhelper = GridHelper(grid)
x, y, z = coordinates(grid)
ϕ =  ScalarField(copy(x), gridhelper)
@testset "Check Interpolation 2D" begin
    ϕ .= x
    @test ϕ(0.1, 0.2) ≈ 0.1
    ϕ .= y
    @test ϕ(0.1, 0.2) ≈ 0.2
end