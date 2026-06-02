using TestItems

@testmodule Volume_Potential begin
    using Test
    using AdaptiveMeshSolver
    using Inti

    function test_volume_potential_constant_density()
        @testset "Constant Density" begin
            # define function and target points
            forcing(x) = x[1]^2 / 4 + x[2]^2 / 4
            density(x) = 1.0 # laplacian of forcing
            target_points = [(0.0, 0.0), (0.5, 0.0), (1.0, 0.0), (1.5, 0.0)]
            op = Inti.Laplace(; dim = 2)

            # now define meshes
            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            forest, meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, density)

            # now compute volume potentials
            quadratures = [AdaptiveMeshSolver.getDomainQuadrature(mesh, 4) for mesh in meshes]
            potentials = AdaptiveMeshSolver.calculateVolumePotential(forest, quadratures, density, target_points, op)
            println(potentials)
            error("hi")
        end
    end
end

@testitem "Volume_Potentials: Constant Density" setup=[Volume_Potential] begin
    Volume_Potential.test_volume_potential_constant_density()
end