using TestItems

@testmodule Volume_Potential begin
    using Test
    using AdaptiveMeshSolver
    using Inti
    using MultivariateSingularIntegrals
    using ClassicalOrthogonalPolynomials
    using LinearAlgebra

    function test_volume_potential_constant_density()
        @testset "Constant Density" begin
            # define function and target points
            forcing(x) = x[1]^2 / 4 + x[2]^2 / 4
            density(x) = 1.0 # laplacian of forcing
            target_points = [(0.0, 0.0), (0.5, 0.0), (1.0, 0.0), (1.5, 0.0)]

            # now define meshes
            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, density)

            # now compute volume potentials
            # quadratures = [AdaptiveMeshSolver.getDomainQuadrature(mesh, 4) for mesh in meshes]
            potentials = AdaptiveMeshSolver.calculateVolumePotential(meshes, density, target_points)
            println(potentials)
            error("hi")
        end
    end

    function test_volume_potential_libraries()
        @testset "Libraries" begin
            # define function and target points
            forcing(x) = x[1]^2 / 4 + x[2]^2 / 4
            density(x) = 1.0 # laplacian of forcing
            density(x, y) = density((x, y))
            
            target = [0.1, 1.1]
            P = ClassicalOrthogonalPolynomials.Legendre()
            p = 4
            grid = ClassicalOrthogonalPolynomials.grid(P, p)
            F = density.(grid, grid')
            C = ClassicalOrthogonalPolynomials.plan_transform(P, (p, p)) * F
            N = MultivariateSingularIntegrals.newtoniansquare(target, p)
            println(N)
            println(C)
            println(N * C)
            println(LinearAlgebra.dot(N, C))
        end
    end
end

@testitem "Volume_Potentials: Constant Density" setup=[Volume_Potential] begin
    Volume_Potential.test_volume_potential_constant_density()
end

@testitem "Volume_Potentials: Libraries" setup=[Volume_Potential] begin
    Volume_Potential.test_volume_potential_libraries()
end