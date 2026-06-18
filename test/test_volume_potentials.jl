using TestItems

@testmodule Volume_Potential begin
    using Test
    using AdaptiveMeshSolver
    using Inti
    using MultivariateSingularIntegrals
    using ClassicalOrthogonalPolynomials
    using LinearAlgebra

    function get_greens_identity_funcs(forcing, partial_x_forcing, partial_y_forcing, density, test_point)
        greens_fn = (r, x) -> 1/(2pi) * log(AdaptiveMeshSolver.distance(x, r))
        partial_x_greens_fn = (r, x) -> 1/(2pi) * ((x[1]-r[1]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))
        partial_y_greens_fn = (r, x) -> 1/(2pi) * ((x[2]-r[2]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))

        normal_derivative_u = (x, n) -> partial_x_forcing(x) * n[1] + partial_y_forcing(x) * n[2]
        normal_derivative_greens_fn = (r, x, n) -> partial_x_greens_fn(r, x) * n[1] + partial_y_greens_fn(r, x) * n[2]

        boundary_function = (q) -> (greens_fn(test_point, q.coords) * normal_derivative_u(q.coords, q.normal) - 
                                    forcing(q.coords) * normal_derivative_greens_fn(test_point, q.coords, q.normal))

        domain_function = (q) -> (greens_fn(test_point, q.coords) * density(q.coords))
        return (boundary_function, domain_function)
    end

    function get_greens_identity_integral(forcing, partial_x_forcing, partial_y_forcing, density, point, meshes, multiplicative_terms, idx)
        boundary_function, domain_function = get_greens_identity_funcs(forcing, partial_x_forcing, partial_y_forcing, density, point)
        boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
            meshes;
            bndry_func = boundary_function, 
            bndry_order = 12
        )
        # println(multiplicative_terms)
        if any(multiplicative_terms[row][idx] == :inside for row in 1:length(meshes))
            forcing_val = forcing(point)
        elseif any(multiplicative_terms[row][idx] == :on for row in 1:length(meshes))
            forcing_val = 1/2*forcing(point)
        else
            forcing_val = 0
        end
        
        return forcing_val+boundary_integral
    end

    function test_triangle_mesh_volume_potential()
        @testset "Triangle mesh" begin
            forcing(x) = x[1]^2 / 4 + x[2]^2 / 4
            partial_x_forcing(x) = x[1] / 2
            partial_y_forcing(x) = x[2] / 2
            density(x) = 1.0 # laplacian of forcing
            # target_points = [(0.0, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0), (1.5, 0.0)]
            target_points = [(0.1, 0.1)]
            # now define meshes
            mesh = AdaptiveMeshSolver.createMesh()

            multiplicative_terms = AdaptiveMeshSolver.getMultiplicativeTerm(target_points, mesh)
            
            # AdaptiveMeshSolver.showMesh(mesh)
            triangle_potentials = AdaptiveMeshSolver.calculateTriangleVolumePotential(mesh, density, target_points, multiplicative_terms)
            println(multiplicative_terms)
            for (i, point) in enumerate(target_points)
                greens_identity_integral = get_greens_identity_integral(forcing, partial_x_forcing, partial_y_forcing, density, point, meshes, multiplicative_terms, i)
                @test greens_identity_integral ≈ triangle_potentials[i] atol=1e-10
            end
        end
    end

    function test_volume_potential_one_density()
        @testset "Constant Density" begin
            # define function and target points
            forcing(x) = x[1]^2/4 + x[2]^2/4
            partial_x_forcing(x) = x[1] / 2
            partial_y_forcing(x) = x[2] / 2
            density(x) = 1.0 # laplacian of forcing
            target_points = [(0.0, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0), (1.5, 0.0)]

            # now define meshes
            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, density)
            AdaptiveMeshSolver.showMeshes(meshes)

            # now compute volume potentials
            multiplicative_terms = AdaptiveMeshSolver.getMultiplicativeTerm(target_points, meshes)
            potentials = AdaptiveMeshSolver.calculateVolumePotential(meshes, density, target_points, multiplicative_terms)

            for (i, point) in enumerate(target_points)
                greens_identity_integral = get_greens_identity_integral(forcing, partial_x_forcing, partial_y_forcing, density, point, meshes, multiplicative_terms, i)                
                # println(greens_identity_integral - potentials[i])
                @test greens_identity_integral ≈ potentials[i] atol=1e-10
            end
        end
    end

    function test_volume_potential_constant_density()
        @testset "Constant Density" begin
            # define function and target points
            forcing(x) = x[1]^2 + x[2]^2
            partial_x_forcing(x) = x[1] * 2
            partial_y_forcing(x) = x[2] * 2
            density(x) = 4.0 # laplacian of forcing
            target_points = [(0.0, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0), (1.5, 0.0)]

            # now define meshes
            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, density)
            AdaptiveMeshSolver.showMeshes(meshes)

            # now compute volume potentials
            multiplicative_terms = AdaptiveMeshSolver.getMultiplicativeTerm(target_points, meshes)
            potentials = AdaptiveMeshSolver.calculateVolumePotential(meshes, density, target_points, multiplicative_terms)

            for (i, point) in enumerate(target_points)
                greens_identity_integral = get_greens_identity_integral(forcing, partial_x_forcing, partial_y_forcing, density, point, meshes, multiplicative_terms, i)                
                # println(greens_identity_integral - potentials[i])
                println(greens_identity_integral)
                @test greens_identity_integral ≈ potentials[i] atol=1e-10
            end
        end
    end

    function test_volume_potential_linear_density()
        @testset "Linear Density" begin
            # define function and target points
            forcing(x) = x[1]^3
            partial_x_forcing(x) = 3*x[1]^2
            partial_y_forcing(x) = 0.0
            density(x) = 6*x[1] # laplacian of forcing
            target_points = [(0.0, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0), (1.5, 0.0)]
            # target_points = [(0.5, 0.0)]

            # now define meshes
            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, density)
            # AdaptiveMeshSolver.showMeshes(meshes)

            # now compute volume potentials
            multiplicative_terms = AdaptiveMeshSolver.getMultiplicativeTerm(target_points, meshes)
            potentials = AdaptiveMeshSolver.calculateVolumePotential(meshes, density, target_points, multiplicative_terms)

            for (i, point) in enumerate(target_points)
                greens_identity_integral = get_greens_identity_integral(forcing, partial_x_forcing, partial_y_forcing, density, point, meshes, multiplicative_terms, i)                
                println(point)                
                println(greens_identity_integral - potentials[i])
                @test greens_identity_integral ≈ potentials[i] atol=1e-10
            end
        end
    end

    function test_volume_potential_sine_density()
        @testset "Complex Density" begin
            # define function and target points
            k = 1.0
            forcing = (x) -> sin(k*pi*x[1])*sin(k*pi*x[2])
            # check these eqns lol
            partial_x_forcing = (x) -> k*pi*cos(k*pi*x[1])*sin(k*pi*x[2])
            partial_y_forcing = (x) -> k*pi*sin(k*pi*x[1])*cos(k*pi*x[2])
            density = (x) -> -2*k^2*pi^2*sin(k*pi*x[1])*sin(k*pi*x[2])
            target_points = [(0.0, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0), (1.5, 0.0)]
            # target_points = [(0.1*i, 0.0) for i in -15:15]
            target_points = [(0.99, 0.0)]
            target_points = [(0.3, 0.4), (0.3, -0.4)]

            # now define meshes
            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, density)
            AdaptiveMeshSolver.showMeshes(meshes)

            # now compute volume potentials
            multiplicative_terms = AdaptiveMeshSolver.getMultiplicativeTerm(target_points, meshes)
            potentials = AdaptiveMeshSolver.calculateVolumePotential(meshes, density, target_points, multiplicative_terms)
            println(multiplicative_terms)
            for (i, point) in enumerate(target_points)
                greens_identity_integral = get_greens_identity_integral(forcing, partial_x_forcing, partial_y_forcing, density, point, meshes, multiplicative_terms, i)
                println("------")
                println(point)                
                println(greens_identity_integral - potentials[i])
                @test greens_identity_integral ≈ potentials[i] atol=1e-10
            end
        end
    end

    function test_volume_potential_complex_density()
        @testset "Complex Density" begin
            # define function and target points
            r0 = (-0.5, 0.0)
            forcing = (x) -> -2*pi^2 * sin(pi*x[1]) * sin(pi*x[2]) + exp(-600*((x[1]- r0[1])^2 + (x[2] - r0[2])^2)) # forcing function (is it sufficiently smooth?)
            # check these eqns lol
            partial_x_forcing = (x) -> -2*pi^3 * cos(pi*x[1]) * sin(pi*x[2]) - 1200*(x[1]-r0[1])*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
            partial_y_forcing = (x) -> -2*pi^3 * sin(pi*x[1]) * cos(pi*x[2]) - 1200*(x[2]-r0[2])*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
            density = (x) -> 4*pi^4 * sin(pi*x[1]) * sin(pi*x[2]) + (1440000*((x[1]-r0[1])^2 + (x[2]-r0[2])^2)-2400)*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
            target_points = [(0.0, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0), (1.5, 0.0)]
            target_points = [(0.1*i, 0.0) for i in -15:15]

            # now define meshes
            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, density)
            AdaptiveMeshSolver.showMeshes(meshes)

            # now compute volume potentials
            multiplicative_terms = AdaptiveMeshSolver.getMultiplicativeTerm(target_points, meshes)
            potentials = AdaptiveMeshSolver.calculateVolumePotential(meshes, density, target_points, multiplicative_terms)

            for (i, point) in enumerate(target_points)
                greens_identity_integral = get_greens_identity_integral(forcing, partial_x_forcing, partial_y_forcing, density, point, meshes, multiplicative_terms, i)
                println("--------")
                println(point)                
                println(greens_identity_integral - potentials[i])
                # @test greens_identity_integral ≈ potentials[i] atol=1e-10
            end
        end
    end

    function test_library()
        @testset "Library test" begin
            f = (x, y) -> 6*x
            point = [0.635925620438, 1.0]
            P = ClassicalOrthogonalPolynomials.Legendre()
            p = 5;
            x = ClassicalOrthogonalPolynomials.grid(P, p)
            F = f.(x, x')
            C = plan_transform(P, (p, p))*F;
            print(F)
            N = Float64.(newtoniansquare(big.(point), p))
            println(dot(N, C)/(2pi))
        end
    end
end


@testitem "Volume_Potentials: Triangle mesh" setup=[Volume_Potential] begin
    Volume_Potential.test_triangle_mesh_volume_potential()
end

@testitem "Volume_Potentials: One Density" setup=[Volume_Potential] begin
    Volume_Potential.test_volume_potential_one_density()
end

@testitem "Volume_Potentials: Constant Density" setup=[Volume_Potential] begin
    Volume_Potential.test_volume_potential_constant_density()
end

@testitem "Volume_Potentials: Linear Density" setup=[Volume_Potential] begin
    Volume_Potential.test_volume_potential_linear_density()
end

@testitem "Volume Potentials: Sine Density" setup = [Volume_Potential] begin
    Volume_Potential.test_volume_potential_sine_density()
end

@testitem "Volume Potentials: Complex Density" setup = [Volume_Potential] begin
    Volume_Potential.test_volume_potential_complex_density()
end

@testitem "Volume Potentials: Library Test" setup=[Volume_Potential] begin
    Volume_Potential.test_library()
end