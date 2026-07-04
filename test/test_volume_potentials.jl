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

    function get_greens_identity_integrals(forcing, partial_x_forcing, partial_y_forcing, density, points, meshes, multiplicative_terms)
        integrals = []
        for (idx, point) in enumerate(points)
            boundary_function, domain_function = get_greens_identity_funcs(forcing, partial_x_forcing, partial_y_forcing, density, point)
            boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
                meshes;
                bndry_func = boundary_function, 
                bndry_order = 5
            )
            if any(multiplicative_terms[row][idx] == :inside for row in 1:length(meshes))
                forcing_val = forcing(point)
            elseif any(multiplicative_terms[row][idx] == :on for row in 1:length(meshes))
                forcing_val = 1/2*forcing(point)
            else
                forcing_val = 0
            end
            push!(integrals, forcing_val+boundary_integral)
        end
        return integrals
    end

    function get_quadrature_greens_identity_integrals(forcing, partial_x_forcing, partial_y_forcing, density, target_points, boundary_quadratures)
        integrals = []
        for point in target_points
            boundary_function, domain_function = get_greens_identity_funcs(forcing, partial_x_forcing, partial_y_forcing, density, point)
            point_sum = forcing(point)
            for quadrature in boundary_quadratures
                boundary_integral = Inti.integrate(boundary_function, quadrature)
                point_sum += boundary_integral
            end
            push!(integrals, point_sum)
        end
        # for (idx, point) in enumerate(target_points)
        #     boundary_function, domain_function = get_greens_identity_funcs(forcing, partial_x_forcing, partial_y_forcing, density, point)
        #     if idx > 1
        #         boundary_integral = Inti.integrate(boundary_function, boundary_quads[idx-1])
        #     else
        #         boundary_integral = 0
        #     end
        #     push!(integrals, forcing(point)+boundary_integral)
        # end
        return integrals

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

            # now define meshes
            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, density)
            # AdaptiveMeshSolver.showMeshes(meshes)

            comp_mesh = AdaptiveMeshSolver.createMesh()
            # AdaptiveMeshSolver.showMesh(comp_mesh)           

            order = 5

            domain_quadratures = [AdaptiveMeshSolver.getDomainQuadrature(mesh, order) for mesh in meshes]
            boundary_quadratures = [AdaptiveMeshSolver.getBoundaryQuadrature(mesh, order) for mesh in meshes[2:end]]
            target_points = []
            for quadrature in domain_quadratures
                points = [(q.coords[1], q.coords[2]) for q in quadrature]
                append!(target_points, points)
            end
            # target_points = [(0.6857194332951928, 0.6864982436202831)]

            # now compute volume potentials
            # we can do this easily now that we know the points are quadrature points
            multiplicative_terms = AdaptiveMeshSolver.getMultiplicativeTerm(target_points, meshes)
            potentials = AdaptiveMeshSolver.calculateVolumePotential(meshes, density, target_points, multiplicative_terms)
            # green_potentials = get_quadrature_greens_identity_integrals(forcing, partial_x_forcing, partial_y_forcing, density, target_points, boundary_quadratures)
            triangle_potentials = AdaptiveMeshSolver.calculateTriangleVolumePotential(comp_mesh, density, target_points, [:inside for _ in target_points])
            diffs = abs.(potentials - triangle_potentials)
            AdaptiveMeshSolver.showErrorMesh(meshes, target_points, diffs)
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
            target_points = [(round(0.1*x, digits=1), 0.0) for x in -15:15]

            # now define meshes
            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, density)
            # AdaptiveMeshSolver.showMeshes(meshes)

            # now compute volume potentials
            multiplicative_terms = AdaptiveMeshSolver.getMultiplicativeTerm(target_points, meshes)
            potentials = AdaptiveMeshSolver.calculateVolumePotential(meshes, density, target_points, multiplicative_terms)

            for (i, point) in enumerate(target_points)
                greens_identity_integral = get_greens_identity_integral(forcing, partial_x_forcing, partial_y_forcing, density, point, meshes, multiplicative_terms, i)  
                println(point, " - ", greens_identity_integral - potentials[i])              
                # @test greens_identity_integral ≈ potentials[i] atol=1e-10
            end
        end
    end

    function test_volume_potential_fancy_density()
        @testset "Fancy Density" begin
            # define function and target points
            k = 1.0
            forcing = (x) -> (1-x[1]^2-x[2]^2)^2
            # check these eqns lol
            partial_x_forcing = (x) -> -4*x[1] + 4*x[1]^3 + 4*x[1]*x[2]^2
            partial_y_forcing = (x) -> -4*x[2] + 4*x[2]^3 + 4*x[2]*x[1]^2
            density = (x) -> 16*(x[1]^2+x[2]^2) - 8

            # now define meshes
            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, density)
            AdaptiveMeshSolver.showMeshes(meshes)

            comp_mesh = AdaptiveMeshSolver.createMesh()
            AdaptiveMeshSolver.showMesh(comp_mesh)           

            order = 17

            domain_quadratures = [AdaptiveMeshSolver.getDomainQuadrature(mesh, order) for mesh in meshes]
            boundary_quadratures = [AdaptiveMeshSolver.getBoundaryQuadrature(mesh, order) for mesh in meshes[2:end]]
            target_points = []
            for quadrature in domain_quadratures
                points = [(q.coords[1], q.coords[2]) for q in quadrature]
                append!(target_points, points)
            end

            # now compute volume potentials
            # we can do this easily now that we know the points are quadrature points
            multiplicative_terms = AdaptiveMeshSolver.getMultiplicativeTerm(target_points, meshes)
            potentials = AdaptiveMeshSolver.calculateVolumePotential(meshes, density, target_points, multiplicative_terms)
            # green_potentials = get_quadrature_greens_identity_integrals(forcing, partial_x_forcing, partial_y_forcing, density, target_points, boundary_quadratures)
            triangle_potentials = AdaptiveMeshSolver.calculateTriangleVolumePotential(comp_mesh, density, target_points, [:inside for _ in target_points])
            diffs = abs.(potentials - triangle_potentials)
            AdaptiveMeshSolver.showErrorMesh(meshes, target_points, diffs)
            # for index in 1:length(target_points)
            #     magnitude = round(sqrt(target_points[index][1]^2 + target_points[index][2]^2), digits=3)
            #     # println(magnitude)
            #     if magnitude == 0.734
            #         println("--------")
            #         println("Target point: ", target_points[index])
            #         println("Target point magnitude: ", )
            #         println("Green potential: ", green_potentials[index])
            #         println("Potential: ", potentials[index])
            #         println("Diff: ", diffs[index])
            #         println("--------")
            #     end
            # end
            

            
            
            println(maximum(diffs))
            println(argmax(diffs))
            println(minimum(diffs))
            println(argmin(diffs))
            println(Float64(sum(diffs)) / Float64(length(diffs)))
            # for diff in diffs
            #     @test diff ≈ 0.0 atol=1e-10
            # end
        end
    end

    function test_volume_potential_sine_density()
        @testset "Sine Density" begin
            # define function and target points
            k = 1.0
            forcing = (x) -> sin(k*pi*x[1])*sin(k*pi*x[2])
            # check these eqns lol
            partial_x_forcing = (x) -> k*pi*cos(k*pi*x[1])*sin(k*pi*x[2])
            partial_y_forcing = (x) -> k*pi*sin(k*pi*x[1])*cos(k*pi*x[2])
            density = (x) -> -2*k^2*pi^2*sin(k*pi*x[1])*sin(k*pi*x[2])

            # now define meshes
            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, density)
            # AdaptiveMeshSolver.showMeshes(meshes)

            comp_mesh = AdaptiveMeshSolver.createMesh()
            # AdaptiveMeshSolver.showMesh(comp_mesh)           

            order = 5

            domain_quadratures = [AdaptiveMeshSolver.getDomainQuadrature(mesh, order) for mesh in meshes]
            boundary_quadratures = [AdaptiveMeshSolver.getBoundaryQuadrature(mesh, order) for mesh in meshes[2:end]]
            target_points = []
            for quadrature in domain_quadratures
                points = [(q.coords[1], q.coords[2]) for q in quadrature]
                append!(target_points, points)
            end
            # target_points = [(-0.7046722991094017, -0.0001280194753559916)]
            # target_points = [target_points[2]]
            # println(target_points[1])

            # now compute volume potentials
            # we can do this easily now that we know the points are quadrature points
            multiplicative_terms = AdaptiveMeshSolver.getMultiplicativeTerm(target_points, meshes)
            println("Calculating volume potentials")
            potentials = AdaptiveMeshSolver.calculateVolumePotential(meshes, density, target_points, multiplicative_terms)
            println("Calculating reference potentials")
            # green_potentials = get_quadrature_greens_identity_integrals(forcing, partial_x_forcing, partial_y_forcing, density, target_points, boundary_quadratures)
            triangle_potentials = AdaptiveMeshSolver.calculateTriangleVolumePotential(comp_mesh, density, target_points, [:inside for _ in target_points])
            diffs = abs.((potentials - triangle_potentials))
            relative_diffs = abs.((potentials - triangle_potentials) ./ (abs.(potentials).+1e-9))
            AdaptiveMeshSolver.showErrorMesh(meshes, target_points, diffs)
            AdaptiveMeshSolver.showErrorMesh(meshes, target_points, relative_diffs)
        end
    end

    function test_volume_potential_complex_density()
        @testset "Complex Density" begin
            # define function and target points
            r0 = (0.0, 0.0)
            forcing = (x) -> -2*pi^2 * sin(pi*x[1]) * sin(pi*x[2]) + exp(-100*((x[1]- r0[1])^2 + (x[2] - r0[2])^2)) # forcing function (is it sufficiently smooth?)
            # check these eqns lol
            partial_x_forcing = (x) -> -2*pi^3 * cos(pi*x[1]) * sin(pi*x[2]) - 200*(x[1]-r0[1])*exp(-100*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
            partial_y_forcing = (x) -> -2*pi^3 * sin(pi*x[1]) * cos(pi*x[2]) - 200*(x[2]-r0[2])*exp(-100*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
            density = (x) -> 4*pi^4 * sin(pi*x[1]) * sin(pi*x[2]) + (40000*((x[1]-r0[1])^2 + (x[2]-r0[2])^2)-2400)*exp(-100*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
            

            # now define meshes
            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, density)
            AdaptiveMeshSolver.showMeshes(meshes)

            comp_mesh = AdaptiveMeshSolver.createMesh()
            # AdaptiveMeshSolver.showMesh(comp_mesh)           

            order = 5

            domain_quadratures = [AdaptiveMeshSolver.getDomainQuadrature(mesh, order) for mesh in meshes]
            boundary_quadratures = [AdaptiveMeshSolver.getBoundaryQuadrature(mesh, order) for mesh in meshes[2:end]]
            target_points = []
            for quadrature in domain_quadratures
                points = [(q.coords[1], q.coords[2]) for q in quadrature]
                append!(target_points, points)
            end
            
            multiplicative_terms = AdaptiveMeshSolver.getMultiplicativeTerm(target_points, meshes)
            println("Calculating volume potentials")
            potentials = AdaptiveMeshSolver.calculateVolumePotential(meshes, density, target_points, multiplicative_terms)
            println("Calculating reference potentials")
            # green_potentials = get_quadrature_greens_identity_integrals(forcing, partial_x_forcing, partial_y_forcing, density, target_points, boundary_quadratures)
            triangle_potentials = AdaptiveMeshSolver.calculateTriangleVolumePotential(comp_mesh, density, target_points, [:inside for _ in target_points])
            diffs = abs.((potentials - triangle_potentials))
            relative_diffs = abs.((potentials - triangle_potentials) ./ (abs.(potentials).+1e-9))
            AdaptiveMeshSolver.showErrorMesh(meshes, target_points, diffs)
            AdaptiveMeshSolver.showErrorMesh(meshes, target_points, relative_diffs)
        end
    end

    function test_library()
        @testset "Library test" begin
            f = (x) -> -2*pi^2*sin(pi*x[1])*sin(pi*x[2])

            # points_1 = [[1.0 + 1e-15, 0.1*x] for x in -10:10]
            # points_2 = [[-1.0 + 1e-15, 0.1*x] for x in -10:10]

            points_1 = [[1.0, 0.1*x] for x in -10:10]
            points_2 = [[-1.0, 0.1*x] for x in -10:10]
            
            P = ClassicalOrthogonalPolynomials.Legendre()
            p = 5;
            x = ClassicalOrthogonalPolynomials.grid(P, p)
            coords_1 = [(-0.275, 0.275), (0.0, 0.275), (0.0, 0.550), (-0.275, 0.550)]
            coords_2 = [(0.0, 0.275), (0.275, 0.275), (0.275, 0.550), (0.0, 0.550)]
            h = 0.275
           
            x_tensor_product = [(x1, x2) for x1 in x, x2 in x]
            u_1 = (x) -> h/2 .* (x .+ 1) .+ coords_1[1]
            u_2 = (x) -> h/2 .* (x .+ 1) .+ coords_2[1]
            
            f_u_1 = f.(u_1.(x_tensor_product))
            f_u_2 = f.(u_2.(x_tensor_product))
            
            C_1 = plan_transform(P, (p, p))*f_u_1;
            C_2 = plan_transform(P, (p, p))*f_u_2

            for i in 1:length(points_1)
                N_1 = Float64.(newtoniansquare(big.(points_1[i]), p))
                N_2 = Float64.(newtoniansquare(big.(points_2[i]), p))

                println("-----")
                println(points_1[i])
                println(dot(N_1, C_1))
                println(points_2[i])
                println(dot(N_2, C_2))
            end
            
            
            # print(big.(point))
            # println(N)
            # println(dot(N, C))
            # println(dot(N, C)/(2pi))
        end
    end

    function test_library_2()
        @testset "Library test 2" begin
            k = 1.0
            f = (x) -> -2*k^2*pi^2*sin(k*pi*x[1])*sin(k*pi*x[2])
            point = [-0.51287112518, 0.5071777696040001]
            P = ClassicalOrthogonalPolynomials.Legendre()
            p = 40;
            x = ClassicalOrthogonalPolynomials.grid(P, p)
            coords = [(-0.6875, -0.1375), (-0.55, -0.1375), (-0.55, 5.551115123125783e-17), (-0.6875, 5.551115123125783e-17)]
            h = coords[3][1] - coords[1][1]

            x_tensor_product = [(x1, x2) for x1 in x, x2 in x]
            u = (x) -> h/2 .* (x .+ 1) .+ coords[1]

            f_u = f.(u.(x_tensor_product))
            
            C = plan_transform(P, (p, p))*f_u
            N = Float64.(newtoniansquare(big.(point), p))

            println(dot(N, C))
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

@testitem "Volume Potentials: Fancy Density" setup = [Volume_Potential] begin
    Volume_Potential.test_volume_potential_fancy_density()
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

@testitem "Volume Potentials: Library Test 2" setup=[Volume_Potential] begin
    Volume_Potential.test_library_2()
end