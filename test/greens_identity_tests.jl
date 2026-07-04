using TestItems

@testmodule Greens_Identity begin
    using Test
    using AdaptiveMeshSolver
    using Inti
    using StaticArrays
    using LinearAlgebra

    function test_simple_quadtree_greens_third_identity_simple_forcing()
        @testset "Greens Third Identity, u=1" begin
            x_test = (1, 1) # point just outside of the region
            u = (x) -> 1.0
            # check these eqns lol
            partial_x_u = (x) -> 0.0
            partial_y_u = (x) -> 0.0
            laplacian_u = (x) -> 0.0

            r0 = (10, 10)

            greens_fn = (r, x) -> 1/(2pi) * log(AdaptiveMeshSolver.distance(x, r))
            partial_x_greens_fn = (r, x) -> 1/(2pi) * ((x[1]-r[1]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))
            partial_y_greens_fn = (r, x) -> 1/(2pi) * ((x[2]-r[2]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))

            normal_derivative_u = (x, n) -> partial_x_u(x) * n[1] + partial_y_u(x) * n[2]
            normal_derivative_greens_fn = (r, x, n) -> partial_x_greens_fn(r, x) * n[1] + partial_y_greens_fn(r, x) * n[2]

            boundary_function = (q) -> (greens_fn(x_test, q.coords) * normal_derivative_u(q.coords, q.normal) - 
                                        u(q.coords) * normal_derivative_greens_fn(x_test, q.coords, q.normal))

            domain_function = (q) -> (greens_fn(x_test, q.coords) * laplacian_u(q.coords))

            forcing_func = (x) -> greens_fn(x_test, x) * laplacian_u(x)

            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_func)

            AdaptiveMeshSolver.showMeshes(meshes)

            println("Done creating meshes!")
            domain_integral, boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
                meshes; 
                dom_func = domain_function,
                bndry_func = boundary_function, 
                dom_order = 4,
                bndry_order = 12
            )
            calculated_u_val = abs(domain_integral - boundary_integral)

            @test 0.0 ≈ calculated_u_val atol=1e-13
        end
    end

    function test_simple_quadtree_greens_third_identity_linear_x()
        @testset "Greens Third Identity, u=x" begin
            x_test = (1, 1) # point just outside of the region
            r0 = (10, 10)
            u = (x) -> x[1]
            # check these eqns lol
            partial_x_u = (x) -> 1.0
            partial_y_u = (x) -> 0.0
            laplacian_u = (x) -> 0.0
            greens_fn = (r, x) -> 1/(2pi) * log(AdaptiveMeshSolver.distance(x, r))
            partial_x_greens_fn = (r, x) -> 1/(2pi) * ((x[1]-r[1]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))
            partial_y_greens_fn = (r, x) -> 1/(2pi) * ((x[2]-r[2]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))

            normal_derivative_u = (x, n) -> partial_x_u(x) * n[1] + partial_y_u(x) * n[2]
            normal_derivative_greens_fn = (r, x, n) -> partial_x_greens_fn(r, x) * n[1] + partial_y_greens_fn(r, x) * n[2]

            boundary_function = (q) -> (greens_fn(x_test, q.coords) * normal_derivative_u(q.coords, q.normal) - 
                                        u(q.coords) * normal_derivative_greens_fn(x_test, q.coords, q.normal))

            domain_function = (q) -> (greens_fn(x_test, q.coords) * laplacian_u(q.coords))

            forcing_func = (x) -> greens_fn(x_test, x) * laplacian_u(x)

            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_func)

            domain_integral, boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
                meshes, r0; 
                dom_func = domain_function,
                bndry_func = boundary_function, 
                dom_order = 4,
                bndry_order = 12
            )

            calculated_u_val = abs(domain_integral - boundary_integral)
            println(calculated_u_val)

            @test 0.0 ≈ calculated_u_val atol=1e-13
        end
    end

    function test_simple_quadtree_greens_third_identity_linear_y()
        @testset "Greens Third Identity, u=y" begin
            x_test = (1, 1) # point just outside of the region
            r0 = (10, 10)
            u = (x) -> x[2]
            # check these eqns lol
            partial_x_u = (x) -> 0.0
            partial_y_u = (x) -> 1.0
            laplacian_u = (x) -> 0.0
            greens_fn = (r, x) -> 1/(2pi) * log(AdaptiveMeshSolver.distance(x, r))
            partial_x_greens_fn = (r, x) -> 1/(2pi) * ((x[1]-r[1]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))
            partial_y_greens_fn = (r, x) -> 1/(2pi) * ((x[2]-r[2]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))

            normal_derivative_u = (x, n) -> partial_x_u(x) * n[1] + partial_y_u(x) * n[2]
            normal_derivative_greens_fn = (r, x, n) -> partial_x_greens_fn(r, x) * n[1] + partial_y_greens_fn(r, x) * n[2]

            boundary_function = (q) -> (greens_fn(x_test, q.coords) * normal_derivative_u(q.coords, q.normal) - 
                                        u(q.coords) * normal_derivative_greens_fn(x_test, q.coords, q.normal))

            domain_function = (q) -> (greens_fn(x_test, q.coords) * laplacian_u(q.coords))

            forcing_func = (x) -> greens_fn(x_test, x) * laplacian_u(x)

            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_func)

            domain_integral, boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
                meshes, r0; 
                dom_func = domain_function,
                bndry_func = boundary_function, 
                dom_order = 4,
                bndry_order = 12
            )

            calculated_u_val = abs(domain_integral - boundary_integral)

            @test 0.0 ≈ calculated_u_val atol=1e-13
        end
    end

    function test_simple_quadtree_greens_third_identity_quadratic()
        @testset "Greens Third Identity, u=x^2+y^2" begin
            op = Inti.Laplace(; dim=2)

            # compute exact solution
            u = (x) -> x[1]^2 + x[2]^2
            du = (x, normal) -> dot((2*x[1], 2*x[2]), normal)

            laplacian_u = (x) -> 4.0

            # now create meshes
            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, laplacian_u)
            # AdaptiveMeshSolver.showMeshes(meshes)

            # create quadratures + get target points
            domain_quadratures = [AdaptiveMeshSolver.getDomainQuadrature(mesh, 4) for mesh in meshes]
            boundary_quadratures = [AdaptiveMeshSolver.getBoundaryQuadrature(mesh, 6) for mesh in meshes[2:end]]
            target = []
            multiplicative_terms = []
            for (i, quadrature) in enumerate(domain_quadratures)
                points = [(q.coords[1], q.coords[2]) for q in quadrature]
                append!(target, points)
            end
            # target = [(0.5370997288165662, 0.5370997288165662)]
            target = [(0.2908767452425044, 0.5754898716685045)]
            count = 1
            for (i, quadrature) in enumerate(domain_quadratures)
                l = length(quadrature)
                terms = [i >= count && i < count+l ? :inside : :outside for i in 1:length(target)]
                push!(multiplicative_terms, terms)
                count += l
            end
            multiplicative_terms[1][1] = :outside
            multiplicative_terms[2][1] = :inside
            println(multiplicative_terms)
            # multiplicative_terms = AdaptiveMeshSolver.getMultiplicativeTerm(target, meshes)
            potentials = AdaptiveMeshSolver.calculateVolumePotential(domain_quadratures, meshes, laplacian_u, target, multiplicative_terms)

            for quad in boundary_quadratures
                S, D = Inti.single_double_layer(;
                    op,
                    target,
                    source = quad,
                    compression = (method = :none, ), 
                    correction = (method = :dim, target_location = :inside, maxdist = 0.4)
                )

                γ₀u = map(q -> u(q.coords), quad)
                γ₁u = map(q -> du(q.coords, q.normal), quad)
                
                contribution = S*γ₁u - D*γ₀u
                potentials += contribution
            end
            errors = abs.(potentials - u.(target))
            replace!(errors, 0.0 => 1e-16) # fix exact errors
            # println(multiplicative_terms)
            println(errors)
            # println(argmax(errors))
            # println(target[argmax(errors)])
            # println(errors[argmax(errors)])
            # println(errors)
            # println(errors)
            # AdaptiveMeshSolver.showErrorMesh(meshes, target, errors)
            # AdaptiveMeshSolver.showErrorMesh(meshes, target, relative_diffs)

            # @test 0.0 ≈ calculated_u_val atol=2e-13
        end
    end

    function test_simple_quadtree_greens_third_identity_sin_term()
        @testset "Greens Third Identity, u=sin(kπx)sin(kπy), k=0.25, 0.5, 1.0, 2.0, 4.0" begin
            for k in [0.25, 0.5, 1.0, 2.0, 4.0]
                x_test = (1, 1) # point just outside of the region
                r0 = (10, 10)
                u = (x) -> sin(k*pi*x[1])*sin(k*pi*x[2])
                # check these eqns lol
                partial_x_u = (x) -> k*pi*cos(k*pi*x[1])*sin(k*pi*x[2])
                partial_y_u = (x) -> k*pi*sin(k*pi*x[1])*cos(k*pi*x[2])
                laplacian_u = (x) -> -2*k^2*pi^2*sin(k*pi*x[1])*sin(k*pi*x[2])
                greens_fn = (r, x) -> 1/(2pi) * log(AdaptiveMeshSolver.distance(x, r))
                partial_x_greens_fn = (r, x) -> 1/(2pi) * ((x[1]-r[1]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))
                partial_y_greens_fn = (r, x) -> 1/(2pi) * ((x[2]-r[2]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))

                normal_derivative_u = (x, n) -> partial_x_u(x) * n[1] + partial_y_u(x) * n[2]
                normal_derivative_greens_fn = (r, x, n) -> partial_x_greens_fn(r, x) * n[1] + partial_y_greens_fn(r, x) * n[2]

                boundary_fn1 = (q) -> greens_fn(x_test, q.coords) * normal_derivative_u(q.coords, q.normal)
                boundary_fn2 = (q) -> u(q.coords) * normal_derivative_greens_fn(x_test, q.coords, q.normal)

                boundary_function = (q) -> boundary_fn1(q) - boundary_fn2(q)

                domain_function = (q) -> (greens_fn(x_test, q.coords) * laplacian_u(q.coords))

                forcing_func = (x) -> greens_fn(x_test, x) * laplacian_u(x)

                parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
                meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_func)

                domain_integral, boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
                    meshes, r0; 
                    dom_func = domain_function,
                    bndry_func = boundary_function, 
                    dom_order = 17,
                    bndry_order = 12
                )

                calculated_u_val = abs(domain_integral - boundary_integral)

                @test 0.0 ≈ calculated_u_val atol=1e-13
            end
        end
    end

    function test_simple_quadtree_greens_third_identity_complex_forcing()
        @testset "Greens Third Identity, complex u" begin
            x_test = (1, 1) # point just outside of the region
            r0 = (-0.5, 0.0)
            u = (x) -> -2*pi^2 * sin(pi*x[1]) * sin(pi*x[2]) + exp(-600*((x[1]- r0[1])^2 + (x[2] - r0[2])^2)) # forcing function (is it sufficiently smooth?)
            # check these eqns lol
            partial_x_u = (x) -> -2*pi^3 * cos(pi*x[1]) * sin(pi*x[2]) - 1200*(x[1]-r0[1])*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
            partial_y_u = (x) -> -2*pi^3 * sin(pi*x[1]) * cos(pi*x[2]) - 1200*(x[2]-r0[2])*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
            laplacian_u = (x) -> 4*pi^4 * sin(pi*x[1]) * sin(pi*x[2]) + (1440000*((x[1]-r0[1])^2 + (x[2]-r0[2])^2)-2400)*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
            greens_fn = (r, x) -> 1/(2pi) * log(AdaptiveMeshSolver.distance(x, r))
            partial_x_greens_fn = (r, x) -> 1/(2pi) * ((x[1]-r[1]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))
            partial_y_greens_fn = (r, x) -> 1/(2pi) * ((x[2]-r[2]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))

            normal_derivative_u = (x, n) -> partial_x_u(x) * n[1] + partial_y_u(x) * n[2]
            normal_derivative_greens_fn = (r, x, n) -> partial_x_greens_fn(r, x) * n[1] + partial_y_greens_fn(r, x) * n[2]

            boundary_function = (q) -> (greens_fn(x_test, q.coords) * normal_derivative_u(q.coords, q.normal) - 
                                        u(q.coords) * normal_derivative_greens_fn(x_test, q.coords, q.normal))

            domain_function = (q) -> (greens_fn(x_test, q.coords) * laplacian_u(q.coords))

            forcing_func = (x) -> greens_fn(x_test, x) * laplacian_u(x)

            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_func)

            
            domain_integral, boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
                meshes, r0; 
                dom_func = domain_function,
                bndry_func = boundary_function, 
                dom_order = 17,
                bndry_order = 12
            )

            # println(domain_integral)
            # println(boundary_integral)

            calculated_u_val = abs(domain_integral - boundary_integral)

            # println("Boundary integral = ", boundary_integral)
            # println("Domain Integral = ", domain_integral)
            # println("Difference = ", calculated_u_val)
            @test 0.0 ≈ calculated_u_val atol=1e-12
        end
    end
end

@testitem "Greens_Identity: Greens Theorem Area" setup=[Greens_Identity] begin
    Greens_Identity.test_simple_quadtree_greens_theorem()
end

@testitem "Greens_Identity: Greens Third Identity, u=1" setup=[Greens_Identity] begin
    Greens_Identity.test_simple_quadtree_greens_third_identity_simple_forcing()
end

@testitem "Greens_Identity: Greens Third Identity, u=x" setup=[Greens_Identity] begin
    Greens_Identity.test_simple_quadtree_greens_third_identity_linear_x()
end

@testitem "Greens_Identity: Greens Third Identity, u=y" setup=[Greens_Identity] begin
    Greens_Identity.test_simple_quadtree_greens_third_identity_linear_y()
end

@testitem "Greens_Identity: Greens Third Identity, u=x^2+y^2" setup=[Greens_Identity] begin
    Greens_Identity.test_simple_quadtree_greens_third_identity_quadratic()
end

@testitem "Greens_Identity: Greens Third Identity, u=sin(kπx)sin(kπy), k=0.25, 0.5, 1.0, 2.0, 4.0" setup=[Greens_Identity] begin
    Greens_Identity.test_simple_quadtree_greens_third_identity_sin_term()
end

@testitem "Greens_Identity: Greens Third Identity, complex u" setup=[Greens_Identity] begin
    Greens_Identity.test_simple_quadtree_greens_third_identity_complex_forcing()
end