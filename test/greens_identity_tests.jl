using TestItems

@testmodule Greens_Identity begin
    using Test
    using AdaptiveMeshSolver
    using Inti
    using StaticArrays
    using LinearAlgebra
    using Gmsh

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
            AdaptiveMeshSolver.showMeshes(meshes)

            dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, meshes[2])
            boundary = Inti.boundary(dom)
            boundary_mesh = Inti.view(meshes[2], boundary)
            boundary_quadrature = Inti.Quadrature(boundary_mesh; qorder = 6)
            # println(boundary_quadrature[])
            # for q in boundary_quadrature
            #     println(q)
            # end
            # error("HI")
            # create quadratures + get target points
            domain_quadratures = [AdaptiveMeshSolver.getDomainQuadrature(mesh, 4) for mesh in meshes]
            boundary_quadratures = [AdaptiveMeshSolver.getBoundaryQuadrature(mesh, 6) for mesh in meshes[2:end]]
            # println(boundary_quadratures[1])
            target = []
            multiplicative_terms = []
            for (i, quadrature) in enumerate(domain_quadratures)
                points = [(q.coords[1], q.coords[2]) for q in quadrature]
                append!(target, points)
            end
            # target = [(0.5370997288165662, 0.5370997288165662)]
            # target = [(0.2908767452425044, 0.5754898716685045)]
            count = 1
            for (i, quadrature) in enumerate(domain_quadratures)
                l = length(quadrature)
                terms = [i >= count && i < count+l ? :inside : :outside for i in 1:length(target)]
                push!(multiplicative_terms, terms)
                count += l
            end
            # multiplicative_terms[1][1] = :outside
            # multiplicative_terms[2][1] = :inside
            # println(multiplicative_terms)
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
            # println(errors)
            # println(argmax(errors))
            # println(target[argmax(errors)])
            # println(errors[argmax(errors)])
            # println(errors)
            # println(errors)
            # println(maximum(errors))
            AdaptiveMeshSolver.showErrorMesh(meshes, target, errors)
            # AdaptiveMeshSolver.showErrorMesh(meshes, target, relative_diffs)

            # @test 0.0 ≈ calculated_u_val atol=2e-13
        end
    end

    function test_simple_quadtree_greens_third_identity_sin_term()
        @testset "Greens Third Identity, u=sin(kπx)sin(kπy), k=0.25, 0.5, 1.0, 2.0, 4.0" begin
            op = Inti.Laplace(; dim=2)

            # compute exact solution
            u = (x) -> exp(x[1])*cos(x[2]) + log((x[1]-2.0)^2 + x[2]^2)
            du = (x, normal) -> dot((
                exp(x[1])*cos(x[2]) + (2*(x[1]-2)) / ((x[1]-2)^2 + x[2]^2),
                -exp(x[1])*sin(x[2]) + (2*x[2]) / ((x[1]-2)^2 + x[2]^2)
            ), normal)
            laplacian_u = (x) -> 0.0

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
            count = 1
            for (i, quadrature) in enumerate(domain_quadratures)
                l = length(quadrature)
                terms = [i >= count && i < count+l ? :inside : :outside for i in 1:length(target)]
                push!(multiplicative_terms, terms)
                count += l
            end
            quad = boundary_quadratures[1]
            S, D = Inti.single_double_layer(;
                op,
                target,
                source = quad,
                compression = (method = :none, ), 
                correction = (method = :dim, target_location = :inside, maxdist = 1.0)
            )

            γ₀u = map(q -> u(q.coords), quad)
            γ₁u = map(q -> du(q.coords, q.normal), quad)
            potentials = S*γ₁u - D*γ₀u
            errors = abs.(potentials - u.(target))
            replace!(errors, 0.0 => 1e-16) # fix exact errors
            AdaptiveMeshSolver.showErrorMesh(meshes, target, errors)
        end
    end

    function test_simple_quadtree_greens_third_identity_complex_forcing()
        @testset "Greens Third Identity, u=sin(kπx)sin(kπy), k=0.25, 0.5, 1.0, 2.0, 4.0" begin
            println("Starting test!")
            σ = 0.05
            op = Inti.Laplace(; dim=2)
            x_0 = 0.0
            y_0 = 0.85

            # compute exact solution
            u = (x) -> exp(-((x[1]-x_0)^2 + (x[2]-y_0)^2)/(2*σ^2))
            du = (x, normal) -> dot((-(x[1]-x_0)/σ^2 * u(x), -(x[2]-y_0)/σ^2 * u(x)), normal)

            laplacian_u = (x) -> (((x[1]-x_0)^2 + (x[2]-y_0)^2)/σ^4 - 2 / σ^2) * u(x)

            # now create meshes
            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, laplacian_u, false)
            println("Generated mesh!")
            # error("bruh")

            domain_quadratures = [AdaptiveMeshSolver.getDomainQuadrature(mesh, 4) for mesh in meshes]
            boundary_quadratures = [AdaptiveMeshSolver.getBoundaryQuadrature(mesh, 6) for mesh in meshes[2:end]]
            # AdaptiveMeshSolver.showMeshes(meshes)
            # AdaptiveMeshSolver.showMesh(meshes[1], [[-0.5805008802525881, -0.7257491197474122]])
            # error("HI")

            target = Vector{SVector{2, Float64}}()
            multiplicative_terms = []
            for (i, quadrature) in enumerate(domain_quadratures)
                points = [SVector{2, Float64}(q.coords[1], q.coords[2]) for q in quadrature]
                append!(target, points)
            end
            
            count = 1
            for (i, quadrature) in enumerate(domain_quadratures)
                l = length(quadrature)
                terms = [i >= count && i < count+l ? :inside : :outside for i in 1:length(target)]
                push!(multiplicative_terms, terms)
                count += l
            end
            potentials = AdaptiveMeshSolver.calculateVolumePotential(domain_quadratures, meshes, laplacian_u, target, multiplicative_terms, true)

            for quad in boundary_quadratures
                S, D = Inti.single_double_layer(;
                    op,
                    target,
                    source = quad,
                    compression = (method = :none, ), 
                    correction = (method = :dim, target_location = :inside, maxdist = 1.0)
                )

                γ₀u = map(q -> u(q.coords), quad)
                γ₁u = map(q -> du(q.coords, q.normal), quad)
                
                contribution = S*γ₁u - D*γ₀u
                potentials += contribution
            end
            errors = abs.(potentials - u.(target))
            n = length(errors)
            L2_error = sqrt(1/n * sum(errors.^2))
            max_error = maximum(errors)
            @show L2_error
            @show max_error
            # AdaptiveMeshSolver.showErrorMesh(meshes, target, errors)
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