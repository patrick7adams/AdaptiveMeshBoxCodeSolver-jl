using TestItems

@testmodule Triangle_Refinement begin
    using Test
    using AdaptiveMeshSolver
    using Inti
    using StaticArrays
    using LinearAlgebra
    using Gmsh

    function test_triangle_refinement()
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
            multiplicative_terms[1][1] = :outside
            multiplicative_terms[2][1] = :inside
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
            AdaptiveMeshSolver.showErrorMesh(meshes, target, errors)
            # AdaptiveMeshSolver.showErrorMesh(meshes, target, relative_diffs)

            # @test 0.0 ≈ calculated_u_val atol=2e-13
        end
    end
end

@testitem "Triangle_Refinement: Greens Theorem Area" setup=[Triangle_Refinement] begin
    Triangle_Refinement.test_simple_quadtree_greens_theorem()
end