using TestItems

@testmodule Layer_Potential begin
    using Test
    using AdaptiveMeshSolver
    using Inti
    using StaticArrays
    using LinearAlgebra
    using Gmsh

    function test_layer_potentials_simple()
        @testset "Layer Potentials: Simple Test" begin
            op = Inti.Laplace(; dim=2)

            # compute exact solution
            u = (x) -> x[1]^2 - x[2]^2
            du = (x, normal) -> dot((
                2*x[1],
                -2*x[2]
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

    function test_layer_potentials_simple_exp()
        @testset "Layer Potentials: Simple Exponential Test" begin
            op = Inti.Laplace(; dim=2)

            # compute exact solution
            u = (x) -> exp(x[1])*cos(x[2])
            du = (x, normal) -> dot((
                exp(x[1])*cos(x[2]),
                -exp(x[1])*sin(x[2])
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

    function test_layer_potentials_simple_log()
        @testset "Layer Potentials: Simple Logarithm Test" begin
            op = Inti.Laplace(; dim=2)

            a = 2.0
            b = 0.0
            # compute exact solution
            u = (x) -> log((x[1]-a)^2 + (x[2]-b)^2)
            du = (x, normal) -> dot((
                (2*(x[1]-a)) / ((x[1]-a)^2 + (x[2]-b)^2),
                (2*(x[2]-b)) / ((x[1]-a)^2 + (x[2]-b)^2)
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

    function test_layer_potentials_origin_log()
        @testset "Layer Potentials: Origin Logarithm Test" begin
            op = Inti.Laplace(; dim=2)

            a = 0.0
            b = 0.0
            # compute exact solution
            u = (x) -> log((x[1]-a)^2 + (x[2]-b)^2)
            du = (x, normal) -> dot((
                (2*(x[1]-a)) / ((x[1]-a)^2 + (x[2]-b)^2),
                (2*(x[2]-b)) / ((x[1]-a)^2 + (x[2]-b)^2)
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
end

@testitem "Layer Potentials: Simple Test" setup=[Layer_Potential] begin
    Layer_Potential.test_layer_potentials_simple()
end

@testitem "Layer Potentials: Simple Exponential Test" setup=[Layer_Potential] begin
    Layer_Potential.test_layer_potentials_simple_exp()
end

@testitem "Layer Potentials: Simple Logarithm Test" setup=[Layer_Potential] begin
    Layer_Potential.test_layer_potentials_simple_log()
end

@testitem "Layer Potentials: Origin Logarithm Test" setup=[Layer_Potential] begin
    Layer_Potential.test_layer_potentials_origin_log()
end