using TestItems

@testmodule Greens_Theorem begin
    using Test
    using AdaptiveMeshSolver
    using Inti
    using StaticArrays

    function test_quadtree_greens_theorem_constants()
        @testset "Greens Theorem Constants" begin
            P = (x, y) -> 1
            Q = (x, y) -> 0
            dPdx = (x, y) -> 0
            dQdy = (x, y) -> 0

            boundary_fn = (q) -> P(q.coords...)*q.normal[1] + Q(q.coords...)*q.normal[2]
            domain_fn = (q) -> dQdy(q.coords...) + dPdx(q.coords...)

            forcing_function = (x) -> dQdy(x...) + dPdx(x...)

            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_function)

            domain_quads = []
            boundary_quads = []
            for (i, mesh) in enumerate(meshes)
                push!(domain_quads, AdaptiveMeshSolver.getDomainQuadrature(mesh, 4))
                if i != 1
                    push!(boundary_quads, AdaptiveMeshSolver.getBoundaryQuadrature(mesh, 6))
                end
            end
            domain_integral = 0.0
            for quad in domain_quads
                domain_integral += Inti.integrate(domain_fn, quad)
            end
            boundary_integral = 0.0
            for quad in boundary_quads
                boundary_integral += Inti.integrate(boundary_fn, quad)
            end
            calculated_u_val = abs(domain_integral - boundary_integral)
            println(domain_integral)
            println(boundary_integral)
            println(calculated_u_val)
            @test 0.0 ≈ calculated_u_val atol=1e-14
        end
    end

    function test_quadtree_greens_theorem_linear()
        @testset "Greens Theorem Area" begin
            P = (x, y) -> x
            Q = (x, y) -> y
            dPdx = (x, y) -> 1
            dQdy = (x, y) -> 1

            boundary_fn = (q) -> P(q.coords...)*q.normal[1] + Q(q.coords...)*q.normal[2]
            domain_fn = (q) -> dQdy(q.coords...) + dPdx(q.coords...)

            forcing_function = (x) -> dQdy(x...) + dPdx(x...)

            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_function)

            domain_quads = []
            boundary_quads = []
            for (i, mesh) in enumerate(meshes)
                push!(domain_quads, AdaptiveMeshSolver.getDomainQuadrature(mesh, 4))
                if i != 1
                    push!(boundary_quads, AdaptiveMeshSolver.getBoundaryQuadrature(mesh, 6))
                end
            end
            domain_integral = 0.0
            for quad in domain_quads
                domain_integral += Inti.integrate(domain_fn, quad)
            end
            boundary_integral = 0.0
            for quad in boundary_quads
                boundary_integral += Inti.integrate(boundary_fn, quad)
            end
            calculated_u_val = abs(domain_integral - boundary_integral)
            println(domain_integral)
            println(boundary_integral)
            println(calculated_u_val)
            @test 0.0 ≈ calculated_u_val atol=1e-14
        end
    end

    function test_quadtree_greens_theorem_linear_x()
        @testset "Greens Theorem Area" begin
            P = (x, y) -> x
            Q = (x, y) -> 0
            dPdx = (x, y) -> 1
            dQdy = (x, y) -> 0

            boundary_fn = (q) -> P(q.coords...)*q.normal[1] + Q(q.coords...)*q.normal[2]
            domain_fn = (q) -> dQdy(q.coords...) + dPdx(q.coords...)

            forcing_function = (x) -> dQdy(x...) + dPdx(x...)

            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_function)

            domain_quads = []
            boundary_quads = []
            for (i, mesh) in enumerate(meshes)
                push!(domain_quads, AdaptiveMeshSolver.getDomainQuadrature(mesh, 4))
                if i != 1
                    push!(boundary_quads, AdaptiveMeshSolver.getBoundaryQuadrature(mesh, 6))
                end
            end
            domain_integral = 0.0
            for quad in domain_quads
                domain_integral += Inti.integrate(domain_fn, quad)
            end
            boundary_integral = 0.0
            for quad in boundary_quads
                boundary_integral += Inti.integrate(boundary_fn, quad)
            end
            calculated_u_val = abs(domain_integral - boundary_integral)
            println(domain_integral)
            println(boundary_integral)
            println(calculated_u_val)
            @test 0.0 ≈ calculated_u_val atol=1e-14
        end
    end

    function test_quadtree_greens_theorem_linear_y()
        @testset "Greens Theorem Area" begin
            P = (x, y) -> 0
            Q = (x, y) -> y
            dPdx = (x, y) -> 0
            dQdy = (x, y) -> 1

            boundary_fn = (q) -> P(q.coords...)*q.normal[1] + Q(q.coords...)*q.normal[2]
            domain_fn = (q) -> dQdy(q.coords...) + dPdx(q.coords...)

            forcing_function = (x) -> dQdy(x...) + dPdx(x...)

            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_function)

            domain_quads = []
            boundary_quads = []
            for (i, mesh) in enumerate(meshes)
                push!(domain_quads, AdaptiveMeshSolver.getDomainQuadrature(mesh, 4))
                if i != 1
                    push!(boundary_quads, AdaptiveMeshSolver.getBoundaryQuadrature(mesh, 6))
                end
            end
            domain_integral = 0.0
            for quad in domain_quads
                domain_integral += Inti.integrate(domain_fn, quad)
            end
            boundary_integral = 0.0
            for quad in boundary_quads
                boundary_integral += Inti.integrate(boundary_fn, quad)
            end
            calculated_u_val = abs(domain_integral - boundary_integral)
            println(domain_integral)
            println(boundary_integral)
            println(calculated_u_val)
            @test 0.0 ≈ calculated_u_val atol=1e-14
        end
    end

    function test_quadtree_greens_theorem_quadratic()
        @testset "Greens Theorem Area" begin
            P = (x, y) -> x^2
            Q = (x, y) -> y^2
            dPdx = (x, y) -> 2*x
            dQdy = (x, y) -> 2*y

            boundary_fn = (q) -> P(q.coords...)*q.normal[1] + Q(q.coords...)*q.normal[2]
            domain_fn = (q) -> dQdy(q.coords...) + dPdx(q.coords...)

            forcing_function = (x) -> dQdy(x...) + dPdx(x...)

            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_function)

            domain_quads = []
            boundary_quads = []
            for (i, mesh) in enumerate(meshes)
                push!(domain_quads, AdaptiveMeshSolver.getDomainQuadrature(mesh, 4))
                if i != 1
                    push!(boundary_quads, AdaptiveMeshSolver.getBoundaryQuadrature(mesh, 6))
                end
            end
            domain_integral = 0.0
            for quad in domain_quads
                domain_integral += Inti.integrate(domain_fn, quad)
            end
            boundary_integral = 0.0
            for quad in boundary_quads
                boundary_integral += Inti.integrate(boundary_fn, quad)
            end
            calculated_u_val = abs(domain_integral - boundary_integral)
            println(domain_integral)
            println(boundary_integral)
            println(calculated_u_val)
            @test 0.0 ≈ calculated_u_val atol=1e-14
        end
    end

    function test_quadtree_greens_theorem_high_order()
        @testset "Greens Theorem Area" begin
            order = 6
            P = (x, y) -> x^order
            Q = (x, y) -> y^order
            dPdx = (x, y) -> order*x^(order-1)
            dQdy = (x, y) -> order*y^(order-1)

            boundary_fn = (q) -> P(q.coords...)*q.normal[1] + Q(q.coords...)*q.normal[2]
            domain_fn = (q) -> dQdy(q.coords...) + dPdx(q.coords...)

            forcing_function = (x) -> dQdy(x...) + dPdx(x...)

            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_function)

            dom_q_order = 5 # one below bc uses derivatives
            bndry_q_order = 6 # exactly degree bc uses normal funcs

            domain_quads = []
            boundary_quads = []
            for (i, mesh) in enumerate(meshes)
                push!(domain_quads, AdaptiveMeshSolver.getDomainQuadrature(mesh, dom_q_order))
                if i != 1
                    push!(boundary_quads, AdaptiveMeshSolver.getBoundaryQuadrature(mesh, bndry_q_order))
                end
            end
            domain_integral = 0.0
            for quad in domain_quads
                domain_integral += Inti.integrate(domain_fn, quad)
            end
            boundary_integral = 0.0
            for quad in boundary_quads
                boundary_integral += Inti.integrate(boundary_fn, quad)
            end
            calculated_u_val = abs(domain_integral - boundary_integral)
            println(domain_integral)
            println(boundary_integral)
            println(calculated_u_val)
            @test 0.0 ≈ calculated_u_val atol=1e-14
        end
    end

    function test_quadtree_greens_theorem_sine_terms()
        @testset "Greens Theorem Area" begin
            P = (x, y) -> sin(x)*sin(y)
            Q = (x, y) -> cos(x)*cos(y)
            dPdx = (x, y) -> cos(x)*sin(y)
            dQdy = (x, y) -> -cos(x)*sin(y)

            boundary_fn = (q) -> P(q.coords...)*q.normal[1] + Q(q.coords...)*q.normal[2]
            domain_fn = (q) -> dQdy(q.coords...) + dPdx(q.coords...)

            forcing_function = (x) -> dQdy(x...) + dPdx(x...)

            parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
            meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_function)

            dom_q_order = 5 # one below bc uses derivatives
            bndry_q_order = 6 # exactly degree bc uses normal funcs

            domain_quads = []
            boundary_quads = []
            for (i, mesh) in enumerate(meshes)
                push!(domain_quads, AdaptiveMeshSolver.getDomainQuadrature(mesh, dom_q_order))
                if i != 1
                    push!(boundary_quads, AdaptiveMeshSolver.getBoundaryQuadrature(mesh, bndry_q_order))
                end
            end
            domain_integral = 0.0
            for quad in domain_quads
                domain_integral += Inti.integrate(domain_fn, quad)
            end
            boundary_integral = 0.0
            for quad in boundary_quads
                boundary_integral += Inti.integrate(boundary_fn, quad)
            end
            calculated_u_val = abs(domain_integral - boundary_integral)
            println(domain_integral)
            println(boundary_integral)
            println(calculated_u_val)
            @test 0.0 ≈ calculated_u_val atol=1e-14
        end
    end

    function test_quadtree_greens_theorem_complex_terms()
        @testset "Greens Theorem Area" begin
            for σ in [1/2^x for x in 4:6]
                x_0 = 0.8
                y_0 = 0.0
                u = (x, y) -> exp(-((x-x_0)^2 + (y-y_0)^2)/(2*σ^2))

                P = (x, y) -> -(x-x_0) / σ^2 * u(x, y)
                Q = (x, y) -> -(y-y_0) / σ^2 * u(x, y)
                Δu = (x, y) -> (((x-x_0)^2 + (y-y_0)^2)/σ^4 - 2 / σ^2) * u(x, y)

                boundary_fn = (q) -> P(q.coords...)*q.normal[1] + Q(q.coords...)*q.normal[2]
                domain_fn = (q) -> Δu(q.coords...)

                forcing_function = (x) -> Δu(x...)

                parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
                meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_function, true)
                AdaptiveMeshSolver.showMeshes(meshes)

                domain_quads = []
                boundary_quads = []
                for (i, mesh) in enumerate(meshes)
                    push!(domain_quads, AdaptiveMeshSolver.getDomainQuadrature(mesh, 12))
                    if i != 1
                        push!(boundary_quads, AdaptiveMeshSolver.getBoundaryQuadrature(mesh, 6)) # order 6 is fine
                    end
                end
                domain_integral = 0.0
                for quad in domain_quads
                    domain_integral += Inti.integrate(domain_fn, quad)
                end
                boundary_integral = 0.0
                for quad in boundary_quads
                    boundary_integral += Inti.integrate(boundary_fn, quad)
                end
                calculated_u_val = abs(domain_integral - boundary_integral)
                # println(domain_integral)
                # println(boundary_integral)
                println(calculated_u_val)
                # @test 0.0 ≈ calculated_u_val atol=1e-14       
            end
        end
    end
end

@testitem "Greens_Theorem: Greens Theorem Constants" setup=[Greens_Theorem] begin
    Greens_Theorem.test_quadtree_greens_theorem_constants()
end

@testitem "Greens_Theorem: Greens Theorem Linear" setup=[Greens_Theorem] begin
    Greens_Theorem.test_quadtree_greens_theorem_linear()
end

@testitem "Greens_Theorem: Greens Theorem Linear X" setup=[Greens_Theorem] begin
    Greens_Theorem.test_quadtree_greens_theorem_linear_x()
end

@testitem "Greens_Theorem: Greens Theorem Linear Y" setup=[Greens_Theorem] begin
    Greens_Theorem.test_quadtree_greens_theorem_linear_y()
end

@testitem "Greens_Theorem: Greens Theorem Quadratic" setup=[Greens_Theorem] begin
    Greens_Theorem.test_quadtree_greens_theorem_quadratic()
end

@testitem "Greens_Theorem: Greens Theorem High Order" setup=[Greens_Theorem] begin
    Greens_Theorem.test_quadtree_greens_theorem_high_order()
end

@testitem "Greens_Theorem: Greens Theorem Sine Terms" setup=[Greens_Theorem] begin
    Greens_Theorem.test_quadtree_greens_theorem_sine_terms()
end

@testitem "Greens_Theorem: Greens Theorem Complex Terms" setup=[Greens_Theorem] begin
    Greens_Theorem.test_quadtree_greens_theorem_complex_terms()
end