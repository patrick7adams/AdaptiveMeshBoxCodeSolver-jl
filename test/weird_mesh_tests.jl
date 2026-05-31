using TestItems

@testmodule Weird_Meshes begin
    using Test
    using AdaptiveMeshSolver

    function test_simple_hole_mesh()
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

        parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))], [(x) -> (0.1*cos(-x*2*pi)+0.4, 0.1*sin(-x*2*pi))]]
        meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_func)

        AdaptiveMeshSolver.showMeshes(meshes)
        domain_integral, boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
            meshes, r0; 
            dom_func = domain_function,
            bndry_func = boundary_function, 
            dom_order = 17,
            bndry_order = 12
        )

        # println(domain_integral)
        # println(boundary_integral)

        calculated_u_val = domain_integral - boundary_integral

        @test 0.0 ≈ calculated_u_val atol=1e-12
    end

    function test_square_hole_mesh()
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

        parametrizations::Vector{Vector{Function}} = [
            [(x) -> (cos(x*2*pi), sin(x*2*pi))], 
            [(t) -> (0.275, -0.125+0.25*t), (t) -> (0.275+0.25*t, 0.125), (t) -> (0.525, 0.125-0.25*t), (t) -> (0.525-0.25*t, -0.125)]]
        meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_func)

        AdaptiveMeshSolver.showMeshes(meshes)
        domain_integral, boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
            meshes, r0; 
            dom_func = domain_function,
            bndry_func = boundary_function, 
            dom_order = 17,
            bndry_order = 12
        )

        # println(domain_integral)
        # println(boundary_integral)

        calculated_u_val = domain_integral - boundary_integral

        @test 0.0 ≈ calculated_u_val atol=1e-12
    end

    function test_hard_ellipse_mesh()
        x_test = (2, 2) # point just outside of the region
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

        parametrizations::Vector{Vector{Function}} = [[(x) -> (3.0*cos(x*2*pi), 1.0*sin(x*2*pi))]]
        meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_func)
        # showMeshes(meshes)
        domain_integral, boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
            meshes, r0; 
            dom_func = domain_function,
            bndry_func = boundary_function, 
            dom_order = 17,
            bndry_order = 12
        )

        # println(domain_integral)
        # println(boundary_integral)

        calculated_u_val = domain_integral - boundary_integral

        @test 0.0 ≈ calculated_u_val atol=1e-12
    end

    function test_mild_ellipse_mesh()
        x_test = (2, 2) # point just outside of the region
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

        parametrizations::Vector{Vector{Function}} = [[(x) -> (1.5*cos(x*2*pi), 1.0*sin(x*2*pi))]]
        meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_func)
        # showMeshes(meshes)
        domain_integral, boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
            meshes, r0; 
            dom_func = domain_function,
            bndry_func = boundary_function, 
            dom_order = 17,
            bndry_order = 12
        )

        # println(domain_integral)
        # println(boundary_integral)

        calculated_u_val = domain_integral - boundary_integral

        @test 0.0 ≈ calculated_u_val atol=1e-12
    end

    function test_med_ellipse_mesh()
        x_test = (2, 2) # point just outside of the region
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

        parametrizations::Vector{Vector{Function}} = [[(x) -> (2.25*cos(x*2*pi), 1.0*sin(x*2*pi))]]
        meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_func)
        # showMeshes(meshes)
        domain_integral, boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
            meshes, r0; 
            dom_func = domain_function,
            bndry_func = boundary_function, 
            dom_order = 17,
            bndry_order = 12
        )

        # println(domain_integral)
        # println(boundary_integral)

        calculated_u_val = domain_integral - boundary_integral

        @test 0.0 ≈ calculated_u_val atol=1e-12
    end

    function test_ellipse_with_simple_hole_mesh()
        x_test = (2, 2) # point just outside of the region
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

        parametrizations::Vector{Vector{Function}} = [
            [(x) -> (3.0*cos(x*2*pi), 1.0*sin(x*2*pi))],
            [(x) -> (0.1*cos(-x*2*pi)+0.4, 0.1*sin(-x*2*pi))]]
        meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_func)
        # showMeshes(meshes)
        domain_integral, boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
            meshes, r0; 
            dom_func = domain_function,
            bndry_func = boundary_function, 
            dom_order = 17,
            bndry_order = 12
        )

        # println(domain_integral)
        # println(boundary_integral)

        calculated_u_val = domain_integral - boundary_integral

        @test 0.0 ≈ calculated_u_val atol=1e-12
    end

    function test_goofy_hole_mesh()
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

        parametrizations::Vector{Vector{Function}} = [
            [(x) -> (cos(x*2*pi), sin(x*2*pi))], 
            [(t) -> (0.25*t, 0.0), 
            (t) -> (0.25, 0.25*t), 
            (t) -> (0.25 + 0.15*t, 0.25-0.25*t), 
            (t) -> (0.4-0.2*t, -0.125*t), 
            (t) -> (0.2-0.2*t, -0.125+0.125*t)]]
        meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_func)

        # AdaptiveMeshSolver.showMeshes(meshes)
        domain_integral, boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
            meshes, r0; 
            dom_func = domain_function,
            bndry_func = boundary_function, 
            dom_order = 17,
            bndry_order = 12
        )

        # println(domain_integral)
        # println(boundary_integral)

        calculated_u_val = domain_integral - boundary_integral

        @test 0.0 ≈ calculated_u_val atol=1e-12
    end

    function test_two_hole_mesh()
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

        parametrizations::Vector{Vector{Function}} = [
            [(x) -> (cos(x*2*pi), sin(x*2*pi))], 
            [(x) -> (0.1*cos(-x*2*pi)+0.4, 0.1*sin(-x*2*pi)+0.25)],
            [(x) -> (0.1*cos(-x*2*pi)+0.4, 0.1*sin(-x*2*pi)-0.25)]]
        meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_func)

        # AdaptiveMeshSolver.showMeshes(meshes)
        domain_integral, boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
            meshes, r0; 
            dom_func = domain_function,
            bndry_func = boundary_function, 
            dom_order = 17,
            bndry_order = 12
        )

        # println(domain_integral)
        # println(boundary_integral)

        calculated_u_val = domain_integral - boundary_integral

        @test 0.0 ≈ calculated_u_val atol=1e-12
    end

    function test_square_boundary_mesh()
        x_test = (0, 1.5) # point just outside of the region
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

        parametrizations::Vector{Vector{Function}} = [
            [(t) -> (1, -1+2*t), (t) -> (1-2*t, 1), (t) -> (-1, 1-2*t), (t) -> (-1+2*t, -1)]
        ]
        meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_func)

        # AdaptiveMeshSolver.showMeshes(meshes)
        domain_integral, boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
            meshes, r0; 
            dom_func = domain_function,
            bndry_func = boundary_function, 
            dom_order = 17,
            bndry_order = 12
        )

        # println(domain_integral)
        # println(boundary_integral)

        calculated_u_val = domain_integral - boundary_integral

        @test 0.0 ≈ calculated_u_val atol=1e-12
    end

    function test_goofy_boundary_mesh()
        x_test = (0, 1.5) # point just outside of the region
        r0 = (0.0, 0.5)
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

        parametrizations::Vector{Vector{Function}} = [
            [(t) -> (1, -1+2*t), (t) -> (1-2*t, 1), (t) -> (-1, 1-2*t), (t) -> (-1+0.5*t, -1), (t) -> (-0.5, -1+t), (t) -> (-0.5+t, 0), (t) -> (0.5, -t), (t) -> (0.5+0.5*t, -1)]
        ]
        meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_func)

        # AdaptiveMeshSolver.showMeshes(meshes)
        domain_integral, boundary_integral = AdaptiveMeshSolver.calculateIntegrals(
            meshes, r0; 
            dom_func = domain_function,
            bndry_func = boundary_function, 
            dom_order = 17,
            bndry_order = 12
        )

        # println(domain_integral)
        # println(boundary_integral)

        calculated_u_val = domain_integral - boundary_integral

        @test 0.0 ≈ calculated_u_val atol=1e-12
    end
end

@testitem "Weird_Meshes: Simple Hole" setup=[Weird_Meshes] begin
    Weird_Meshes.test_simple_hole_mesh()
end

@testitem "Weird_Meshes: Square Hole" setup=[Weird_Meshes] begin
    Weird_Meshes.test_square_hole_mesh()
end

@testitem "Weird_Meshes: Mild Ellipse" setup=[Weird_Meshes] begin
    Weird_Meshes.test_mild_ellipse_mesh()
end

@testitem "Weird_Meshes: Medium Ellipse" setup=[Weird_Meshes] begin
    Weird_Meshes.test_med_ellipse_mesh()
end

@testitem "Weird_Meshes: Hard Ellipse" setup=[Weird_Meshes] begin
    Weird_Meshes.test_hard_ellipse_mesh()
end

@testitem "Weird_Meshes: Simple Hole Ellipse" setup=[Weird_Meshes] begin
    Weird_Meshes.test_ellipse_with_simple_hole_mesh()
end

@testitem "Weird_Meshes: Goofy Hole" setup=[Weird_Meshes] begin
    Weird_Meshes.test_goofy_hole_mesh()
end

@testitem "Weird_Meshes: Two Hole" setup=[Weird_Meshes] begin
    Weird_Meshes.test_two_hole_mesh()
end

@testitem "Weird_Meshes: Square Boundary" setup=[Weird_Meshes] begin
    Weird_Meshes.test_square_boundary_mesh()
end

@testitem "Weird_Meshes: Goofy Boundary" setup=[Weird_Meshes] begin
    Weird_Meshes.test_goofy_boundary_mesh()
end