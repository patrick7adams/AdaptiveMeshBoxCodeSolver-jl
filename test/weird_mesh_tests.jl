include("../src/main.jl")

using Test

function test_simple_meshsizes()
    meshsizes = [(0.5)^i for i in 2:6]
    for i in 1:length(meshsizes)
        x_test = (1, 1) # point just outside of the region
        r0 = (-0.5, 0.0)
        u = (x) -> -2*pi^2 * sin(pi*x[1]) * sin(pi*x[2]) + exp(-600*((x[1]- r0[1])^2 + (x[2] - r0[2])^2)) # forcing function (is it sufficiently smooth?)
        # check these eqns lol
        partial_x_u = (x) -> -2*pi^3 * cos(pi*x[1]) * sin(pi*x[2]) - 1200*(x[1]-r0[1])*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
        partial_y_u = (x) -> -2*pi^3 * sin(pi*x[1]) * cos(pi*x[2]) - 1200*(x[2]-r0[2])*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
        laplacian_u = (x) -> 4*pi^4 * sin(pi*x[1]) * sin(pi*x[2]) + (1440000*((x[1]-r0[1])^2 + (x[2]-r0[2])^2)-2400)*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
        greens_fn = (r, x) -> 1/(2pi) * log(distance(x, r))
        partial_x_greens_fn = (r, x) -> 1/(2pi) * ((x[1]-r[1]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))
        partial_y_greens_fn = (r, x) -> 1/(2pi) * ((x[2]-r[2]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))

        normal_derivative_u = (x, n) -> partial_x_u(x) * n[1] + partial_y_u(x) * n[2]
        normal_derivative_greens_fn = (r, x, n) -> partial_x_greens_fn(r, x) * n[1] + partial_y_greens_fn(r, x) * n[2]

        boundary_function = (q) -> (greens_fn(x_test, q.coords) * normal_derivative_u(q.coords, q.normal) - 
                                    u(q.coords) * normal_derivative_greens_fn(x_test, q.coords, q.normal))

        domain_function = (q) -> (greens_fn(x_test, q.coords) * laplacian_u(q.coords))

        expected_u_val = u(x_test)

        gmsh.initialize()
        parametrizations::Vector{Function} = [(x) -> (cos(x*2*pi), sin(x*2*pi))]
        quadtree_mesh = createQuadtreeMesh(parametrizations, u, meshsizes[i], meshsizes[i])

        showMesh(quadtree_mesh)

        dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
        boundary = Inti.boundary(dom)

        dom_mesh = Inti.view(quadtree_mesh, dom)
        boundary_mesh = Inti.view(quadtree_mesh, boundary)

        dom_quad = Inti.Quadrature(dom_mesh; qorder = 17)
        boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 12)

        boundary_integral = Inti.integrate(boundary_function, boundary_quad)
        domain_integral = Inti.integrate(domain_function, dom_quad)
        calculated_u_val = domain_integral - boundary_integral

        @test 0.0 ≈ calculated_u_val atol=1e-13
    end
end

function test_simple_hole_mesh()
    meshsize = 0.1
    x_test = (1, 1) # point just outside of the region
    r0 = (-0.5, 0.0)
    u = (x) -> -2*pi^2 * sin(pi*x[1]) * sin(pi*x[2]) + exp(-600*((x[1]- r0[1])^2 + (x[2] - r0[2])^2)) # forcing function (is it sufficiently smooth?)
    # check these eqns lol
    partial_x_u = (x) -> -2*pi^3 * cos(pi*x[1]) * sin(pi*x[2]) - 1200*(x[1]-r0[1])*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
    partial_y_u = (x) -> -2*pi^3 * sin(pi*x[1]) * cos(pi*x[2]) - 1200*(x[2]-r0[2])*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
    laplacian_u = (x) -> 4*pi^4 * sin(pi*x[1]) * sin(pi*x[2]) + (1440000*((x[1]-r0[1])^2 + (x[2]-r0[2])^2)-2400)*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
    greens_fn = (r, x) -> 1/(2pi) * log(distance(x, r))
    partial_x_greens_fn = (r, x) -> 1/(2pi) * ((x[1]-r[1]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))
    partial_y_greens_fn = (r, x) -> 1/(2pi) * ((x[2]-r[2]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))

    normal_derivative_u = (x, n) -> partial_x_u(x) * n[1] + partial_y_u(x) * n[2]
    normal_derivative_greens_fn = (r, x, n) -> partial_x_greens_fn(r, x) * n[1] + partial_y_greens_fn(r, x) * n[2]

    boundary_function = (q) -> (greens_fn(x_test, q.coords) * normal_derivative_u(q.coords, q.normal) - 
                                u(q.coords) * normal_derivative_greens_fn(x_test, q.coords, q.normal))

    domain_function = (q) -> (greens_fn(x_test, q.coords) * laplacian_u(q.coords))

    expected_u_val = u(x_test)

    gmsh.initialize()
    parametrizations::Vector{Function} = [(x) -> (cos(x*2*pi), sin(x*2*pi)), (x) -> (0.1*cos(x*2*pi)+0.4, 0.1*sin(x*2*pi))]
    quadtree_mesh = createQuadtreeMesh(parametrizations, u)

    showMesh(quadtree_mesh)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 17)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 12)

    boundary_integral = Inti.integrate(boundary_function, boundary_quad)
    # domain_integral = Inti.integrate(domain_function, dom_quad)
    domain_integral = integrate(domain_function, dom_quad, quadtree_mesh, 17, r0)
    calculated_u_val = domain_integral - boundary_integral

    @test 0.0 ≈ calculated_u_val atol=1e-13
end

function test_square_hole_mesh()
    x_test = (1, 1) # point just outside of the region
    r0 = (-0.5, 0.0)
    u = (x) -> -2*pi^2 * sin(pi*x[1]) * sin(pi*x[2]) + exp(-600*((x[1]- r0[1])^2 + (x[2] - r0[2])^2)) # forcing function (is it sufficiently smooth?)
    # check these eqns lol
    partial_x_u = (x) -> -2*pi^3 * cos(pi*x[1]) * sin(pi*x[2]) - 1200*(x[1]-r0[1])*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
    partial_y_u = (x) -> -2*pi^3 * sin(pi*x[1]) * cos(pi*x[2]) - 1200*(x[2]-r0[2])*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
    laplacian_u = (x) -> 4*pi^4 * sin(pi*x[1]) * sin(pi*x[2]) + (1440000*((x[1]-r0[1])^2 + (x[2]-r0[2])^2)-2400)*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
    greens_fn = (r, x) -> 1/(2pi) * log(distance(x, r))
    partial_x_greens_fn = (r, x) -> 1/(2pi) * ((x[1]-r[1]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))
    partial_y_greens_fn = (r, x) -> 1/(2pi) * ((x[2]-r[2]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))

    normal_derivative_u = (x, n) -> partial_x_u(x) * n[1] + partial_y_u(x) * n[2]
    normal_derivative_greens_fn = (r, x, n) -> partial_x_greens_fn(r, x) * n[1] + partial_y_greens_fn(r, x) * n[2]

    boundary_function = (q) -> (greens_fn(x_test, q.coords) * normal_derivative_u(q.coords, q.normal) - 
                                u(q.coords) * normal_derivative_greens_fn(x_test, q.coords, q.normal))

    domain_function = (q) -> (greens_fn(x_test, q.coords) * laplacian_u(q.coords))

    expected_u_val = u(x_test)

    gmsh.initialize()

    function polygonal_hole(t)
        if t < 0.25
            return (0.275, -0.125+t)
        elseif t < 0.5 && t >= 0.25
            return (0.025+t, 0.125)
        elseif t < 0.75 && t >= 0.5
            return (0.525, 0.625-t)
        else
            return (1.275-t, -0.125)
        end
    end
    parametrizations::Vector{Function} = [
        (x) -> (cos(x*2*pi), sin(x*2*pi)), 
        (x) -> polygonal_hole(x)]
    quadtree_mesh = createQuadtreeMesh(parametrizations, u)

    showMesh(quadtree_mesh)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 17)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 12)

    boundary_integral = Inti.integrate(boundary_function, boundary_quad)
    domain_integral = Inti.integrate(domain_function, dom_quad)
    calculated_u_val = domain_integral - boundary_integral

    @test 0.0 ≈ calculated_u_val atol=1e-13
end

function test_ellipse_mesh()
    x_test = (1, 1) # point just outside of the region
    r0 = (-0.5, 0.0)
    u = (x) -> -2*pi^2 * sin(pi*x[1]) * sin(pi*x[2]) + exp(-600*((x[1]- r0[1])^2 + (x[2] - r0[2])^2)) # forcing function (is it sufficiently smooth?)
    # check these eqns lol
    partial_x_u = (x) -> -2*pi^3 * cos(pi*x[1]) * sin(pi*x[2]) - 1200*(x[1]-r0[1])*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
    partial_y_u = (x) -> -2*pi^3 * sin(pi*x[1]) * cos(pi*x[2]) - 1200*(x[2]-r0[2])*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
    laplacian_u = (x) -> 4*pi^4 * sin(pi*x[1]) * sin(pi*x[2]) + (1440000*((x[1]-r0[1])^2 + (x[2]-r0[2])^2)-2400)*exp(-600*((x[1]-r0[1])^2 + (x[2]-r0[2])^2))
    greens_fn = (r, x) -> 1/(2pi) * log(distance(x, r))
    partial_x_greens_fn = (r, x) -> 1/(2pi) * ((x[1]-r[1]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))
    partial_y_greens_fn = (r, x) -> 1/(2pi) * ((x[2]-r[2]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))

    normal_derivative_u = (x, n) -> partial_x_u(x) * n[1] + partial_y_u(x) * n[2]
    normal_derivative_greens_fn = (r, x, n) -> partial_x_greens_fn(r, x) * n[1] + partial_y_greens_fn(r, x) * n[2]

    boundary_function = (q) -> (greens_fn(x_test, q.coords) * normal_derivative_u(q.coords, q.normal) - 
                                u(q.coords) * normal_derivative_greens_fn(x_test, q.coords, q.normal))

    domain_function = (q) -> (greens_fn(x_test, q.coords) * laplacian_u(q.coords))

    expected_u_val = u(x_test)

    gmsh.initialize()

    parametrizations::Vector{Function} = [(x) -> (2.2*cos(x*2*pi), sin(x*2*pi))]
    quadtree_mesh = createQuadtreeMesh(parametrizations, u)

    showMesh(quadtree_mesh)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 17)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 12)

    boundary_integral = Inti.integrate(boundary_function, boundary_quad)
    domain_integral = Inti.integrate(domain_function, dom_quad)
    calculated_u_val = domain_integral - boundary_integral

    @test 0.0 ≈ calculated_u_val atol=1e-13
end

function test_ellipse_mesh_area()
    forcing_function = (x) -> 1

    a = 2
    b = 1
    gmsh.initialize()
    parametrizations::Vector{Function} = [(x) -> (a*cos(x*2*pi), b*sin(x*2*pi))]
    quadtree_mesh = createQuadtreeMesh(parametrizations, forcing_function)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    dom_mesh = Inti.view(quadtree_mesh, dom)
    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)

    showMesh(quadtree_mesh)
    
    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)
    
    @test result ≈ a*b*pi atol=1e-15
end