include("../src/main.jl")

using Test

function test_simple_quadtree_greens_theorem()
    meshsize = 0.25
    mesh = createMesh(meshsize, meshsize)
    showMesh(mesh)

    P = (x, y) -> -y
    Q = (x, y) -> x
    dPdy = (x, y) -> -1
    dQdx = (x, y) -> 1

    boundary_fn = (q) -> P(q.coords...)*q.normal[2]*-1 + Q(q.coords...)*q.normal[1]
    domain_fn = (q) -> dQdx(q.coords...) - dPdy(q.coords...)

    forcing_function = (x) -> dQdx(x...) - dPdy(x...)

    gmsh.initialize()
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

    # println(Inti.entities(quadtree_mesh))
    showMesh(quadtree_mesh)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)
    # println(dom)
    # println(boundary)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)
    # println(dom_mesh)
    # println(boundary_mesh)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    computed_boundary_int = Inti.integrate(boundary_fn, boundary_quad)
    computed_domain_int = Inti.integrate(domain_fn, dom_quad)

    println(0.5*computed_boundary_int - pi)
    println(0.5*computed_domain_int - pi)
    println("Error in values: ", abs(computed_boundary_int - computed_domain_int))
end

function test_simple_quadtree_greens_third_identity_linear_x()
    meshsize = 0.05
    mesh = createMesh(meshsize, meshsize)
    # showMesh(mesh)

    x_test = (1, 1) # point just outside of the region
    u = (x) -> x[1]
    # check these eqns lol
    partial_x_u = (x) -> 1.0
    partial_y_u = (x) -> 0.0
    laplacian_u = (x) -> 0.0
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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, u, meshsize, meshsize, parametrizations)

    # showMesh(quadtree_mesh)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    boundary_integral = Inti.integrate(boundary_function, boundary_quad)
    domain_integral = Inti.integrate(domain_function, dom_quad)
    calculated_u_val = domain_integral - boundary_integral

    println("Boundary integral = ", Inti.integrate(boundary_function, boundary_quad))
    println("Domain Integral = ", Inti.integrate(domain_function, dom_quad))
    println("Difference = ", calculated_u_val)
    @test 0.0 ≈ calculated_u_val atol=1e-15
end

function test_simple_quadtree_greens_third_identity_linear_y()
    meshsize = 0.05
    mesh = createMesh(meshsize, meshsize)
    # showMesh(mesh)

    x_test = (1, 1) # point just outside of the region
    u = (x) -> x[2]
    # check these eqns lol
    partial_x_u = (x) -> 0.0
    partial_y_u = (x) -> 1.0
    laplacian_u = (x) -> 0.0
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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, u, meshsize, meshsize, parametrizations)

    # showMesh(quadtree_mesh)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    boundary_integral = Inti.integrate(boundary_function, boundary_quad)
    domain_integral = Inti.integrate(domain_function, dom_quad)
    calculated_u_val = domain_integral - boundary_integral

    println("Boundary integral = ", Inti.integrate(boundary_function, boundary_quad))
    println("Domain Integral = ", Inti.integrate(domain_function, dom_quad))
    println("Difference = ", calculated_u_val)
    @test 0.0 ≈ calculated_u_val atol=1e-15
end

function test_simple_quadtree_greens_third_identity_simple_forcing()
    meshsize = 0.05
    mesh = createMesh(meshsize, meshsize)
    # showMesh(mesh)

    x_test = (1, 1) # point just outside of the region
    u = (x) -> 1.0
    # check these eqns lol
    partial_x_u = (x) -> 0.0
    partial_y_u = (x) -> 0.0
    laplacian_u = (x) -> 0.0
    greens_fn = (r, x) -> 1/(2pi) * log(distance(x, r))
    partial_x_greens_fn = (r, x) -> 1/(2pi) * ((x[1]-r[1]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))
    partial_y_greens_fn = (r, x) -> 1/(2pi) * ((x[2]-r[2]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))

    normal_derivative_u = (x, n) -> partial_x_u(x) * n[1] + partial_y_u(x) * n[2]
    normal_derivative_greens_fn = (r, x, n) -> partial_x_greens_fn(r, x) * n[1] + partial_y_greens_fn(r, x) * n[2]

    boundary_function = (q) -> (greens_fn(x_test, q.coords) * normal_derivative_u(q.coords, q.normal) - 
                                u(q.coords) * normal_derivative_greens_fn(x_test, q.coords, q.normal))

    domain_function = (q) -> (greens_fn(x_test, q.coords) * laplacian_u(q.coords))

    gmsh.initialize()
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, u, meshsize, meshsize, parametrizations)

    # showMesh(quadtree_mesh)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    boundary_integral = Inti.integrate(boundary_function, boundary_quad)
    domain_integral = Inti.integrate(domain_function, dom_quad)
    calculated_u_val = domain_integral - boundary_integral

    println("Boundary integral = ", Inti.integrate(boundary_function, boundary_quad))
    println("Domain Integral = ", Inti.integrate(domain_function, dom_quad))
    println("Difference = ", calculated_u_val)
    @test 0.0 ≈ calculated_u_val atol=1e-15
end

function test_simple_quadtree_greens_third_identity_quadratic()
    meshsize = 0.05
    mesh = createMesh(meshsize, meshsize)
    # showMesh(mesh)

    x_test = (1, 1) # point just outside of the region
    u = (x) -> x[1]^2 + x[2]^2
    # check these eqns lol
    partial_x_u = (x) -> 2*x[1]
    partial_y_u = (x) -> 2*x[2]
    laplacian_u = (x) -> 4.0
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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, u, meshsize, meshsize, parametrizations)

    # showMesh(quadtree_mesh)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 17)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 12)

    boundary_integral = Inti.integrate(boundary_function, boundary_quad)
    domain_integral = Inti.integrate(domain_function, dom_quad)
    calculated_u_val = domain_integral - boundary_integral

    println("Boundary integral = ", Inti.integrate(boundary_function, boundary_quad))
    println("Domain Integral = ", Inti.integrate(domain_function, dom_quad))
    println("Difference = ", calculated_u_val)
    @test 0.0 ≈ calculated_u_val atol=1e-13
end

function test_simple_quadtree_greens_third_identity_sin_term()
    meshsize = 0.05
    mesh = createMesh(meshsize, meshsize)
    # showMesh(mesh)

    for k in [0.25, 0.5, 1.0, 2.0, 4.0]
    # for k in [4.0]
        x_test = (1, 1) # point just outside of the region
        u = (x) -> sin(k*pi*x[1])*sin(k*pi*x[2])
        # check these eqns lol
        partial_x_u = (x) -> k*pi*cos(k*pi*x[1])*sin(k*pi*x[2])
        partial_y_u = (x) -> k*pi*sin(k*pi*x[1])*cos(k*pi*x[2])
        laplacian_u = (x) -> -2*k^2*pi^2*sin(k*pi*x[1])*sin(k*pi*x[2])
        greens_fn = (r, x) -> 1/(2pi) * log(distance(x, r))
        partial_x_greens_fn = (r, x) -> 1/(2pi) * ((x[1]-r[1]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))
        partial_y_greens_fn = (r, x) -> 1/(2pi) * ((x[2]-r[2]) / ((x[1]-r[1])^2 + (x[2]-r[2])^2))

        normal_derivative_u = (x, n) -> partial_x_u(x) * n[1] + partial_y_u(x) * n[2]
        normal_derivative_greens_fn = (r, x, n) -> partial_x_greens_fn(r, x) * n[1] + partial_y_greens_fn(r, x) * n[2]

        boundary_fn1 = (q) -> greens_fn(x_test, q.coords) * normal_derivative_u(q.coords, q.normal)
        boundary_fn2 = (q) -> u(q.coords) * normal_derivative_greens_fn(x_test, q.coords, q.normal)
        # boundary_function = (q) -> (greens_fn(x_test, q.coords) * normal_derivative_u(q.coords, q.normal) - 
        #                             u(q.coords) * normal_derivative_greens_fn(x_test, q.coords, q.normal))
        boundary_function = (q) -> boundary_fn1(q) - boundary_fn2(q)

        domain_function = (q) -> (greens_fn(x_test, q.coords) * laplacian_u(q.coords))

        # expected_u_val = u(x_test)
        forcing_func = (x) -> greens_fn(x_test, x) * laplacian_u(x)

        gmsh.initialize()
        parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
        quadtree_mesh = createQuadtreeMesh(mesh, forcing_func, meshsize, meshsize, parametrizations)

        # showMesh(quadtree_mesh)

        dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
        boundary = get_boundary(quadtree_mesh)

        dom_mesh = Inti.view(quadtree_mesh, dom)
        boundary_mesh = Inti.view(quadtree_mesh, boundary)

        dom_quad = Inti.Quadrature(dom_mesh; qorder = 17)
        boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 12)

        triangle_dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
        triangle_dom_mesh = Inti.view(mesh, triangle_dom)
        triangle_quad = Inti.Quadrature(triangle_dom_mesh; qorder = 17) # 17 is highest order

        triangle_boundary = Inti.boundary(triangle_dom)
        triangle_boundary_mesh = Inti.view(mesh, triangle_boundary)
        triangle_boundary_quad = Inti.Quadrature(triangle_boundary_mesh; qorder = 12)

        boundary_integral = Inti.integrate(boundary_function, boundary_quad)
        domain_integral = Inti.integrate(domain_function, dom_quad)
        calculated_u_val = domain_integral - boundary_integral

        exact_triangle_domain_integral = Inti.integrate(domain_function, triangle_quad)
        exact_triangle_boundary_integral = Inti.integrate(boundary_function, triangle_boundary_quad)

        domain_boundary_error = abs(calculated_u_val)
        quadtree_triangle_error = abs(domain_integral - exact_triangle_domain_integral)
        # this is the current dominant error, the error between the integral taken over the triangle mesh and
        # the error on the proper quadtree domain. 
        println(abs(exact_triangle_boundary_integral - boundary_integral))
        println(abs(boundary_integral - domain_integral))
        println(abs(exact_triangle_boundary_integral - exact_triangle_domain_integral))
        println(abs(exact_triangle_domain_integral - domain_integral))
        println("---------")
        println("k = ", k)
        println("Boundary integral = ", boundary_integral)
        println("Domain integral = ", domain_integral)
        println("Triangle Domain integral = ", exact_triangle_domain_integral)
        println("Domain boundary error = ", domain_boundary_error)
        println("Quadtree vs. Triangle mesh error = ", quadtree_triangle_error)
        @test 0.0 ≈ calculated_u_val atol=1e-13
    # @test expected_u_val ≈ calculated_u_val atol=1e-15
    end
end

function test_simple_quadtree_greens_third_identity_complex_forcing()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)
    # showMesh(mesh)

    x_test = (1, 1) # point just outside of the region
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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, u, meshsize, meshsize, parametrizations)

    showMesh(quadtree_mesh)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 17)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 12)

    boundary_integral = Inti.integrate(boundary_function, boundary_quad)
    domain_integral = Inti.integrate(domain_function, dom_quad)
    calculated_u_val = domain_integral - boundary_integral

    println("Boundary integral = ", Inti.integrate(boundary_function, boundary_quad))
    println("Domain Integral = ", Inti.integrate(domain_function, dom_quad))
    println("Difference = ", calculated_u_val)
    @test 0.0 ≈ calculated_u_val atol=1e-13
end