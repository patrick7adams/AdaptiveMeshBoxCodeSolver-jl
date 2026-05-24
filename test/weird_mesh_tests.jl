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
    quadtree_mesh = createQuadtreeMesh(parametrizations, u, meshsize, meshsize)

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