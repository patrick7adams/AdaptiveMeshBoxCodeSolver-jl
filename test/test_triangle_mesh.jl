include("../src/main.jl")
using Test

function test_simple_triangle_mesh_area()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 5)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> 1

    result = Inti.integrate(test_function, dom_quad)

    println(abs(result - pi))
    @test result ≈ pi atol=1e-15
end

function test_simple_triangle_mesh_zero()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> 0

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_triangle_mesh_negative_area()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> -1

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ -pi atol=1e-6
end

function test_simple_triangle_mesh_linear_x()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> q.coords[1]

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_triangle_mesh_linear_y()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> q.coords[2]

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_triangle_mesh_quadratic_x()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> q.coords[1]^2

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ pi/4 atol=1e-6
end

function test_simple_triangle_mesh_quadratic_y()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> q.coords[2]^2

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ pi/4 atol=1e-6
end

function test_simple_triangle_mesh_quadratic_xy()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> q.coords[1]*q.coords[2]

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_triangle_mesh_high_order_x()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> q.coords[1]^100

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0.00490268 atol=1e-6
end

function test_simple_triangle_mesh_high_order_y()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> q.coords[2]^100

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0.00490268 atol=1e-6
end

function test_simple_triangle_mesh_high_order_xy()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> q.coords[1]^51*q.coords[2]^51

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_triangle_mesh_sin_x()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> sin(q.coords[1])

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_triangle_mesh_sin_y()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> sin(q.coords[2])

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_triangle_mesh_sin_xy()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> sin(q.coords[1]*q.coords[2])

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_triangle_mesh_sin_x_cos_y()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> sin(q.coords[1])*cos(q.coords[2])

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_triangle_mesh_exp_xy()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> exp(q.coords[1]+q.coords[2])

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 3.9952371 atol=1e-6
end

function test_simple_triangle_mesh_spike()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> exp(-100 * (q.coords[1]^2+q.coords[2]^2))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0.0314159 atol=1e-6
end