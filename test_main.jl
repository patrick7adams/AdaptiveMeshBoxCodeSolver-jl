include("mesh_tools.jl")
include("main.jl")

using Test

function test_simple_triangle_mesh()
    @testset "Simple Triangle Mesh" begin
        test_simple_triangle_mesh_area()
        test_simple_triangle_mesh_zero()
        test_simple_triangle_mesh_negative_area()
        test_simple_triangle_mesh_linear_x()
        test_simple_triangle_mesh_linear_y()
        test_simple_triangle_mesh_quadratic_x()
        test_simple_triangle_mesh_quadratic_y()
        test_simple_triangle_mesh_quadratic_xy()
        test_simple_triangle_mesh_high_order_x()
        test_simple_triangle_mesh_high_order_y()
        test_simple_triangle_mesh_high_order_xy()
        test_simple_triangle_mesh_sin_x()
        test_simple_triangle_mesh_sin_y()
        test_simple_triangle_mesh_sin_xy()
        test_simple_triangle_mesh_sin_x_cos_y()
        test_simple_triangle_mesh_exp_xy()
        test_simple_triangle_mesh_spike()
    end;
end

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

function test_simple_quadtree_mesh()
    @testset "Simple Quadtree Mesh" begin
        test_simple_quadtree_mesh_area()
        test_simple_quadtree_mesh_zero()
        test_simple_quadtree_mesh_negative_area()
        test_simple_quadtree_mesh_linear_x()
        test_simple_quadtree_mesh_linear_y()
        test_simple_quadtree_mesh_quadratic_x()
        test_simple_quadtree_mesh_quadratic_y()
        test_simple_quadtree_mesh_quadratic_xy()
        test_simple_quadtree_mesh_high_order_x()
        test_simple_quadtree_mesh_high_order_y()
        test_simple_quadtree_mesh_high_order_xy()
        test_simple_quadtree_mesh_sin_x()
        test_simple_quadtree_mesh_sin_y()
        test_simple_quadtree_mesh_sin_xy()
        test_simple_quadtree_mesh_sin_x_cos_y()
        test_simple_quadtree_mesh_exp_xy()
        test_simple_quadtree_mesh_spike()
    end;
end

function test_simple_quadtree_mesh_area()
    meshsize = 0.1
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> 1

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    # for element in Inti.elements(Inti.view(quadtree_mesh, dom))
    #     # println(element)
    #     println(element)
    # end

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)
    
    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    # showMeshWithQuadrature(quadtree_mesh, dom_quad)
    showMesh(quadtree_mesh, showBoundary = false)
    
    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)
    println(abs(result - pi))
    # order 1 - 0.020147640693671143
    # order 2 - 9.880807441575712e-6
    # order 3 - 1.2988655120338422e-6
    # order 4 - 1.0035255693097156e-7
    # order 5 - 1.9941442896964645e-7
    # order 6 - 2.0248774212916487e-7
    @test result ≈ pi atol=1e-15
end

function test_simple_quadtree_mesh_zero()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> 0

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_quadtree_mesh_negative_area()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> -1

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ -pi atol=1e-6
end

function test_simple_quadtree_mesh_linear_x()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> x[1]

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_quadtree_mesh_linear_y()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> x[2]

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_quadtree_mesh_quadratic_x()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> x[1]^2

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ pi/4 atol=1e-6
end

function test_simple_quadtree_mesh_quadratic_y()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> x[2]^2

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ pi/4 atol=1e-6
end

function test_simple_quadtree_mesh_quadratic_xy()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> x[1]*x[2]

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_quadtree_mesh_high_order_x()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> x[1]^100

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0.00490268 atol=1e-6
end

function test_simple_quadtree_mesh_high_order_y()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> x[2]^100

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0.00490268 atol=1e-6
end

function test_simple_quadtree_mesh_high_order_xy()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> x[1]^51 * x[2]^51

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_quadtree_mesh_sin_x()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> sin(x[1])

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_quadtree_mesh_sin_y()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> sin(x[2])

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_quadtree_mesh_sin_xy()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> sin(x[1]*x[2])

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_quadtree_mesh_sin_x_cos_y()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> sin(x[1]) * cos(x[2])

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0 atol=1e-6
end

function test_simple_quadtree_mesh_exp_xy()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> exp(x[1]+x[2])

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 3.9952371 atol=1e-6
end

function test_simple_quadtree_mesh_spike()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> exp(-100 * (x[1]^2+x[2]^2))

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]])

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0.0314159 atol=1e-6
end

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(2*pi*x), sin(2*pi*x)]])

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

function test_simple_quadtree_mesh_quadratures(meshsize)
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> 1

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x) -> [cos(x*2*pi), sin(x*2*pi)]]; return_quadtree_mesh=true)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    center_to_elements = Dict()
    for element in Inti.elements(Inti.view(quadtree_mesh, dom))
        verts = Inti.vertices(element)
        area = (verts[3][1] - verts[1][1])*(verts[3][2] - verts[1][2])
        center_to_elements[Inti.center(element)] = area
    end

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)
    
    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    area = 0
    for i in 0:Int32(length(dom_quad)/9)-1
        q_nodes = [dom_quad[9*i+k] for k in 1:9]
        center_coords = q_nodes[5].coords
        mesh_area = center_to_elements[center_coords]
        area = 0
        for node in q_nodes
            area = area + node.weight
        end
        @test area ≈ mesh_area atol=1e-15
    end
end

function test_simple_triangle_mesh_quadratures(meshsize)
    mesh = createMesh(meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = get_boundary(mesh)

    center_to_elements = Dict()
    elements = []
    for element in Inti.elements(Inti.view(mesh, dom))
        verts = Vector([(point[1], point[2]) for point in element.vals])
        area = order3_triangle_area_from_nodes(verts)
        push!(elements, (Inti.center(element), area))
    end

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)
    
    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    area = 0
    for i in 0:Int32(length(dom_quad)/6)-1
        q_nodes = [dom_quad[6*i+k] for k in 1:6]
        q_poly = [(node.coords[1], node.coords[2]) for node in q_nodes]
        mesh_area = nothing
        for element in elements
            if PolygonAlgorithms.contains(q_poly, (element[1][1], element[1][2]))
                mesh_area = element[2]
                break
            end
        end
        area = 0
        for node in q_nodes
            area = area + node.weight
        end
        @test area ≈ mesh_area atol=1e-15
    end
end

function test_simple_quadtree()
    # mesh will be simple, no holes, but small triangles
    # just testing the internal quadtree part against the triangle representation
    meshsize = 0.025
    tmp_mesh = createMesh(meshsize, meshsize)
    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(tmp_mesh, meshsize, meshsize)
    showMesh(quadtree_mesh)
    # quadtree_views = extractQuadtreeMesh(quadtree_mesh)
    # quadtree_submesh = quadtree_views[1] # because only one quadtree here
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2 && "Quadtree" in Inti.global_get_entity(e).labels, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    quads = []
    for e in Inti.elements(dom_mesh)
        verts = Inti.vertices(e)
        t_verts = [(item[1], item[2]) for item in verts]
        quad = [t_verts[1], t_verts[2], t_verts[3], t_verts[4]]
        push!(quads, quad)
    end

    quad_union = PolygonAlgorithms.union_geometry([quads[i] for i in 1:length(quads)]...)
    println(quad_union)

    triangle_mesh = createMesh(meshsize, meshsize; outside_curve = quad_union[1])
    showMesh(triangle_mesh)
    
    triangle_domain = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, triangle_mesh)
    triangle_boundary = Inti.boundary(triangle_domain)

    triangle_domain_mesh = Inti.view(triangle_mesh, triangle_domain)
    triangle_boundary_mesh = Inti.view(triangle_mesh, triangle_boundary)

    dom_triangle_quad = Inti.Quadrature(triangle_domain_mesh; qorder = 4)
    boundary_triangle_quad = Inti.Quadrature(triangle_boundary_mesh; qorder = 6)
    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)
    
    # now test quadratures

    @testset "Quadtree Mesh" begin
        @test area(dom_quad, boundary_quad) ≈ area(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test zero(dom_quad, boundary_quad) ≈ zero(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test negative_area(dom_quad, boundary_quad) ≈ negative_area(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test linear_x(dom_quad, boundary_quad) ≈ linear_x(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test linear_y(dom_quad, boundary_quad) ≈ linear_y(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test quadratic_x(dom_quad, boundary_quad) ≈ quadratic_x(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test quadratic_y(dom_quad, boundary_quad) ≈ quadratic_y(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test quadratic_xy(dom_quad, boundary_quad) ≈ quadratic_xy(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test high_order_x(dom_quad, boundary_quad) ≈ high_order_x(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test high_order_y(dom_quad, boundary_quad) ≈ high_order_y(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test high_order_xy(dom_quad, boundary_quad) ≈ high_order_xy(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test sin_x(dom_quad, boundary_quad) ≈ sin_x(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test sin_y(dom_quad, boundary_quad) ≈ sin_y(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test sin_xy(dom_quad, boundary_quad) ≈ sin_xy(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test sin_x_cos_y(dom_quad, boundary_quad) ≈ sin_x_cos_y(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test exp_xy(dom_quad, boundary_quad) ≈ exp_xy(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
        @test spike(dom_quad, boundary_quad) ≈ spike(dom_triangle_quad, boundary_triangle_quad) atol=1e-6
    end

    # this test fails on high_order_x, high_order_y, and spike. but why?????
    # not an issue with the quadtree being built for the function or not; it is built for the function here
    # crank up tolerance when returning to this
end

# @testset verbose = true "All Tests" begin
#     test_simple_triangle_mesh()
#     test_simple_quadtree_mesh()
# end;
@testset "hi" begin
    test_simple_triangle_mesh_area()
    test_simple_quadtree_mesh_area()
    test_simple_quadtree_greens_theorem()
    # test_simple_quadtree_mesh_quadratures(0.4)
    # test_simple_quadtree_mesh_quadratures(0.3)
    # test_simple_quadtree_mesh_quadratures(0.2)
    test_simple_quadtree_mesh_quadratures(0.1)
    # test_simple_quadtree_mesh_quadratures(0.075)
    # test_simple_quadtree_mesh_quadratures(0.05)
    # test_simple_quadtree_mesh_quadratures(0.025)
    # test_simple_triangle_mesh_quadratures(0.4)
    # test_simple_triangle_mesh_quadratures(0.3)
    # test_simple_triangle_mesh_quadratures(0.2)
    test_simple_triangle_mesh_quadratures(0.1)
    # test_simple_triangle_mesh_quadratures(0.075)
    # test_simple_triangle_mesh_quadratures(0.05)
    # test_simple_triangle_mesh_quadratures(0.025)
end;

# test_simple_quadtree_greens_theorem()

# meshsize = 0.4
# mesh = createMesh(meshsize, meshsize)
# # showMesh(mesh)
# parametrizations = [(x) -> [cos(2*x*pi), sin(2*x*pi)]]
# boundary_mesh, quadtree_mesh = createQuadtreeMesh(mesh, x -> 1.0, 0.4, 0.4, parametrizations; curve_mesh = true)
# println(Inti.entities(boundary_mesh))

# dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, boundary_mesh)
# quadtree_dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
# dom_quad = Inti.Quadrature(boundary_mesh[dom]; qorder = 10)
# quadtree_dom_quad = Inti.Quadrature(quadtree_mesh[quadtree_dom]; qorder = 10)
# area = Inti.integrate(x -> 1.0, dom_quad) + Inti.integrate(x -> 1.0, quadtree_dom_quad)
# println(area - pi)


# crv_mesh = Inti.curve_mesh(boundary_mesh, parametrizations[1], 6)
# crv_dom_quad = Inti.Quadrature(crv_mesh[dom]; qorder = 10)
# area = Inti.integrate(x -> 1.0, crv_dom_quad) + Inti.integrate(x -> 1.0, quadtree_dom_quad)
# println(area - pi)

# forcing_function = (x) -> 1
# gmsh.initialize()
# quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize; output_split_meshes=true)
# quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize; output_split_meshes=false)

# showMesh(quadtree_mesh)
# meshes = getMeshList()
# for mesh in meshes
    # showMesh(mesh, showBoundary = false)
# end

# meshsize = 0.025
# tmp_mesh = createMesh(meshsize, meshsize)
# gmsh.initialize()
# quadtree_mesh = createQuadtreeMesh(tmp_mesh, meshsize, meshsize)
# # construct parametrization_dict manually

# for e in Inti.entities(quadtree_mesh)
#     if Inti.geometric_dimension(e) == 2
#         labels = Inti.labels(e)
#         for label in labels
#             if label != "Quadtree"
#                 curveMesh(quadtree_mesh, e, parametrization_dict[label])
#             end
#         end
#     end
# end
# curve_mesh = curveMesh(mesh)

# problem:
# must split mesh into separate meshes for each entity
# then curve each mesh (that needs to be curved)
# then combine these meshes into one


# todo for tmr:
# first, write mesh splitting and mesh merging code
# then, write mesh curving code
# next, write a new test function for curved mesh (with higher tols)
# then change test suite to be individual functions with mesh creation within each test (with proper forcing function for quadtrees)
# then, write another test to replicate conditions in other file with inaccuracy
# then dig in further