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
    boundary = Inti.boundary(dom)

    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> 1

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ pi atol=1e-6
end

function test_simple_triangle_mesh_zero()
    mesh = createMesh(0.025, 0.025)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary = Inti.boundary(dom)

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
    boundary = Inti.boundary(dom)

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
    boundary = Inti.boundary(dom)

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
    boundary = Inti.boundary(dom)

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
    boundary = Inti.boundary(dom)

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
    boundary = Inti.boundary(dom)

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
    boundary = Inti.boundary(dom)

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
    boundary = Inti.boundary(dom)

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
    boundary = Inti.boundary(dom)

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
    boundary = Inti.boundary(dom)

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
    boundary = Inti.boundary(dom)

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
    boundary = Inti.boundary(dom)

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
    boundary = Inti.boundary(dom)

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
    boundary = Inti.boundary(dom)

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
    boundary = Inti.boundary(dom)

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
    boundary = Inti.boundary(dom)

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
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> 1

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)
    
    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)
    showMeshWithQuadrature(quadtree_mesh, dom_quad)
    
    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)
    
    @test result ≈ pi atol=1e-6
end

function test_simple_quadtree_mesh_zero()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> 0

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

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
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)

    @test result ≈ 0.0314159 atol=1e-6
end

function test_simple_quadtree_greens_theorem()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    P = (x, y) -> x^2*y
    Q = (x, y) -> x^3*y^2
    dPdy = (x, y) -> x^2
    dQdx = (x, y) -> 3*x^2*y^2

    boundary_fn = (q) -> P(q.coords...)*q.normal[2]*-1 + Q(q.coords...)*q.normal[1]
    domain_fn = (q) -> dQdx(q.coords...) - dPdy(q.coords...)

    forcing_function = (x) -> dQdx(x...) - dPdy(x...)

    gmsh.initialize()
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, [(x -> [cos(2*pi*x), sin(2*pi*x)])])

    showMesh(quadtree_mesh)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.external_boundary(dom)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    computed_boundary_int = Inti.integrate(boundary_fn, boundary_quad)
    computed_domain_int = Inti.integrate(domain_fn, dom_quad)

    println(computed_boundary_int - (-pi/8))
    println(computed_domain_int - (-pi/8))
    println("Error in values: ", abs(computed_boundary_int - computed_domain_int))
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
    boundary = Inti.boundary(dom)

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

# test_simple_quadtree_mesh_area()

test_simple_quadtree_greens_theorem()

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