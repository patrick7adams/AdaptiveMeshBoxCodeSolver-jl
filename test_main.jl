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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

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

function area_from_boundary_mesh(boundary_mesh)
    area = 0
    for e in Inti.elements(boundary_mesh)
        verts = Vector([(point[1], point[2]) for point in e.vals])
        area += verts[1][1]*verts[2][2] - verts[2][1]*verts[1][2]
    end
    return 0.5 * abs(area)
end

function test_simple_quadtree_vs_triangle_mesh_area()
    meshsize = 0.2
    mesh = createMesh(meshsize, meshsize*2)

    forcing_function = (x) -> 1

    gmsh.initialize()
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize*2, parametrizations)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    # for element in Inti.elements(dom_mesh)
    #     if typeof(element) == Inti.LagrangeElement{Inti.ReferenceSimplex{2}, 3, SVector{2, Float64}}
    #         println("PP")
    #     elseif typeof(element) == Inti.LagrangeElement{Inti.ReferenceHyperCube{2}, 4, SVector{2, Float64}}
    #         println("QQ")
    #     end
    # end
    E = Inti.LagrangeElement{Inti.ReferenceHyperCube{2}, 4, SVector{2, Float64}}
    polygon_points = map((x) -> Inti.vertices(x), Inti.elements(dom_mesh, E))
    typed_polygon_points = Vector{Vector{Tuple{Float64, Float64}}}()
    x_points = []
    for polygon in polygon_points
        points = Vector{Tuple{Float64, Float64}}()
        for point in polygon
            push!(points, (point[1], point[2]))
            push!(x_points, point[2])
        end
        push!(typed_polygon_points, points)
    end
    boundary = PolygonAlgorithms.union_geometry(typed_polygon_points...)[1]
    # println(boundary)
    # println(PolygonAlgorithms.y_coords(boundary))
    # println(maximum(x_points))
    # println(maximum(PolygonAlgorithms.y_coords(boundary)) - maximum(x_points))
    # println(PolygonAlgorithms.y_coords(boundary))
    # println(PolygonAlgorithms.y_coords(boundary)[1][1] + PolygonAlgorithms.y_coords(boundary)[1][2])
    E = Inti.LagrangeElement{Inti.ReferenceSimplex{2}, 3, SVector{2, Float64}}
    polygon_points = map((x) -> Inti.vertices(x), Inti.elements(dom_mesh, E))
    for polygon in polygon_points
        for point in polygon
            # if PolygonAlgorithms.contains(boundary, (point[1], point[2]), atol=1e-15)
            #     println(point)
            # end
        end
    end

    # typed_polygon_points = [PolygonAlgorithms.Polygon([Tuple(point[1], point[2]) for point in polygon]) for polygon in polygon_points] 
    # println(typed_polygon_points)
    # println(PolygonAlgorithms.union_geometry(typed_polygon_points))
    # boundary_points = [PolygonAlgorithms.union_geometry(typed_polygon_points[i]...) for polygon in polygon_points if length(remove_points[i]) > 0]

    boundary_mesh_area = area_from_boundary_mesh(boundary_mesh)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    triangle_dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    triangle_boundary = Inti.boundary(triangle_dom)

    triangle_dom_mesh = Inti.view(mesh, triangle_dom)
    triangle_boundary_mesh = Inti.view(mesh, triangle_boundary)

    triangle_boundary_mesh_area = area_from_boundary_mesh(triangle_boundary_mesh)

    triangle_dom_quad = Inti.Quadrature(triangle_dom_mesh; qorder = 4)
    triangle_boundary_quad = Inti.Quadrature(triangle_boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    domain_mesh_area = Inti.integrate(test_function, dom_quad)
    triangle_domain_mesh_area = Inti.integrate(test_function, triangle_dom_quad)

    showMesh(mesh)
    showMesh(quadtree_mesh)

    println("Triangle mesh domain area = ", triangle_domain_mesh_area)
    println("Quadtree mesh domain area = ", domain_mesh_area)
    println("Difference in mesh domain areas = ", abs(triangle_domain_mesh_area - domain_mesh_area))
    println("Triangle mesh boundary area = ", triangle_boundary_mesh_area)
    println("Quadtree mesh boundary area = ", boundary_mesh_area)
    println("Difference in mesh boundary areas = ", abs(triangle_boundary_mesh_area - boundary_mesh_area))
    println("Triangle mesh domain and boundary area difference = ", abs(triangle_domain_mesh_area - triangle_boundary_mesh_area))
    println("Quadtree mesh domain and boundary area difference = ", abs(domain_mesh_area - boundary_mesh_area))
end

function test_simple_quadtree_vs_triangle_mesh_linear_x()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> x[1]

    gmsh.initialize()
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    triangle_dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    triangle_boundary = Inti.boundary(triangle_dom)

    triangle_dom_mesh = Inti.view(mesh, triangle_dom)
    triangle_boundary_mesh = Inti.view(mesh, triangle_boundary)

    triangle_dom_quad = Inti.Quadrature(triangle_dom_mesh; qorder = 4)
    triangle_boundary_quad = Inti.Quadrature(triangle_boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)
    triangle_result = Inti.integrate(test_function, triangle_dom_quad)

    @test result ≈ triangle_result atol=1e-15
end

function test_simple_quadtree_vs_triangle_mesh_linear_y()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> x[2]

    gmsh.initialize()
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    triangle_dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    triangle_boundary = Inti.boundary(triangle_dom)

    triangle_dom_mesh = Inti.view(mesh, triangle_dom)
    triangle_boundary_mesh = Inti.view(mesh, triangle_boundary)

    triangle_dom_quad = Inti.Quadrature(triangle_dom_mesh; qorder = 4)
    triangle_boundary_quad = Inti.Quadrature(triangle_boundary_mesh; qorder = 6)

    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)
    triangle_result = Inti.integrate(test_function, triangle_dom_quad)

    @test result ≈ triangle_result atol=1e-15
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

function test_simple_quadtree_greens_third_identity_sin_term()
    meshsize = 0.05
    mesh = createMesh(meshsize, meshsize)
    # showMesh(mesh)

    # for k in [0.25, 0.5, 1.0, 2.0, 4.0]
    for k in [4.0]
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

        showMesh(quadtree_mesh)

        dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
        boundary = get_boundary(quadtree_mesh)

        dom_mesh = Inti.view(quadtree_mesh, dom)
        println(dom_mesh)
        boundary_mesh = Inti.view(quadtree_mesh, boundary)

        dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
        boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

        triangle_dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
        triangle_dom_mesh = Inti.view(mesh, triangle_dom)
        triangle_quad = Inti.Quadrature(triangle_dom_mesh; qorder = 17) # 17 is highest order

        # boundary_integral = Inti.integrate(boundary_function, boundary_quad)
        boundary_int1 = Inti.integrate(boundary_fn1, boundary_quad)
        boundary_int2 = Inti.integrate(boundary_fn2, boundary_quad)
        boundary_integral = Inti.integrate(boundary_function, boundary_quad)
        domain_integral = Inti.integrate(domain_function, dom_quad)
        calculated_u_val = domain_integral - boundary_integral

        # exact_domain_integral = Inti.integrate(domain_function, triangle_quad)
        dom_quad_highorder = Inti.Quadrature(dom_mesh; qorder = 17)
        exact_domain_integral = Inti.integrate(domain_function, dom_quad_highorder)
        exact_triangle_domain_integral = Inti.integrate(domain_function, triangle_quad)

        # vals_lowres = Dict{DataType, Float64}()
        # vals_highres = Dict{DataType, Float64}()
        # for E in Inti.element_types(dom_mesh)
        #     vals_lowres[E] = 0.0
        #     vals_highres[E] = 0.0
        #     qtags_lowres = dom_quad.etype2qtags[E]
        #     qtags_highres = dom_quad_highorder.etype2qtags[E]
        #     for row in qtags_lowres
        #         for val in row
        #             qtag = val
        #             q = dom_quad[qtag]
        #             vals_lowres[E] += domain_function(q) * q.weight
        #         end
        #     end
        #     for row in qtags_highres
        #         for val in row
        #             qtag = val
        #             q = dom_quad_highorder[qtag]
        #             vals_highres[E] += domain_function(q) * q.weight
        #         end
        #     end
        # end
        # for E in Inti.element_types(dom_mesh)
        #     println("Element type = ", E)
        #     println("Lowres value for type = ", vals_lowres[E])
        #     println("highres value for type = ", vals_highres[E])
        #     println("Error in values = ", abs(vals_highres[E] - vals_lowres[E]))
        # end

        total_error = abs(calculated_u_val)
        quadrature_error = abs(domain_integral - exact_domain_integral)
        # this quadrature error can get down very low when degree of interpolating polynomial is low
        # thus, when the quadtree is highly refined, increasing quadrature points significantly does not decrease error
        geometry_error = abs(exact_triangle_domain_integral - exact_domain_integral)
        # this is the current dominant error, the error between the integral taken over the triangle mesh and
        # the error on the proper quadtree domain. 

        println("---------")
        println("k = ", k)
        println("boundary_int1 = ", boundary_int1)
        println("boundary_int2 = ", boundary_int2)
        println("Whole boundary integral = ", boundary_integral)
        println("Domain integral = ", domain_integral)
        println("Quadrature error = ", quadrature_error)
        println("Geometry error = ", geometry_error)
        println("Total error = ", total_error)
        println("Quadrature error + geometry error = ", quadrature_error + geometry_error)
        println("Boundary integral = ", Inti.integrate(boundary_function, boundary_quad))
        println("Domain Integral = ", Inti.integrate(domain_function, dom_quad))
        println("Difference = ", calculated_u_val)
        @test 0.0 ≈ calculated_u_val atol=1e-15
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

function test_simple_quadtree_mesh_quadratures(meshsize)
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> 1

    gmsh.initialize()
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations; return_quadtree_mesh=true)

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
        area = 0.5 * (verts[1][1] * (verts[2][2] - verts[3][2]) + verts[2][1] * (verts[3][2] - verts[1][2]) + verts[3][1] * (verts[1][2] - verts[2][2]))
        # if [Inti.center(element)[1], Inti.center(element)[2]] == [-0.35404532580101067, 0.7382086437260025]
        #     println(verts)
        # end
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
            if PolygonAlgorithms.contains(q_poly, (element[1][1], element[1][2]), atol=0.0)
                mesh_area = element[2]
                # if [element[1][1], element[1][2]] == [-0.35404532580101067, 0.7382086437260025]
                #     println(element[1])
                #     println(q_poly)
                #     # println(q_nodes)
                # end
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

# tests take hella long but they all pass!

# @testset "hi" begin 
#     # test_simple_triangle_mesh_area()
#     # test_simple_quadtree_mesh_area()
    # test_simple_quadtree_greens_theorem()
    # test_simple_quadtree_greens_third_identity_simple_forcing()
    # test_simple_quadtree_greens_third_identity_linear_x()
    # test_simple_quadtree_greens_third_identity_linear_y()
    # test_simple_quadtree_greens_third_identity_quadratic()
    # test_simple_quadtree_greens_third_identity_sin_term()
    # test_simple_quadtree_vs_triangle_mesh_area()
    # test_simple_quadtree_vs_triangle_mesh_linear_x()
    # test_simple_quadtree_vs_triangle_mesh_linear_y()
    # test_simple_quadtree_greens_third_identity_complex_forcing()
    # test_simple_quadtree_mesh_quadratures(0.4)
    # test_simple_quadtree_mesh_quadratures(0.3)
    # test_simple_quadtree_mesh_quadratures(0.2)
    # test_simple_quadtree_mesh_quadratures(0.1)
    # test_simple_quadtree_mesh_quadratures(0.075)
    # test_simple_quadtree_mesh_quadratures(0.05)
    # test_simple_quadtree_mesh_quadratures(0.025)
    # test_simple_triangle_mesh_quadratures(0.4)
    # test_simple_triangle_mesh_quadratures(0.3)
    # test_simple_triangle_mesh_quadratures(0.2)
    # test_simple_triangle_mesh_quadratures(0.1)
    # test_simple_triangle_mesh_quadratures(0.075)
    # test_simple_triangle_mesh_quadratures(0.05)
    # test_simple_triangle_mesh_quadratures(0.025)
# end;

@testset "test_greens_identity" begin
    # test_simple_quadtree_greens_third_identity_simple_forcing()
    # test_simple_quadtree_greens_third_identity_linear_x()
    # test_simple_quadtree_greens_third_identity_linear_y()
    # test_simple_quadtree_greens_third_identity_quadratic()
    # test_simple_quadtree_greens_third_identity_sin_term()
    # test_simple_quadtree_greens_third_identity_complex_forcing()
    test_simple_quadtree_vs_triangle_mesh_area()
end;

# test_simple_quadtree_greens_theorem()

# meshsize = 0.4
# mesh = createMesh(meshsize, meshsize)
# # showMesh(mesh)
# gmsh.initialize()
# parametrizations::Vector{Function} = [(x) -> [cos(2*x*pi), sin(2*x*pi)]]
# quadtree_mesh = createQuadtreeMesh(mesh, x -> 1.0, 0.4, 0.4, parametrizations)
# println(Inti.entities(quadtree_mesh))

# dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
# # quadtree_dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
# # dom_quad = Inti.Quadrature(boundary_mesh[dom]; qorder = 10)
# # quadtree_dom_quad = Inti.Quadrature(quadtree_mesh[quadtree_dom]; qorder = 10)
# # area = Inti.integrate(x -> 1.0, dom_quad) + Inti.integrate(x -> 1.0, quadtree_dom_quad)
# # println(area - pi)
# dict = Dict{Inti.EntityKey, Function}()
# for entity in Inti.entities(quadtree_mesh)
#     if Inti.geometric_dimension(entity) == 2
#         l = Inti.labels(entity)
#         if "Quadtree" in l
#             dict[entity] = (x) -> [x, x]
#         elseif "Boundary 1" in l
#             dict[entity] = parametrizations[1]
#         end
#     end
# end

# crv_mesh = Inti.curve_mesh(quadtree_mesh, dict, 6)
# crv_dom_quad = Inti.Quadrature(crv_mesh[dom]; qorder = 10)
# area = Inti.integrate(x -> 1.0, crv_dom_quad)
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