include("../src/main.jl")

using Test

function test_simple_quadtree_mesh_area()
    meshsize = 0.025

    forcing_function = (x) -> 1

    gmsh.initialize()
    parametrizations::Vector{Function} = [(x) -> (cos(x*2*pi), sin(x*2*pi))]
    quadtree_mesh = createQuadtreeMesh(parametrizations, forcing_function, meshsize, meshsize)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    dom_mesh = Inti.view(quadtree_mesh, dom)
    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)

    showMesh(quadtree_mesh)
    
    test_function = (q) -> forcing_function((q.coords[1], q.coords[2]))

    result = Inti.integrate(test_function, dom_quad)
    
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