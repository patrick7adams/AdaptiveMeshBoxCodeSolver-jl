include("../src/main.jl")

using Test

function area_from_boundary_mesh(boundary_mesh)
    area = 0
    for e in Inti.elements(boundary_mesh)
        verts = Vector([(point[1], point[2]) for point in e.vals])
        area += verts[1][1]*verts[2][2] - verts[2][1]*verts[1][2]
    end
    return 0.5 * abs(area)
end

function test_simple_quadtree_vs_triangle_mesh_area()
    meshsize = 0.025
    mesh = createMesh(meshsize, meshsize*2)

    forcing_function = (x) -> 1

    gmsh.initialize()
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize*2, parametrizations)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = get_boundary(quadtree_mesh)

    dom_mesh = Inti.view(quadtree_mesh, dom)
    boundary_mesh = Inti.view(quadtree_mesh, boundary)

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

    # showMesh(mesh)
    # showMesh(quadtree_mesh)

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

    @test result ≈ triangle_result atol=1e-13
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

    @test result ≈ triangle_result atol=1e-13
end