include("../src/main.jl")

using Test

function test_volume_potential_far_x()
    meshsize = 0.1
    mesh = createMesh(meshsize, meshsize)

    u = (x) -> x[1]
    x = (100.0, 0.0)

    gmsh.initialize()
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, u, meshsize, meshsize, parametrizations)

    dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary_dom = Inti.boundary(dom)

    dom_mesh = view(quadtree_mesh, dom)
    boundary_mesh = view(quadtree_mesh, boundary_dom)

    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder = 6)

    volume_potential = calculateVolumePotential(dom_mesh, u)

    println(volume_potential(x))
    
    estimated_result = log(x[1]^2+x[2]^2)/(2pi) * Inti.integrate((q) -> u(q.coords), dom_quad)
    println(estimated_result)
end

test_volume_potential_far_x()