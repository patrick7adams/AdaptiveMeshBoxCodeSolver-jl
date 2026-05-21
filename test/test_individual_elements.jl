include("../src/main.jl")

using Test

function test_simple_quadtree_mesh_quadratures(meshsize)
    mesh = createMesh(meshsize, meshsize)

    forcing_function = (x) -> 1

    gmsh.initialize()
    parametrizations::Vector{Function} = [(x) -> [cos(x*2*pi), sin(x*2*pi)]]
    quadtree_mesh = createQuadtreeMesh(mesh, forcing_function, meshsize, meshsize, parametrizations)

    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, quadtree_mesh)
    boundary = Inti.boundary(dom)
        
    center_to_elements = Dict()
    for element in Inti.elements(Inti.view(quadtree_mesh, dom))
        if typeof(element) == Inti.LagrangeElement{Inti.ReferenceHyperCube{2}, 4, SVector{2, Float64}}
            verts = Inti.vertices(element)
            area = (verts[3][1] - verts[1][1])*(verts[3][2] - verts[1][2])
            center_to_elements[Inti.center(element)] = area
        end
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