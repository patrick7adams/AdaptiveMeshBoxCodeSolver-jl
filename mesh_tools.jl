# a bunch of functions that create a variety of meshes for testing
using Inti, Meshes, CairoMakie
using Gmsh

order = 3

function distance(a, b)
    return sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)
end

function showMesh(msh; showBoundary = true)
    println("Showing mesh!!!")
    Ω = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, msh)
    Ω_msh = @views msh[Ω]
    fig = viz(
        Ω_msh;
        segmentsize = 1,
        showsegments = true,
        axis = (aspect = DataAspect(),),
        figure = (; size = (500, 400)),
    )
    if showBoundary
        Γ = Inti.external_boundary(Ω)
        Γ_msh = @views msh[Γ]
        viz!(Γ_msh; color = :red, segmentsize = 1)
    end
    display(fig)
end

function showMeshWithQuadrature(msh, quad)
    println("Showing mesh!!!")
    Ω = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, msh)
    Γ = Inti.boundary(Ω)
    Ω_msh = @views msh[Ω]
    Γ_msh = @views msh[Γ]
    domain_quad_coords = Meshes.PointSet(map((q) -> Meshes.Point(q.coords...), quad))
    fig = viz(
        Ω_msh;
        segmentsize = 1,
        showsegments = true,
        axis = (aspect = DataAspect(),),
        figure = (; size = (500, 400)),
    )
    viz!(Γ_msh; color = :red, segmentsize = 1)
    viz!(domain_quad_coords; color = :blue, pointsize = 2)
    display(fig)
end

function createCurveLoop(points, max_triangle_size)
    pointTags = []
    for i in 1:length(points)
        if distance(points[i], points[mod1(i+1, end)]) > 1e-6
            push!(pointTags, gmsh.model.occ.addPoint(points[i]..., 0.0, max_triangle_size)) # TMP
        end
    end
    lineTags = []
    for i in 1:length(pointTags)-1
        push!(lineTags, gmsh.model.occ.addLine(pointTags[i], pointTags[i+1]))
    end
    push!(lineTags, gmsh.model.occ.addLine(pointTags[end], pointTags[1]))
    return lineTags, gmsh.model.occ.addCurveLoop(lineTags)
end

function createMesh(meshsize, max_triangle_size; outside_curve = [], holes = [])
    # this is good so far. I can now create inscribed mesh polygons from a
    # list of points. Very good, now will work on refinement criteria
    # and return once proper internal mesh is done.
    # Creates and returns a gmsh mesh to use
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 3) # turn to 4/5 for info or debug, 3 is all that is necessary I think though

    if outside_curve == []
        circle = Inti.gmsh_curve(0, 2pi; npts=200, meshsize=max_triangle_size) do s # this could be the problem? maybe not approximating this curve fully
            return Inti.Point2D(cos(s), sin(s))
        end
        # circle = gmsh.model.occ.addCircle(0.0, 0.0, 0.0, 1.0)
        curve_loop = gmsh.model.occ.addCurveLoop([circle])
    else
        curve_loop = createCurveLoop(outside_curve, max_triangle_size)
    end

    loops = [curve_loop]
    if holes != []
        for hole in holes
            hole_curve = createCurveLoop(hole, max_triangle_size)
            push!(loops, hole_curve)
        end
    end
    surf = gmsh.model.occ.addPlaneSurface(loops)
    
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", max_triangle_size)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 1)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(order)
    mesh = Inti.import_mesh(; dim = 2)
    gmsh.finalize()
    return mesh
end

function getLabels()
    labels = ["Quadtree"]
    for i in 1:num_boundaries
        push!(labels, string("Boundary ", i))
    end
    return labels
end

function getMeshList()
    gmsh.initialize()
    meshes = []
    for name in getLabels()
        gmsh.open(string(name, ".msh"))
        msh = Inti.import_mesh(; dim = 2)
        push!(meshes, msh)
    end
    gmsh.finalize()
    return meshes
end

function get_boundary(mesh)
    # gets the boundary of the mesh according to my label system
    labels = getLabels()[2:end]
    return Inti.Domain((e) -> Inti.geometric_dimension(e) == 1 && length(Inti.labels(e)) > 0 && Inti.labels(e)[1] in labels, mesh)
end

using StaticArrays

function cubic_lagrange_basis(t)
    ts = (0.0, 1/3, 2/3, 1.0)

    L = zeros(4)
    dL = zeros(4)

    for i in 1:4
        # L_i(t)
        num = 1.0
        den = 1.0
        for j in 1:4
            if j != i
                num *= t - ts[j]
                den *= ts[i] - ts[j]
            end
        end
        L[i] = num / den

        # derivative L_i'(t)
        s = 0.0
        for k in 1:4
            if k != i
                prod = 1.0
                for j in 1:4
                    if j != i && j != k
                        prod *= t - ts[j]
                    end
                end
                s += prod
            end
        end
        dL[i] = s / den
    end

    return L, dL
end

function cubic_edge_area_contribution(edge)
    # Exact enough: 5-point Gauss-Legendre integrates degree <= 9 exactly.
    # x(t)y'(t)-y(t)x'(t) is degree <= 4 for cubic edges.
    qx = [
        0.5 - 0.5*sqrt(5 + 2sqrt(10/7))/3,
        0.5 - 0.5*sqrt(5 - 2sqrt(10/7))/3,
        0.5,
        0.5 + 0.5*sqrt(5 - 2sqrt(10/7))/3,
        0.5 + 0.5*sqrt(5 + 2sqrt(10/7))/3,
    ]

    qw = [
        (322 - 13sqrt(70)) / 1800,
        (322 + 13sqrt(70)) / 1800,
        64 / 225,
        (322 + 13sqrt(70)) / 1800,
        (322 - 13sqrt(70)) / 1800,
    ]

    val = 0.0

    for (t, w) in zip(qx, qw)
        L, dL = cubic_lagrange_basis(t)

        x = sum(L[i]  * edge[i][1] for i in 1:4)
        y = sum(L[i]  * edge[i][2] for i in 1:4)
        dx = sum(dL[i] * edge[i][1] for i in 1:4)
        dy = sum(dL[i] * edge[i][2] for i in 1:4)

        val += w * (x * dy - y * dx)
    end

    return val
end

function order3_triangle_area_from_nodes(pts)
    # calculates the exact order 3 triangle area with green's thm
    edge1 = [pts[1],  pts[2], pts[3], pts[4]]
    edge2 = [pts[4],  pts[7], pts[9], pts[10]]
    edge3 = [pts[10], pts[8], pts[5], pts[1]]

    A = 0.5 * (
        cubic_edge_area_contribution(edge1) +
        cubic_edge_area_contribution(edge2) +
        cubic_edge_area_contribution(edge3)
    )

    return abs(A)
end