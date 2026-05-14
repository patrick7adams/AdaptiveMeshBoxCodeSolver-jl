# a bunch of functions that create a variety of meshes for testing
using Inti, Meshes, CairoMakie
using Gmsh

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
        Γ = Inti.boundary(Ω)
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
        circle = Inti.gmsh_curve(0, 2pi;meshsize) do s
            return Inti.Point2D(cos(s), sin(s))
        end
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
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(1)
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