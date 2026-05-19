# a bunch of functions that create a variety of meshes for testing
using Inti, Meshes, CairoMakie
using Gmsh

order = 1

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
            push!(pointTags, gmsh.model.geo.addPoint(points[i]..., 0.0, max_triangle_size)) # TMP
        end
    end
    lineTags = []
    for i in 1:length(pointTags)
        push!(lineTags, gmsh.model.geo.addLine(pointTags[i], pointTags[mod1(i+1, length(pointTags))]))
    end
    curve_loop = gmsh.model.geo.addCurveLoop(lineTags)
    return curve_loop, lineTags
end

function createInternalCurveLoop(points, max_triangle_size)
    pointTags = []
    for i in 1:length(points)
        if distance(points[i], points[mod1(i+1, end)]) > 1e-6
            # first determine the size of the quad the point is associated with
            size = min(distance(points[i], points[mod1(i+1, end)]), distance(points[i], points[mod1(i-1, end)]))
            push!(pointTags, gmsh.model.geo.addPoint(points[i]..., 0.0, size)) # TMP
        end
    end
    lineTags = []
    for i in 1:length(pointTags)
        lineTag = gmsh.model.geo.addLine(pointTags[i], pointTags[mod1(i+1, length(pointTags))])
        gmsh.model.geo.mesh.setTransfiniteCurve(lineTag, 2)
        push!(lineTags, lineTag)
    end
    curve_loop = gmsh.model.geo.addCurveLoop(lineTags)
    return curve_loop, lineTags
end

function createMesh(meshsize, max_triangle_size; outside_curve = [], holes = [])
    # this is good so far. I can now create inscribed mesh polygons from a
    # list of points. Very good, now will work on refinement criteria
    # and return once proper internal mesh is done.
    # Creates and returns a gmsh mesh to use
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 3) # turn to 4/5 for info or debug, 3 is all that is necessary I think though
    loops = []
    if outside_curve == []
        points = []
        for t in range(0, stop = 2pi, length=Int32(ceil(100 / max_triangle_size)))[1:end-1]
            point = (cos(t), sin(t))
            push!(points, gmsh.model.geo.addPoint(point..., 0.0, meshsize))
        end
        push!(points, points[1])
        polyloop = gmsh.model.geo.addPolyline(points)
        curve_loop = gmsh.model.geo.addCurveLoop([polyloop])
        push!(loops, curve_loop)
    else
        curve_loop = createCurveLoop(outside_curve, max_triangle_size)
        push!(loops, curve_loop)
    end

    if holes != []
        for hole in holes
            hole_curve = createCurveLoop(hole, max_triangle_size)
            push!(loops, hole_curve)
        end
    end
    surf = gmsh.model.geo.addPlaneSurface(loops)
    
    gmsh.model.geo.synchronize()
    # gmsh.option.setNumber("Mesh.MeshSizeMax", max_triangle_size)
    # gmsh.option.setNumber("Mesh.HighOrderOptimize", 1)
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