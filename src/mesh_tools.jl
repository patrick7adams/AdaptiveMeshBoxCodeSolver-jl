order = 1
greens_fn = (r, x) -> 1/(2pi) * log(AdaptiveMeshSolver.distance(x, r))

function distance(a, b)
    return sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)
end

function showMesh(msh)
    println("Showing mesh!!!")

    Ω = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, msh)
    Ω_msh = @views msh[Ω]
    fig = viz(
        Ω_msh;
        segmentsize = 1,
        showsegments = true,
        axis = (aspect = DataAspect(),),
        figure = (; size = (400, 400)),
    )
    Γ = get_boundary(msh)
    if length(keys(Γ)) > 0 # boundary exists
        Γ_msh = @views msh[Γ]
        viz!(Γ_msh; color = :red, segmentsize = 1)
    end
    display(fig)
end

function showSeparatedMesh(msh, singular_elems, near_singular_elems, quad)
    println("Showing Separated mesh!!!")

    Ω = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, msh)
    Ω_msh = @views msh[Ω]
    fig = viz(
        Ω_msh;
        segmentsize = 1,
        showsegments = true,
        axis = (aspect = DataAspect(),),
        figure = (; size = (400, 400)),
    )
    if length(near_singular_elems) > 0
        near_singular_coords = []
        for elem in near_singular_elems
            append!(near_singular_coords, [(x[1], x[2]) for x in Inti.vertices(elem)])
        end
        near_singular_pointset = Meshes.PointSet(near_singular_coords)
        viz!(near_singular_pointset; color = :green, pointsize = 6, showpoints = true)
    end
    if length(singular_elems) > 0
        singular_coords = []
        for elem in singular_elems
            append!(singular_coords, [(x[1], x[2]) for x in Inti.vertices(elem)])
        end
        singular_pointset = Meshes.PointSet(singular_coords)
        viz!(singular_pointset; color = :red, pointsize = 6, showpoints = true)
    end
    domain_quad_coords = Meshes.PointSet(map((q) -> Meshes.Point(q.coords...), quad))
    viz!(domain_quad_coords; color = :blue, pointsize = 2)
    display(fig)
end

function showMeshes(meshes)
    println("Showing combined mesh!!!")
    quadtree_dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, meshes[1])
    quadtree_mesh = view(meshes[1], quadtree_dom)
    fig = viz(
        quadtree_mesh;
        segmentsize = 1,
        showsegments = true,
        axis = (aspect = DataAspect(),),
        figure = (; size = (1000, 1000)),
    )
    for mesh in meshes[2:end]
        strip_dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, mesh)
        strip_boundary = get_boundary(mesh)

        strip_dom_mesh = view(mesh, strip_dom)
        strip_boundary_mesh = view(mesh, strip_boundary)

        viz!(strip_dom_mesh; showsegments = true)
        viz!(strip_boundary_mesh; color = :red, segmentsize = 1)
    end
    display(fig)
end


function showMeshWithQuadrature(msh, quad)
    println("Showing mesh!!!")
    Ω = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, msh)
    Ω_msh = @views msh[Ω]
    domain_quad_coords = Meshes.PointSet(map((q) -> Meshes.Point(q.coords...), quad))
    fig = viz(
        Ω_msh;
        segmentsize = 1,
        showsegments = true,
        axis = (aspect = DataAspect(),),
        figure = (; size = (500, 400)),
    )
    viz!(domain_quad_coords; color = :blue, pointsize = 2)
    display(fig)
end

function showEdgeList(points::Vector{Tuple{Float64, Float64}})
    println("Showing edge vector")
    line1 = Meshes.Segment(Meshes.Point(points[1][1], points[1][2]), Meshes.Point(points[2][1], points[2][2]))
    fig = viz(
        line1;
        color = :red,
        segmentsize = 1,
        showsegments = true,
        axis = (aspect = DataAspect(),),
        figure = (; size = (500, 400)),
    )
    for i in 2:length(points)-1
        line = Meshes.Segment(Meshes.Point(points[i][1], points[i][2]), Meshes.Point(points[i+1][1], points[i+1][2]))
        # println(line)
        viz!(line; color = :red, segmentsize = 1, showpoints = true, pointsize = 6)
    end
    display(fig)
end

function showEdgeList(edges::Set{Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}})
    println("Showing edge set")
    first_edge = first(edges)
    line1 = Meshes.Segment(Meshes.Point(first_edge[1][1], first_edge[1][2]), Meshes.Point(first_edge[2][1], first_edge[2][2]))
    fig = viz(
        line1;
        color = :red,
        segmentsize = 1,
        showsegments = true,
        axis = (aspect = DataAspect(),),
        figure = (; size = (500, 400)),
    )
    for edge in edges
        line = Meshes.Segment(Meshes.Point(edge[1][1], edge[1][2]), Meshes.Point(edge[2][1], edge[2][2]))
        # println(line)
        viz!(line; color = :red, segmentsize = 1, showpoints = true, pointsize = 6)
    end
    display(fig)
end

function createCurveLoop(points)
    pointTags = []
    for i in 1:length(points)
        if distance(points[i], points[mod1(i+1, end)]) > 1e-6
            # first determine the size of the quad the point is associated with
            size = min(distance(points[i], points[mod1(i+1, end)]), distance(points[i], points[mod1(i-1, end)]))
            push!(pointTags, gmsh.model.occ.addPoint(points[i]..., 0.0, size)) # TMP
        end
    end
    lineTags = []
    for i in 1:length(pointTags)
        lineTag = gmsh.model.occ.addLine(pointTags[i], pointTags[mod1(i+1, length(pointTags))])
        push!(lineTags, lineTag)
    end
    gmsh.model.occ.synchronize()
    curve_loop = gmsh.model.occ.addCurveLoop(lineTags)
    return curve_loop, lineTags
end

function createMesh(meshsize; outside_curve = [], holes = [])
    # this is good so far. I can now create inscribed mesh polygons from a
    # list of points. Very good, now will work on refinement criteria
    # and return once proper internal mesh is done.
    # Creates and returns a gmsh mesh to use
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 3) # turn to 4/5 for info or debug, 3 is all that is necessary I think though
    loops = []
    if outside_curve == []
        # circle = Inti.gmsh_curve((x) -> Inti.Point2D(cos(2*pi*x), sin(2*pi*x)), 0.0, 1.0; meshsize = meshsize)
        points = []
        for t in range(0, stop = 1, length=Int32(ceil(100 / meshsize)))[1:end-1]
            point = (cos(t*2*pi), sin(t*2*pi))
            push!(points, gmsh.model.geo.addPoint(point..., 0.0, meshsize))
        end
        push!(points, points[1])
        polyloop = gmsh.model.geo.addPolyline(points)
        curve_loop = gmsh.model.geo.addCurveLoop([polyloop])
        push!(loops, curve_loop)
    else
        curve_loop = createCurveLoop(outside_curve)
        push!(loops, curve_loop)
    end

    if holes != []
        for hole in holes
            hole_curve = createCurveLoop(hole)
            push!(loops, hole_curve)
        end
    end
    surf = gmsh.model.geo.addPlaneSurface(loops)
    
    gmsh.model.geo.synchronize()
    # gmsh.option.setNumber("Mesh.HighOrderOptimize", 1)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(order)
    mesh = Inti.import_mesh(; dim = 2)
    gmsh.finalize()
    return mesh
end

function get_boundary(mesh)
    return Inti.Domain((e) -> Inti.geometric_dimension(e) == 1 && "Boundary" in Inti.labels(e), mesh) 
end