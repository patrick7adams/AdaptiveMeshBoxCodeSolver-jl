order = 1
greens_fn = (r, x) -> 1/(2pi) * log(AdaptiveMeshSolver.distance(x, r))

function distance(a, b)
    return sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)
end

function showMesh(msh, points=[])
    println("Showing mesh!!!")

    Ω = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, msh)
    Ω_msh = @views msh[Ω]
    # fig = Inti.plot(Ω_msh; linewidth=3)
    fig, ax, plt = plot(Ω_msh; strokewidth = 1, axis = (aspect = DataAspect(),), figure=(; size = (1000, 1000)))
    Γ = Inti.boundary(Ω)
    # Γ = Inti.boundary(Ω)
    if length(keys(Γ)) > 0 # boundary exists
        Γ_msh = @views msh[Γ]
        plot!(Γ_msh; color = :red, strokewidth = 1)
    end
    if length(points) > 0
        # pointset = Meshes.PointSet([Meshes.Point((point[1], point[2])) for point in points])
        pointset = Point2f.(points)
        scatter!(pointset; color = :red)
    end
    display(fig)
end

function showErrorMesh(meshes, quadrature_points, errors)
    println("Showing error mesh!!!")
    quadtree_dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, meshes[1])
    quadtree_mesh = view(meshes[1], quadtree_dom)

    quadrature_pointset = Point2f.(quadrature_points)
    xs, ys = ([q[i] for q in quadrature_points] for i in 1:2)
    colors = log10.(errors)

    fig, ax, plt = tricontourf(
        xs, ys, colors; 
        colormap = :inferno, 
        levels = 20, 
        axis = (aspect = DataAspect(),), 
        figure = (; size = (1000, 1000))
    )
    Colorbar(fig[1, 2]; colormap = :inferno, colorrange = (minimum(colors), maximum(colors)), label = "log_{10} of error")
    
    display(fig)

    plot!(quadtree_mesh; strokewidth = 1, color = (:gray, 0.0), shading = NoShading)

    for mesh in meshes[2:end]
        strip_dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, mesh)
        strip_boundary = get_boundary(mesh)

        strip_dom_mesh = view(mesh, strip_dom)
        strip_boundary_mesh = view(mesh, strip_boundary)

        plot!(strip_dom_mesh; strokewidth = 1, color = (:gray, 0.0), shading = NoShading)
        plot!(strip_boundary_mesh; color = :red, strokewidth = 1)
    end

    display(fig)
end

function showSeparatedMesh(meshes, singular_elems, quad, point)
    println("Showing Separated mesh!!!")

    Ω = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, meshes[1])
    Ω_msh = @views meshes[1][Ω]
    fig = viz(
        Ω_msh;
        segmentsize = 1,
        showsegments = true,
        axis = (aspect = DataAspect(),),
        figure = (; size = (400, 400)),
    )
    for mesh in meshes[2:end]
        strip_dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, mesh)
        strip_boundary = get_boundary(mesh)

        strip_dom_mesh = view(mesh, strip_dom)
        strip_boundary_mesh = view(mesh, strip_boundary)

        viz!(strip_dom_mesh; showsegments = true)
        viz!(strip_boundary_mesh; color = :red, segmentsize = 1)
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
    viz!(Meshes.Point(point); color = :green, pointsize = 10)

    point_list = [(0.2908767452425044, 0.6752146835957239),(0.2908767452425044, 0.5754898716685045), (0.4166184434922614, 0.6776414331925313), (0.35231495232434595, 0.627674221001969), (0.35231495232434595, 0.6741279735048046), (0.2937420293285783, 0.626543793949986)]
    for pnt in point_list
        viz!(Meshes.Point(pnt); color = :red, pointsize = 4)
    end
    display(fig)
end

function showMeshes(meshes)
    println("Showing combined mesh!!!")
    quadtree_dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, meshes[1])
    quadtree_mesh = view(meshes[1], quadtree_dom)
    fig, ax, plt = plot(quadtree_mesh; strokewidth = 1, axis = (aspect = DataAspect(),), figure=(; size = (1000, 1000)))
    for mesh in meshes[2:end]
        strip_dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, mesh)
        # strip_boundary = get_boundary(mesh)
        strip_boundary = Inti.boundary(strip_dom)

        strip_dom_mesh = view(mesh, strip_dom)
        strip_boundary_mesh = view(mesh, strip_boundary)

        plot!(strip_dom_mesh; strokewidth = 1)
        plot!(strip_boundary_mesh; color = :red, strokewidth = 1)
    end
    display(fig)
end

# function showMeshesHeatmap(meshes, )

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
        viz!(line; color = :red, segmentsize = 1, showpoints = true, pointsize = 3)
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
        figure = (; size = (1000, 1000)),
    )
    for edge in edges
        line = Meshes.Segment(Meshes.Point(edge[1][1], edge[1][2]), Meshes.Point(edge[2][1], edge[2][2]))
        # println(line)
        viz!(line; color = :red, segmentsize = 1, showpoints = true, pointsize = 3)
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

function getMultiplicativeTerm(target_points, meshes)
    terms = []
    
    # check if point in quad mesh
    quad_mesh = meshes[1]
    # quad_terms = [:outside for x in target_points]
    quad_terms = [:outside for x in target_points]
    for (i, point) in enumerate(target_points)
        for quad in Inti.elements(quad_mesh)
            coords = [(x[1], x[2]) for x in Inti.vertices(quad)]
            if (coords[1][1] <= point[1] && point[1] <= coords[3][1]) && (coords[1][2] <= point[2] && point[2] <= coords[3][2])
                quad_terms[i] = :inside
                break
            end
        end
    end
    push!(terms, quad_terms)

    # first check if on boundary
    for mesh in meshes[2:end]
        mesh_terms = [:outside for x in target_points]
        # quadrature = getBoundaryQuadrature(mesh, 12)
        domain = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, mesh)
        boundary_dom = Inti.boundary(domain)
        boundary_mesh = Inti.view(mesh, boundary_dom)
        quadrature = Inti.Quadrature(boundary_mesh; qorder = order)
        for (i, point) in enumerate(target_points)
            mult = Inti._green_multiplier(point, quadrature)
            if mult < -0.5
                mesh_terms[i] = :inside
            end
        end
        push!(terms, mesh_terms)
    end
    # now 1 if inside
    return terms
end

function createMesh()
    gmsh.initialize()
    gmsh.clear()
    parametrizations::Vector{Function} = [(x) -> (cos(2*x*pi), sin(2*x*pi))]
    splines = createSplines(parametrizations)
    curve = gmsh.model.occ.addCurveLoop(splines)
    surface = gmsh.model.occ.addPlaneSurface([curve])
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(1, splines, -1, "Boundary")
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 30)
    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.mesh.generate(2)
    mesh = Inti.import_mesh(; dim=2)
    gmsh.finalize()
    return mesh
end

function get_boundary(mesh)
    return Inti.Domain((e) -> Inti.geometric_dimension(e) == 1 && "Boundary" in Inti.labels(e), mesh) 
end