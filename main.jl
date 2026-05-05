# Project Requirements:
# - Look into p4est for managing the quadtree
# - Create a simple test for the solver with a known analytical solution 
#   to compare to on the region [0,1] x [0,1]
# - Create algorithm for seeing if quad should be refined, based on 
#   function values with chebyshev points
# - Create heatmap display for integral error at different levels of 
#   quadtree depth, generally create multiple methods of display
using P4estTypes;
using Polynomials
using FastChebInterp
using Inti, Meshes, CairoMakie, StaticArrays, Gmsh
using PolygonAlgorithms
using GeometryBasics
using Colorfy
using FMM2D
using IterativeSolvers, LinearAlgebra

@enum GEOS null=0 cut=1 regular=2 contains=3 outside=4

p = 12 # deg of interpolating polys for function refinement
len = 2^30# size of default quadtree
meshsize = 0.025
max_triangle_size = meshsize # temporary constraint
contain_tol = 1e-15

mesh_len = -1
num_boundaries = -1
bndry_mesh_points = []
bndry_mesh_orientations = []


const max_num_boundaries = 16

struct geo_dom_data
    GEO::GEOS
    DOM::NTuple{max_num_boundaries, Int32}
    PARENT_GEO::GEOS
    PARENT_DOM::NTuple{max_num_boundaries, Int32}
end

r0 = (-0.5, 0)
forcing(x) = -2*pi^2 * sin(pi*x[1]) * sin(pi*x[2]) + exp(-600*((x[1]- r0[1])^2 + (x[2] - r0[2])^2))

function analytic_solution(x, y)
    sin(pi*x) * sin(pi*y)
end

tol = 1e-13

function treeToCoords(x, y, bounds, lengths)
    # right now just converts to unit square. change later lol
    bounds[1][1]+x/len*lengths[1], bounds[1][2]+y/len*lengths[2]
end

function ResolvedSourceCheby(Q, f, tol)
    # first, get coefficients of interpolant into c
    coord = P4estTypes.coordinates(Q)
    width = len / 2^P4estTypes.level(Q)
    lb = treeToCoords(coord[1], coord[2], bounds, mesh_len)
    ub = treeToCoords(coord[1] + width, coord[2] + width, bounds, mesh_len)
    x = chebpoints((p, p), lb, ub)
    chebpoly = chebinterp(f.(x), lb, ub, tol=0)
    c = chebpoly.coefs
    c_size = size(c)[1]
    # then, calculate error
    E = 0
    for i  = 1:c_size
        E += abs(c[i, c_size]) + abs(c[c_size, i])
    end
    if E / (2p) < tol
        true
    else
        false
    end
end

function QuadToCoords(quadrant)
    coord = P4estTypes.coordinates(quadrant)
    level = P4estTypes.level(quadrant)
    width = mesh_len / 2^level
    tl = treeToCoords(coord[1], coord[2], bounds, mesh_len)
    tr = (tl[1]+width[1], tl[2])
    bl = (tl[1], tl[2]+width[2])
    br = (tl[1]+width[1], tl[2]+width[2])
    return [tl, tr, br, bl]
end

function classifyQuad(quadrant, boundaries, orientations, delta)
    GEO = regular
    DOM = Set()
    BiQ = [[] for boundary in boundaries]
    QiB = [[] for boundary in boundaries]
    quad = QuadToCoords(quadrant)
    for (i, boundary) in enumerate(boundaries)
        for boundary_point in boundary
            bndry_contains = PolygonAlgorithms.contains(quad, (boundary_point[1], boundary_point[2]), atol=contain_tol)
            push!(BiQ[i], bndry_contains)
        end
        for quad_point in quad
            bndry_contains = PolygonAlgorithms.contains(boundary, quad_point, atol=contain_tol)
            orient_bool = orientations[i] > 0 ? false : true
            # PROBLEM: main outer boundary does not contain (0, 0), this is weird!!
            push!(QiB[i], xor(bndry_contains, orient_bool))
        end
    end
    for (i, boundary) in enumerate(boundaries)
        if !any(QiB[i]) && !all(BiQ[i])
            GEO = outside
            DOM = i
            return GEO, DOM
        end
    end
    # now deal with quadrants that are partially inside bdry
    for (i, boundary) in enumerate(boundaries)
        # if any boundary point is in the quadrant, it is cut
        BiQ_cond = any(BiQ[i])
        # if any quadrant point is not in the boundary, it is cut
        QiB_cond = !all(QiB[i])
        if BiQ_cond || QiB_cond
            GEO = cut
            DOM = i
            if all(BiQ[i])
                GEO = contains
            end
            break
        end
    end
    if P4estTypes.level(quadrant) == 0
        return GEO, DOM
    end
    (parent_GEO, parent_DOM) = getParentData(quadrant)
    if length(parent_DOM) > 0
        if GEO == regular
            DOM = classifyQuadDomainClose(quadrant, boundaries, delta)
        end
    end
    return GEO, DOM
end

function getParentData(quadrant)
    if P4estTypes.level(quadrant) == 0
        error("Level 0 quad is asking for ParentData, ensure starting quad encompasses the whole mesh boundary")
    end
    geo_dom_struct = P4estTypes.unsafe_loaduserdata(quadrant, geo_dom_data)
    GEO = geo_dom_struct.PARENT_GEO
    DOM = Set(geo_dom_struct.PARENT_DOM)
    delete!(DOM, -1)
    return GEO, DOM
end

function classifyQuadDomainClose(quadrant, boundaries, delta)
    DOM = []
    quad = QuadToCoords(quadrant)
    max_quad_len = distance(quad[1], quad[3]) # corner to corner distance
    for (i, boundary) in enumerate(boundaries)
        # going to approximate h gamma locally at each boundary point
        # this is slow but works better for later?

        # first, initial pass for the first quadrant point
        for (j, point) in enumerate(boundary)
            dist = distance(quad[1], point)
            # now calc h
            h = get_h(boundary, j)
            if dist < (1+delta)*h
                push!(DOM, i)
                break
            end
            # if the point is almost close to the quad, check the other corners
            if dist < (1+delta)*h + max_quad_len
                found = false
                for k in 2:4
                    dist = distance(quad[k], point)
                    if dist < (1+delta)*h
                        push!(DOM, i)
                        found = true
                        break
                    end
                end
                if found
                    break
                end
            end
        end
    end
    
    return DOM
end

function distance(a, b)
    return sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)
end

function get_h(boundary, j)
    boundary_len = length(boundary)
    if j > 1 && j < boundary_len
        return (distance(boundary[j-1], boundary[j]) + distance(boundary[j], boundary[j+1]))/2
    elseif j == 1
        return (distance(boundary[boundary_len], boundary[j]) + distance(boundary[j], boundary[j+1]))/2
    elseif j == boundary_len
        return (distance(boundary[j-1], boundary[j]) + distance(boundary[j], boundary[1]))/2
    end
end

function get_GEO_DOM(quadrant)
    geo_dom_struct = P4estTypes.unsafe_loaduserdata(quadrant, geo_dom_data)
    GEO = geo_dom_struct.GEO
    DOM = Set(geo_dom_struct.DOM)
    delete!(DOM, -1)
    return GEO, DOM
end

function getQuadCoordsFromInternalCoord(forest, coordinate)
    for quad in forest[1]
        quad_coords = QuadToCoords(quad)
        if PolygonAlgorithms.contains(quad_coords, coordinate, atol=contain_tol)
            return P4estTypes.coordinates(quad), P4estTypes.level(quad)
        end
    end
end

# !!! ADD TO THIS LATER AS WE DEAL WITH HARDER BOUNDARIES !!! #
function refinementStage(forest, treeid, quadrant)
    delta = 0.1
    eps = 0.3
    GEO, DOM = get_GEO_DOM(quadrant)
    quad_level = P4estTypes.level(quadrant)
    quad_coords = QuadToCoords(quadrant)
    formatted_coords = [quad_coords[1][1], quad_coords[3][1], quad_coords[1][2], quad_coords[3][2]]
    if P4estTypes.coordinates(quadrant) == (167772160, 528482304) && P4estTypes.level(quadrant) == 7
        println("HI")
    end
    if GEO == regular
        if ~ResolvedSourceCheby(quadrant, forcing, tol)
            return true
        end
        if length(DOM) == 0
            return false
        end
    elseif GEO == outside
        return false
    elseif GEO == contains
        return true
    end
    # now the quad is either cut or regular and close to a boundary
    # so we check further
    for j in DOM
        boundary = bndry_mesh_points[j]
        # loop through close boundaries
        for (i, point) in enumerate(boundary)
            dist = minimum(distance(point, quad_coords[k]) for k in 1:4)
            local_h = get_h(boundary, i)
            quad_h = distance(quad_coords[1], quad_coords[2])
            if quad_h > (1+eps) * local_h
                if GEO == regular && dist < (1+delta)*local_h
                    return true
                elseif GEO != regular
                    return true
                end
            end
        end
    end
    return false
end

function coarsenStage(forest, treeid, siblings)
    # first, we must get the GEO value of each of the siblings
    for quadrant in siblings
        GEO, DOM = get_GEO_DOM(quadrant)
        if GEO == regular
            return false
        end
    end
    return true
end

function initializeQuadrant(forest, treeid, quadrant)
    # classifyQuad should NEVER RUN if quadrant doesn't have parent data defined!!!
    GEO, DOM = classifyQuad(quadrant, bndry_mesh_points, bndry_mesh_orientations, 0.1)
    # parent is (939524096, 536870912) level 4
    # child is (973078528, 536870912) level 5
    DOM_list = [x for x in DOM] 
    DOM_list = vcat(DOM_list, [-1 for i in 1:max_num_boundaries-length(DOM_list)])
    null_list = NTuple{max_num_boundaries, Int32}([-1 for i in 1:max_num_boundaries])
    DOM_NTuple = NTuple{max_num_boundaries, Int32}(DOM_list)
    geo_dom_val = geo_dom_data(GEO, DOM_NTuple, null, null_list)
    P4estTypes.unsafe_storeuserdata!(quadrant, geo_dom_val)
end

function replaceQuadrant(forest, treeid, outgoing, incoming)
    if length(outgoing) == 1 && length(incoming) > length(outgoing)
        parent_GEO, parent_DOM = get_GEO_DOM(outgoing[1])
        DOM_list = [x for x in parent_DOM] 
        DOM_list = vcat(DOM_list, [-1 for i in 1:max_num_boundaries-length(DOM_list)])
        null_list = NTuple{max_num_boundaries, Int32}([-1 for i in 1:max_num_boundaries])
        parent_DOM_NTuple = NTuple{max_num_boundaries, Int32}(DOM_list)
        for quadrant in incoming
            # geo_dom_struct = P4estTypes.unsafe_loaduserdata(quadrant, geo_dom_data)
            # GEO = geo_dom_struct.GEO
            # DOM_NTuple = geo_dom_struct.DOM
            # if P4estTypes.coordinates(quadrant) == (973078528, 536870912)
            #     println("-------------------")
            #     println(P4estTypes.level(quadrant))
            #     println(GEO)
            #     println(DOM_NTuple)
            #     println(parent_GEO)
            #     println(parent_DOM_NTuple)
            #     println(P4estTypes.coordinates(outgoing[1]))
            #     println(P4estTypes.level(outgoing[1]))
            # end
            geo_dom_val_tmp = geo_dom_data(null, null_list, parent_GEO, parent_DOM_NTuple)
            P4estTypes.unsafe_storeuserdata!(quadrant, geo_dom_val_tmp)
            GEO, DOM = classifyQuad(quadrant, bndry_mesh_points, bndry_mesh_orientations, 0.1)
            DOM_list = [x for x in DOM] 
            DOM_list = vcat(DOM_list, [-1 for i in 1:max_num_boundaries-length(DOM_list)])
            DOM_NTuple = NTuple{max_num_boundaries, Int32}(DOM_list)
            geo_dom_val = geo_dom_data(GEO, DOM_NTuple, parent_GEO, parent_DOM_NTuple)
            P4estTypes.unsafe_storeuserdata!(quadrant, geo_dom_val)
            # println(geo_dom_val)
        end
    end
end

function getInternalBoundaryPoints(tree)
    # first, construct polygon to represent all quadrants in quadtree
    quadtree_points = map(QuadToCoords, tree)
    quadtree_poly = PolygonAlgorithms.union_geometry(quadtree_points...)
    remove_points = [[] for i in 1:num_boundaries]
    for quadrant in tree
        GEO, DOM = get_GEO_DOM(quadrant)
        for boundary in 1:num_boundaries
            if ((GEO == regular && boundary in DOM) || (GEO == cut)) && boundary == minimum(DOM)
            # if ((GEO == regular && length(DOM)>0))
                # coords = QuadToCoords(quadrant)
                # if all(coords[i][1] < 0 && coords[i][2] < 0 && coords[i][1]^2 + coords[i][2]^2 < 1 for i in 1:4)
                #     println(quadrant)
                # end
                push!(remove_points[boundary], QuadToCoords(quadrant))
            end
        end
    end
    internalBoundaryPoints = [PolygonAlgorithms.union_geometry(remove_points[i]...) for i in 1:num_boundaries if length(remove_points[i]) > 0]

    if length(internalBoundaryPoints) == 0
        error("No internal boundaries found, try refining more!")
    end

    return internalBoundaryPoints
end

function createMesh()
    # this is good so far. I can now create inscribed mesh polygons from a
    # list of points. Very good, now will work on refinement criteria
    # and return once proper internal mesh is done.
    # Creates and returns a gmsh mesh to use
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 3) # turn to 4/5 for info or debug, 3 is all that is necessary I think though

    circle = Inti.gmsh_curve(0, 2pi;meshsize) do s
        return Inti.Point2D(cos(s), sin(s))
    end
    curve_loop = gmsh.model.occ.addCurveLoop([circle])

    # inscribedPoints = [(0, 0), (0.5, 0), (0.6, 0.6), (0, 0.5), (-0.05, 0.25), (-0.025, -0.5), (0, 0)]
    inscribedPoints = [(0, 0), (0.5, 0.05), (0, 0.1), (-0.05, 0.4), (-0.1, 0.1), (-0.4, 0.05), (-0.1, 0), (-0.05, -0.4), (0, 0)]

    point_tags = []
    for point in inscribedPoints
        tag = gmsh.model.occ.addPoint(point[1], point[2], 0, max_triangle_size)
        push!(point_tags, tag)
    end
    curve_tags = []
    for i in 1:length(point_tags)-1
        tag = gmsh.model.occ.addLine(point_tags[i], point_tags[i+1])
        push!(curve_tags, tag)
    end
    inner_loop = gmsh.model.occ.addCurveLoop(curve_tags)
    surf = gmsh.model.occ.addPlaneSurface([curve_loop, inner_loop])

    
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", max_triangle_size)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(2)
    gmsh.write("testing.msh")
    gmsh.finalize()
end

function showMesh(msh)
    println("Showing mesh!!!")
    Ω = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, msh)
    Γ = Inti.boundary(Ω)
    Ω_msh = @views msh[Ω]
    Γ_msh = @views msh[Γ]
    fig = viz(
        Ω_msh;
        segmentsize = 1,
        showsegments = true,
        axis = (aspect = DataAspect(),),
        figure = (; size = (500, 400)),
    )
    viz!(Γ_msh; color = :red, segmentsize = 1)
    display(fig)
end

function createMeshByGeoFromQuadtree(forest)
    # first is regular + empty dom
    # second is regular + elements in dom
    # third is cut
    # fourth is contains
    # fifth is outside
    # sixth is null
    box_lists = [[] for i in 1:(5 + num_boundaries)]
    for quad in forest[1]
        # println(quad)
        coords = QuadToCoords(quad)
        GEO, DOM = get_GEO_DOM(quad)
        box = Meshes.Box(coords[1], coords[3])
        if GEO == regular && length(DOM) == 0
            push!(box_lists[1], box)
        elseif GEO == cut
            push!(box_lists[2], box)
        elseif GEO == contains
            push!(box_lists[3], box)
        elseif GEO == outside
            push!(box_lists[4], box)
        elseif GEO == null
            push!(box_lists[5], box)
        elseif GEO == regular && length(DOM) > 0
            push!(box_lists[5+minimum(DOM)], box)
        end
    end
    geometries = []
    for (i, box_list) in enumerate(box_lists)
        if length(box_list) > 0
            push!(geometries, Meshes.GeometrySet(box_list))
        end
    end
    return geometries
end

function showGeoQuadtreeMesh(msh, forest)
    geometries = createMeshByGeoFromQuadtree(forest)
    Ω = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, msh)
    Γ = Inti.boundary(Ω)
    # Ω_msh = @views msh[Ω]
    Γ_msh = @views msh[Γ]
    fig = viz(
        Γ_msh;
        segmentsize = 1,
        showsegments = true,
        color= :red,
        axis = (aspect = DataAspect(),),
        figure = (; size = (500, 400)),
    )
    # print(geometries)
    # viz!(Γ_msh; color = :red, segmentsize = 1)
    color_vals = [0.8, 0.5, 0, 0.2, 1]
    for i = 1:num_boundaries
        push!(color_vals, 0.8 - i*(0.2)/num_boundaries)
    end
    colorfier = Colorfy.Colorfier(color_vals)
    colors = Colorfy.colors(colorfier)
    for i in 1:length(geometries)
        println(i)
        viz!(geometries[i]; color=colors[i], showsegments = true, segmentsize = 1)
        # viz!(geometries[i]; color=i, showsegments = true, segmentsize = 1)
    end

    display(fig)
end

function createMeshFromQuadtree(forest)
    box_list = []
    for quad in forest[1]
        coords = QuadToCoords(quad)
        # println(coords)
        GEO, DOM = get_GEO_DOM(quad)
        if GEO == regular && length(DOM) == 0
        # if (GEO == regular && 2 in DOM) || (GEO == cut && 2 in DOM)
            box = Meshes.Box(coords[1], coords[3])
            push!(box_list, box)
        end
    end
    geometry = Meshes.GeometrySet(box_list)
    return geometry
end

function getBoundaryMesh(msh)
    # returns the boundary of the mesh
    domain = Inti.Domain(e -> Inti.geometric_dimension(e)==2, msh)
    bndry_mesh = view(msh, Inti.boundary(domain))
    bndry_mesh
end

function getBoundaryPoints(bndry_mesh)
    verts = Vector{Vector{Tuple{Float64, Float64}}}()
    orientations = []
    function compute(iter, connectivity, orientation)
        compute_verts = Vector{Vector{Tuple{Float64, Float64}}}()
        push!(compute_verts, Vector{Tuple{Float64, Float64}}())
        compute_orientations = []
        j = 1
        original_connectivity = -1
        for (i, el) in enumerate(iter)
            v1 = Inti.vertices(el)[1]
            # boundary idx indices, change if using non-triangular mesh
            c1 = connectivity[1, i]
            c2 = connectivity[3, i]
            if original_connectivity == -1
                # triggers on first index of new boundary, so add orientation
                or = orientation[i]
                push!(compute_orientations, or)
                original_connectivity = c1
            end            
            # now, to fix the issue with point containment around y=0, we must
            # lock points close to y=0 directly to y=0.
            if abs(v1[2]) < 1e-15
                push!(compute_verts[j], (v1[1], 0.0))
            else
                push!(compute_verts[j], (v1[1], v1[2]))
            end
            # now decide if placing in verts or moving to next
            if (c2 == original_connectivity)
                j += 1
                original_connectivity = -1
                push!(compute_verts, [])
            end
        end
        return compute_verts, compute_orientations
    end
    for (i, E) in enumerate(Inti.element_types(bndry_mesh))
        connectivity = Inti.connectivity(bndry_mesh, E)
        orientation = Inti.orientation(bndry_mesh, E)
        tmp_verts, tmp_orientations = compute(Inti.elements(bndry_mesh, E), connectivity, orientation)
        append!(verts, tmp_verts)
        append!(orientations, tmp_orientations)
    end
    verts[1:end-1], orientations
end

function getMeshBounds(msh)
    # gets the bounds of a gmsh mesh (assumed to be 2d)
    # and returns the lower/upper bound for each dimension

    # first get boundary
    bndry_mesh = getBoundaryMesh(msh)
    function compute!(iter, lb, ub)
        for el in iter
            # el is of type LagrangeElement (an approximating polynomial, representing a line)
            verts = Inti.vertices(el)
            for i in 1:2
                for dim in 1:2
                    if verts[i][dim] < lb[dim]
                        lb[dim] = verts[i][dim]
                    end 
                    if verts[i][dim] > ub[dim]
                        ub[dim] = verts[i][dim]
                    end
                end
            end
        end
    end
    lb = [Inf, Inf]
    ub = [-Inf, -Inf]
    for E in Inti.element_types(bndry_mesh)
        compute!(Inti.elements(bndry_mesh, E), lb, ub)
    end
    
    lb, ub
end

function createCurveLoop(points)
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
    return gmsh.model.occ.addCurveLoop(lineTags)
end

function createMeshWithInternalBoundary(boundary_points, internal_points)
    # assumes boundary_points aligns with internal_points, which it had better by default. If it doesn't
    # I can totally create a shitty band aid solution but I don't want to do that
    # also assumes first set of boundary points is exterior
    # meshes = []
    internal_curves = []
    external_curves = []
    for i in 1:num_boundaries
        boundary_curve_loop = createCurveLoop(boundary_points[i])
        internal_curve_loop = createCurveLoop(internal_points[i])
        if bndry_mesh_orientations[i] == 1
            surf = gmsh.model.occ.addPlaneSurface([boundary_curve_loop, internal_curve_loop])
        else
            surf = gmsh.model.occ.addPlaneSurface([internal_curve_loop, boundary_curve_loop])
        end
        push!(internal_curves, internal_curve_loop)
        push!(external_curves, boundary_curve_loop)
    end
    return internal_curves, external_curves
end

function createQuadtree(maxlevel = -1)
    conn = P4estTypes.Connectivity{4}(:unitsquare)
    forest = pxest(conn; min_level=0, init_function = initializeQuadrant, data_type = geo_dom_data)
    refine!(forest; refine = refinementStage, replace = replaceQuadrant, recursive = true, maxlevel = maxlevel)
    coarsen!(forest; coarsen = coarsenStage, init = initializeQuadrant, recursive = true)
    balance!(forest; init = initializeQuadrant)
    return forest
end

function initializeMeshVariables(mesh)
    global bounds = getMeshBounds(mesh)

    mesh_len = [bounds[2][1] - bounds[1][1], bounds[2][2] - bounds[1][2]]
    for i in 1:2
        bounds[1][i] -= mesh_len[i] * 0.05
        bounds[2][i] += mesh_len[i] * 0.05
    end
    global mesh_len = [bounds[2][1] - bounds[1][1], bounds[2][2] - bounds[1][2]]
    global original_mesh = mesh
    global bndry_mesh = getBoundaryMesh(mesh)
    global bndry_mesh_points, bndry_mesh_orientations = getBoundaryPoints(bndry_mesh)
    global num_boundaries = length(bndry_mesh_points)

    if num_boundaries > max_num_boundaries
        error("Too many boundaries! Should be < 256")
    end
end

function filterBoundaries(internalBoundaryPoints)
    internalBoundaries = []

    for boundary in 1:num_boundaries
        for j in 1:length(internalBoundaryPoints[boundary])
            point = internalBoundaryPoints[boundary][j][1]
            point_contains = PolygonAlgorithms.contains(bndry_mesh_points[boundary], point) 
            if bndry_mesh_orientations[boundary] == 1
                orientation = true
            else
                orientation = false
            end
            point_within = xor(point_contains, orientation)
            if ~point_within
                # then second boundary side is within
                push!(internalBoundaries, internalBoundaryPoints[boundary][j])
            end
        end
    end
    
    return internalBoundaries
end

function cullBoundaryPoints(internalBoundaries)
    culledInternalBoundaries = []

    for (j, boundary) in enumerate(internalBoundaries)
        push!(culledInternalBoundaries, [])
        for i in 1:length(boundary)
            if distance(boundary[i], boundary[mod1(i+1, length(boundary))]) >= 1e-8
                # if i != length(boundary) || j != 2
                #     push!(culledInternalBoundaries[j], boundary[i])
                # end
                # really weird issues happen when the max refinement is 4 here. gonna focus on this later but this is pretty problematic
                modifiedPoint = (boundary[i][1] + 1e-7*(rand()-0.5), boundary[i][2] + 1e-7*(rand()-0.5))
                push!(culledInternalBoundaries[j], modifiedPoint)
            end
        end
    end

    return culledInternalBoundaries
end

function createCombinedMesh(forest, surface_tag)
    element_type = 3
    # gmsh.model.add("quadtreeRegions")
    numQuads = 0
    node_coords = []
    for quadrant in forest[1]
        quad_coords = QuadToCoords(quadrant)
        geo, dom = get_GEO_DOM(quadrant)
        if geo == regular && length(dom) == 0
            numQuads += 1
            for coord in quad_coords
                push!(node_coords, coord[1], coord[2], 0.0)
            end
        end
    end
    node_tags = [x for x in 1:numQuads*4]
    element_tags = [x for x in 1:numQuads]

    gmsh.model.mesh.addNodes(2, surface_tag, node_tags, node_coords)
    gmsh.model.mesh.addElementsByType(surface_tag, element_type, element_tags, node_tags)
end

function createQuadtreeMesh()
    # Creates the meshes describing the boundary region and the quadtree
    gmsh.open("testing.msh")
    mesh = Inti.import_mesh(; dim = 2)

    # showMesh(mesh)

    initializeMeshVariables(mesh)
    
    forest = createQuadtree()

    internalBoundaryPoints = getInternalBoundaryPoints(forest[1])
    
    internalBoundaries = filterBoundaries(internalBoundaryPoints)
        
    culledInternalBoundaries = cullBoundaryPoints(internalBoundaries)

    gmsh.model.add("BoundaryRegions")

    internal_curves, external_curves = createMeshWithInternalBoundary(bndry_mesh_points, culledInternalBoundaries)

    gmsh.model.occ.synchronize()
    surface_tag = gmsh.model.addDiscreteEntity(2)

    createCombinedMesh(forest, surface_tag)

    gmsh.model.occ.synchronize()

    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.option.setNumber("Mesh.MeshSizeMax", max_triangle_size)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(2)
    gmsh.write("combined.msh")
    combined_msh = Inti.import_mesh(; dim = 2)
    gmsh.finalize()

    return combined_msh
end

function testArea(domain_quad)
    expected_area = pi - 0.085 # pi is circle area, 0.085 is star area
    quadeval_f = (q) -> 1
    computed_area = Inti.integrate(quadeval_f, domain_quad)
    print("Difference between areas: ")
    println(abs(expected_area - computed_area))

    # note so that this doesn't have to be run more: with meshsize = 0.01, 4651 triangles and 50874 quadrature nodes,

    # hardcoding values bc they take a while to compute
    meshsizes = [0.08, 0.04, 0.02, 0.01]
    area_diffs = [3.31e-3, 8.28e-4, 2.08e-4, 5.24e-5]
    p = 0
    for i in 1:length(area_diffs)-1
        p += log(area_diffs[i+1] / area_diffs[i]) / log(1/2)
    end
    p /= length(area_diffs) - 1
    println("Area converges at rate p=", p) # quadratic convergence! small sample size but still good to know
end

function testBoundaryLength(boundary_quad)
    star_len = sqrt(0.3^2 + 0.05^2)*4 + sqrt(0.4^2 + 0.05^2)*2 + sqrt(0.5^2 + 0.05^2)*2
    expected_length = 2*pi + star_len
    quadeval_f = (q) -> 1
    computed_length = Inti.integrate(quadeval_f, boundary_quad)
    # qorder 1 - 121, qorder 2,3 - 242, qorder 4,5 - 363, qorder 6 - 484 (num quadrature points)
    # convergence also slows extremely quickly (once it gets down to around 4e-7, close to tol though so should be ok? also maybe expected)
    # convergence is much more exact (for this at least) when reducing meshsize rather than qorder. Makes sense as the qorder doesn't do
    # anything at a point, they just further approximate a bad line (which converges to 2pi when meshsize is reduced instead).
    println(expected_length)
    println(computed_length)
    print("Difference between lengths: ")
    println(abs(expected_length - computed_length))
end

function testGreenThirdIdentity(domain_quad, boundary_quad)
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

    boundary_integral = Inti.integrate(boundary_function, boundary_quad)
    domain_integral = Inti.integrate(domain_function, domain_quad)
    calculated_u_val = domain_integral - boundary_integral

    println(Inti.integrate(boundary_function, boundary_quad))
    println(Inti.integrate(domain_function, domain_quad))
    println(expected_u_val)
    println(calculated_u_val)
    # changing qorder changes very little; meshsize is the thing to change here to get convergence
    # but the convergence feels much much slower than it should be right now. We shouldn't have to go down to 0.025 size triangles
    # with 23000 domain quadrature nodes and 1500 boundary quadrature nodes to get an error of 1e-3.
    # but with 11000 domain quadrature nodes we get an equivalent error
    # and with 3300 domain quadrature nodes we only have an error of 2e-2. 
    return
end

function calculateMeshvals(mesh)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    boundary_dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 1, bndry_mesh)
    # boundary = Inti.external_boundary(dom)
    dom_mesh = Inti.view(mesh, dom)
    boundary_mesh = Inti.view(original_mesh, boundary_dom)

    domain_quad = Inti.Quadrature(dom_mesh; qorder=4)
    boundary_quad = Inti.Quadrature(boundary_mesh; qorder=6)
    
    println(domain_quad)
    println(boundary_quad)

    # greenfn = (x) -> 1/(2pi) * log(distance(x, r0))
    
    # testArea(domain_quad)
    # testBoundaryLength(boundary_quad)
    testGreenThirdIdentity(domain_quad, boundary_quad)

    # test area with Quadrature GOOD
    # test length with Quadrature GOOD-ish
    # greens third identity test GOOD-ish
end

# driver code
gmsh.initialize()
gmsh.option.setNumber("General.Verbosity", 3) # turn to 4/5 for info or debug, 3 is all that is necessary I think though
msh = createQuadtreeMesh()
showMesh(msh)
calculateMeshvals(msh)



# getElementQualities
# visualize hybrid mesh

# PROBLEM WITH MESHING:
# Some quads are associated with more than one boundary. Which should it mesh for?
# temporarily going to choose the mesh with lower value

# create inti mesh from forest, msh

# TODO
# - fix issues when boundaries are close
# - SPEED UP PLEASE
# - change small triangles system
#    - right now we just cap the size of triangle formed for every boundary to meshsize. While the external boundary should have
#      meshsize-sized triangles, we shouldn't require every boundary to be like this, especially inner boundaries
# - fix issues with boundaries sometimes ignoring quads that are closer than other quads that get cut (prob issue with quad classification)

# TO DISCUSS
# - fixed boundary strip length to be relative to meshsize
# - area and length testing
# - qorder vs meshsize convergence
# - greens third identity tests (and why its still bad)