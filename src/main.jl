# Project Requirements:
# - Look into p4est for managing the quadtree
# - Create a simple test for the solver with a known analytical solution 
#   to compare to on the region [0,1] x [0,1]
# - Create algorithm for seeing if quad should be refined, based on 
#   function values with chebyshev points
# - Create heatmap display for integral error at different levels of 
#   quadtree depth, generally create multiple methods of display

@enum GEOS null=0 cut=1 regular=2 contains=3 outside=4

p = 12 # deg of interpolating polys for function refinement
len = 2^30# size of default quadtree
contain_tol = 1e-15
curvature_size = 30

mesh_len = -1.0
num_boundaries = -1
bndry_mesh_points = []
bndry_mesh_orientations = []
bndry_mesh_splines = []
global delta = 0.1
global eps = 0.25
global const forcing = Ref{Function}((x) -> 0)


const max_num_boundaries = 16

struct geo_dom_data
    GEO::GEOS
    DOM::NTuple{max_num_boundaries, Int32}
    PARENT_GEO::GEOS
    PARENT_DOM::NTuple{max_num_boundaries, Int32}
end

tol = 1e-13

function treeToCoords(x, y)
    # right now just converts to unit square. change later lol
    bounds[1][1]+x/len*mesh_len, bounds[1][2]+y/len*mesh_len
end

function ResolvedSourceCheby(Q, f, tol)
    # first, get coefficients of interpolant into c
    coord = P4estTypes.coordinates(Q)
    width = len / 2^P4estTypes.level(Q)
    lb = treeToCoords(coord[1], coord[2])
    ub = treeToCoords(coord[1] + width, coord[2] + width)
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
    tl = treeToCoords(coord[1], coord[2])
    tr = (tl[1]+width, tl[2])
    bl = (tl[1], tl[2]+width)
    br = (tl[1]+width, tl[2]+width)
    return [tl, tr, br, bl]
end

function classifyQuad(quadrant, boundaries, orientations)
    GEO = regular
    DOM = Set()
    BiQ = [[] for boundary in boundaries]
    QiB = [[] for boundary in boundaries]
    quad = QuadToCoords(quadrant)
    for (i, boundary) in enumerate(boundaries)
        for boundary_point in boundary
            bndry_contains = PolygonAlgorithms.contains(quad, boundary_point, atol=0.0)
            push!(BiQ[i], bndry_contains)
        end
        for quad_point in quad
            bndry_contains = PolygonAlgorithms.contains(boundary, quad_point, atol=0.0)
            orient_bool = orientations[i] > 0 ? false : true
            push!(QiB[i], xor(bndry_contains, orient_bool))
        end
    end
    for (i, boundary) in enumerate(boundaries)
        # if no quad points are in the boundary and some boundary point is not in the quad
        if !any(QiB[i])
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
            DOM = classifyQuadDomainClose(quadrant, boundaries)
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

function classifyQuadDomainClose(quadrant, boundaries)
    DOM = []
    quad = QuadToCoords(quadrant)
    max_quad_len = distance(quad[1], quad[3]) # corner to corner distance
    for (i, boundary) in enumerate(boundaries)
        # first, initial pass for the first quadrant point
        for (j, point) in enumerate(boundary)
            dist = distance(quad[1], point)
            # now calc h
            h = get_h(boundary, j)
            # h = 0.19603428065912154
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
                        # if P4estTypes.coordinates(quadrant) == (536870912, 805306368) && P4estTypes.level(quadrant) == 4
                        #     println(quadrant)
                        #     println((1+delta)*h)
                        #     println(dist)
                        #     println(point)
                        #     println(quad)
                        #     println(k)
                        # end
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

# !!! ADD TO THIS LATER AS WE DEAL WITH HARDER BOUNDARIES !!! #
function refinementStage(forest, treeid, quadrant)
    GEO, DOM = get_GEO_DOM(quadrant)
    quad_level = P4estTypes.level(quadrant)
    quad_coords = QuadToCoords(quadrant)
    formatted_coords = [quad_coords[1][1], quad_coords[3][1], quad_coords[1][2], quad_coords[3][2]]
    if GEO == regular
        if ~ResolvedSourceCheby(quadrant, forcing[], tol)
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
            max_quad_h = distance(quad_coords[1], quad_coords[3])
            if dist < (1+delta)*local_h + max_quad_h
                if quad_h > (1+eps) * local_h
                    if GEO == regular && dist < (1+delta)*local_h
                        return true
                    elseif GEO != regular
                        return true
                    end
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
    GEO, DOM = classifyQuad(quadrant, bndry_mesh_points, bndry_mesh_orientations)
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
            geo_dom_val_tmp = geo_dom_data(null, null_list, parent_GEO, parent_DOM_NTuple)
            P4estTypes.unsafe_storeuserdata!(quadrant, geo_dom_val_tmp)
            GEO, DOM = classifyQuad(quadrant, bndry_mesh_points, bndry_mesh_orientations)
            DOM_list = [x for x in DOM] 
            DOM_list = vcat(DOM_list, [-1 for i in 1:max_num_boundaries-length(DOM_list)])
            DOM_NTuple = NTuple{max_num_boundaries, Int32}(DOM_list)
            geo_dom_val = geo_dom_data(GEO, DOM_NTuple, parent_GEO, parent_DOM_NTuple)
            P4estTypes.unsafe_storeuserdata!(quadrant, geo_dom_val)
        end
    elseif length(incoming) == 1 && length(incoming) < length(outgoing)
        GEO, DOM = classifyQuad(incoming[1], bndry_mesh_points, bndry_mesh_orientations)
        DOM_list = [x for x in DOM] 
        DOM_list = vcat(DOM_list, [-1 for i in 1:max_num_boundaries-length(DOM_list)])
        DOM_NTuple = NTuple{max_num_boundaries, Int32}(DOM_list)
        # now just associate parent_DOM with all boundaries of the outgoing doms
        # TO BE CLEAR: THIS IS NOT FINISHED! NEED BETTER SYSTEM FOR GEO ASSOCIATION
        # BUT DOM IS THE ONLY THING NEEDED NOW SO ONLY DOING THAT
        parent_DOM = Set{Int32}()
        for quadrant in outgoing
            child_GEO, child_DOM = get_GEO_DOM(quadrant)
            union!(parent_DOM, child_DOM)
        end
        parent_DOM_list = [x for x in parent_DOM] 
        parent_DOM_list = vcat(parent_DOM_list, [-1 for i in 1:max_num_boundaries-length(parent_DOM_list)])
        parent_DOM_NTuple = NTuple{max_num_boundaries, Int32}(parent_DOM_list)
        geo_dom_val = geo_dom_data(GEO, DOM_NTuple, null, parent_DOM_NTuple)
        P4estTypes.unsafe_storeuserdata!(incoming[1], geo_dom_val)
    end
end

function getInternalBoundaryPoints(tree)
    # new algorithm idea:
    # 1. get all edges of quads in the boundary strip, add to a set (won't be too many)
    # 2. filter such that we only have edges of quads that are also in the regular region of the quadtree
    # 3. 
    
    function roundEdge(p1, p2)::Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}
        r_p1 = (round(p1[1], digits=12)+0.0, round(p1[2], digits=12)+0.0)
        r_p2 = (round(p2[1], digits=12)+0.0, round(p2[2], digits=12)+0.0)
        return (r_p1, r_p2)
    end

    function getPossibleBoundaryEdges(coords)
        short_edges = Vector{Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}}()
        long_edges = Vector{Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}}()
        # only reversed edges are possible for these boundaries
        for i in 1:4
            p1 = coords[i]
            p2 = coords[mod1(i+1, 4)]
            x_len = p2[1] - p1[1]
            y_len = p2[2] - p1[2]
            midpoint = ((p1[1]+p2[1])/2.0, (p1[2]+p2[2])/2.0)
            left_outside = (p1[1]-x_len, p1[2]-y_len)
            right_outside = (p2[1]+x_len, p2[2]+y_len)
            # note: could be errors based on floating point approx here
            tmp_short_edges = [roundEdge(p2, p1), roundEdge(p2, midpoint), roundEdge(midpoint, p1)]
            tmp_long_edges = [roundEdge(right_outside, p1), roundEdge(p2, left_outside)]
            # println(tmp_edges)
            # tmp_edges = [roundEdge(p2, p1), roundEdge(p2, midpoint), roundEdge(midpoint, p1)]
            append!(short_edges, tmp_short_edges)
            append!(long_edges, tmp_long_edges)
        end
        return short_edges, long_edges
    end

    function getPossibleDomainEdges(coords)
        edges = Vector{Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}}()
        # only reversed edges are possible for these boundaries
        for i in 1:4
            p1 = coords[i]
            p2 = coords[mod1(i+1, 4)]
            midpoint = ((p1[1]+p2[1])/2.0, (p1[2]+p2[2])/2.0)
            tmp_edges = [roundEdge(p1, p2), roundEdge(p1, midpoint), roundEdge(midpoint, p2)]
            append!(edges, tmp_edges)
        end
        return edges
    end

    function addNextEdge(good_edges, edge_connections; remove_elements = true)::Bool
        # returns the next point
        found_edge = false
        for edge in good_edges
            if distance(edge[1], edge_connections[end]) < 1e-10 && distance(edge[2], edge_connections[end-1]) >= 1e-10
                push!(edge_connections, edge[2])
                if remove_elements
                    pop!(good_edges, edge)
                end
                found_edge = true
                break
            elseif distance(edge[2], edge_connections[end]) < 1e-10 && distance(edge[1], edge_connections[end-1]) >= 1e-10
                push!(edge_connections, edge[1])
                if remove_elements
                    pop!(good_edges, edge)
                end
                found_edge = true
                break
            end
        end
        return found_edge
    end
    
    function getFirstPoints(good_edges)
        start_edge = first(good_edges)
        edge_connections = [start_edge[2], start_edge[1]] # so this goes in reverse now
        found_edge = addNextEdge(good_edges, edge_connections; remove_elements = false)
        while found_edge
            found_edge = addNextEdge(good_edges, edge_connections; remove_elements = false)
            if edge_connections[1] == edge_connections[end]
                # then this is a loop, any point works, just return something
                break 
            end
        end
        # then final element will be the start
        return edge_connections[end], edge_connections[end-1]
    end

    edge_list = Vector{Vector{Tuple{Float64, Float64}}}()
    # first pass, get regular edges
    regular_edges = Set{Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}}()
    for quad in tree
        GEO, DOM = get_GEO_DOM(quad)
        if GEO == regular && length(DOM) == 0
            coords = QuadToCoords(quad)
            for i in 1:4
                p1 = coords[i]
                p2 = coords[mod1(i+1, 4)]
                edge = roundEdge(p1, p2)
                # rev_edge = roundEdge(p2, p1)
                push!(regular_edges, edge)
                # push!(regular_edges, rev_edge)
            end
        end
    end
    edge_list_associations = []
    for boundary in 1:num_boundaries
        # second pass, filter by boundary info
        good_edges = Set{Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}}()
        for quad in tree
            GEO, DOM = get_GEO_DOM(quad)
            if (GEO == regular && boundary in DOM || GEO == cut && boundary in DOM)
            # if GEO == regular && length(boundary) > 0 || GEO == cut
                coords = QuadToCoords(quad)
                short_edges, long_edges = getPossibleBoundaryEdges(coords)
                skip_tag = false
                for edge in short_edges
                    if edge in regular_edges
                        push!(good_edges, edge)
                    end
                end
                for (j, edge) in enumerate(long_edges)
                    side = Int64(ceil(j/2))-1
                    if edge in regular_edges
                        push!(good_edges, short_edges[3*side+1])
                    end
                end
            end
        end
        # final pass, reconnect large edges. should be very quick
        for quad in tree
            GEO, DOM = get_GEO_DOM(quad)
            if GEO == regular && length(DOM) == 0
                coords = QuadToCoords(quad)
                edges = getPossibleDomainEdges(coords)
                for j in 0:3
                    if edges[3*j+2] in good_edges && edges[3*j+3] in good_edges
                        pop!(good_edges, edges[3*j+2])
                        pop!(good_edges, edges[3*j+3])
                        push!(good_edges, edges[3*j+1])
                    end
                end
            end
        end
        # then connect edges
        
        # ts is too long!!
        n = length(good_edges)
        # println(n)
        edge_connection_list = Vector{Vector{Tuple{Float64, Float64}}}()
        while length(good_edges) > 0
            og_len = length(good_edges)
            k = 0
            # println(getFirstPoints(good_edges))
            edge_connections = [getFirstPoints(good_edges)...]
            if Tuple(edge_connections) in good_edges
                # println(length(good_edges))
                pop!(good_edges, Tuple(edge_connections))
                # println(length(good_edges))
            end

            found_edge = true
            while found_edge
                # println("----")
                # println(edge_connections)
                found_edge = addNextEdge(good_edges, edge_connections)
                k += 1
                if k > n
                    println(edge_connections)
                    error("stop infinite loop")
                end
            end
            push!(edge_connection_list, edge_connections)
            if length(good_edges) == og_len
                println(length(good_edges))
                error("no change in length")
            end
        end
        append!(edge_list_associations, [boundary for x in 1:length(edge_connection_list)])
        append!(edge_list, edge_connection_list)
    end
    # now associate edge sets with each other if they are open
    final_edge_list_associations = []
    final_edge_list = []
    for (i, edge_set) in enumerate(edge_list)
        if distance(edge_set[1], edge_set[end]) <= 1e-10
            push!(final_edge_list_associations, [edge_list_associations[i]])
            push!(final_edge_list, edge_set[1:end-1])
            continue
        end
        # otherwise, set must be open

        # suppose separated quadtree. then we have three edge sets, each associated with the same boundary
        # suppose hole close to quadtree. Then we have two edge sets, one associated with each boundary, that should be combined into one.
        found_before = false
        for association_set in final_edge_list_associations
            if i in association_set && length(association_set) > 1
                found_before = true
            end
        end
        if found_before
            continue
        end
        edge_list_association = [i]
        for (j, other_edge_set) in enumerate(edge_list)
            if distance(edge_set[end], other_edge_set[1]) <= 1e-10
                push!(edge_list_association, edge_list_associations[j])
                append!(edge_set, other_edge_set)
                if distance(edge_set[1], edge_set[end]) <= 1e-10
                    break
                end
            end
        end
        push!(final_edge_list_associations, edge_list_association)
        push!(final_edge_list, edge_set[1:end-1])
    end
    # for edge_set in final_edge_list
    #     showEdgeList(edge_set)
    # end
    return final_edge_list_associations, final_edge_list
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
        else
            push!(geometries, [])
        end
    end
    return geometries
end

function showGeoQuadtreeMesh(meshes, forest)
    println("Showing GEO Quadtree mesh!")
    geometries = createMeshByGeoFromQuadtree(forest)
    boundary_meshes = [view(mesh, get_boundary(mesh)) for mesh in meshes[2:end]]
    fig = viz(
        boundary_meshes[1];
        segmentsize = 1,
        showsegments = true,
        color= :red,
        axis = (aspect = DataAspect(),),
        figure = (; size = (400, 400)),
    )
    
    color_vals = [0.8, 0.5, 0, 0.2, 1]
    for i = 1:num_boundaries
        push!(color_vals, 0.8 - i*(0.2)/num_boundaries)
    end
    colorfier = Colorfy.Colorfier(color_vals)
    colors = Colorfy.colors(colorfier)
    # println(geometries)
    for i in 1:length(geometries)
        if length(geometries[i]) > 0
            viz!(geometries[i]; color=colors[i], showsegments = true, segmentsize = 1)
        end
    end
    for boundary_mesh in boundary_meshes
        println("Showing a boundary")
        viz!(boundary_mesh; color = :red, segmentsize = 1)
    end
    display(fig)
end

function showGeoQuadtreeMesh(forest)
    println("Showing GEO Quadtree mesh!")
    geometries = createMeshByGeoFromQuadtree(forest)
    color_vals = [0.8, 0.5, 0, 0.2, 1]
    for i = 1:num_boundaries
        push!(color_vals, 0.8 - i*(0.2)/num_boundaries)
    end
    colorfier = Colorfy.Colorfier(color_vals)
    colors = Colorfy.colors(colorfier)
    fig = viz(
        geometries[1];
        segmentsize = 1,
        showsegments = true,
        color= colors[1],
        axis = (aspect = DataAspect(),),
        figure = (; size = (1000, 1000)),
    )
    for i in 2:length(geometries)
        if geometries[i] != []
            viz!(geometries[i]; color=colors[i], showsegments = true, segmentsize = 1)
        end
    end
    for boundary in bndry_mesh_points
        for point in boundary
            viz!(Meshes.Point(point...), color=:red, pointsize=6, showpoints=true)
        end
    end
    display(fig)
end

function getBoundaryPoints(parametrization::Function, num_points::Int64)::Vector{Tuple{Float64, Float64}}
    boundary = Vector{Tuple{Float64, Float64}}()
    for t in range(0, stop=1, length=num_points+1)[1:end-1]
        point = parametrization(t)
        push!(boundary, point)
    end
    return boundary
end

function getMeshBounds(boundary_parametrizations::Vector{Function})::Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}
    lb = [Inf, Inf]
    ub = [-Inf, -Inf]
    # using 1000 sample points, this doesn't have to be perfect:
    for parametrization in boundary_parametrizations
        for point in getBoundaryPoints(parametrization, 128)
            for j in 1:2
                if point[j] < lb[j]
                    lb[j] = point[j]
                elseif point[j] > ub[j]
                    ub[j] = point[j]
                end
            end
        end
    end
    return ((lb[1], lb[2]), (ub[1], ub[2]))
end

function check_line_segment(points)
    for point in points[2:end-1]
        if !PolygonAlgorithms.on_segment(point, (points[1], points[end]), atol=1e-12)
            return "Curved"
        end
    end
    return "Straight"
end

function classify_boundary(points, num_points)
    num_segments = Int64(length(points) / num_points)
    classifications = []
    for i in 0:num_segments-1
        tmp_points = points[i*num_points+1:(i+1)*num_points]
        push!(tmp_points, points[mod1((i+1)*num_points+1, length(points))])

        push!(classifications, check_line_segment(tmp_points))
    end
    return classifications
end

function createBoundaryStripMeshes(boundary_associations, internal_points, parametrizations)
    meshes = []
    gmsh.clear()
    boundary_curves = []
    boundary_splines = []
    for i in 1:num_boundaries
        splines = createSplines(parametrizations[i])
        boundary_curve_loop = gmsh.model.occ.addCurveLoop(splines)
        append!(boundary_splines, splines)
        push!(boundary_curves, boundary_curve_loop)
    end
    surfaces = [[boundary_curves[1]]]
    linetags = []
    println(boundary_associations)
    for i in 1:length(boundary_associations)
        association = boundary_associations[i]
        internal_curve_loop, internal_linetags = createCurveLoop(internal_points[i])
        append!(linetags, internal_linetags)
        if 1 in association # main boundary
            push!(surfaces[1], internal_curve_loop)
            for j in association
                if !(boundary_curves[j] in surfaces[1])
                    push!(surfaces[1], boundary_curves[j])
                end
            end
        else
            # then the internal curve must be outside
            boundaries = [boundary_curves[j] for j in association]
            push!(surfaces, [internal_curve_loop, boundaries...])
        end
    end
    for surf in surfaces
        gmsh.model.occ.addPlaneSurface(surf)
        gmsh.model.occ.synchronize()
        for tag in linetags
            gmsh.model.mesh.setTransfiniteCurve(tag, 2)
        end
        gmsh.model.addPhysicalGroup(1, boundary_splines, -1, "Boundary")
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", curvature_size)
        gmsh.model.mesh.removeDuplicateNodes()
        gmsh.model.mesh.generate(2)
        
        msh = Inti.import_mesh(; dim = 2)

        # showMesh(msh)
        
        push!(meshes, msh)
    end

    return meshes
end

function createQuadtree(maxlevel = -1)
    conn = P4estTypes.Connectivity{4}(:unitsquare)
    forest = pxest(conn; min_level=2, init_function = initializeQuadrant, data_type = geo_dom_data)
    refine!(forest; refine = refinementStage, replace = replaceQuadrant, recursive = true, maxlevel = maxlevel)
    coarsen!(forest; coarsen = coarsenStage, replace = replaceQuadrant, recursive = true)
    balance!(forest; replace = replaceQuadrant)
    return forest
end

function get_circle_radius(p1::Tuple{Float64, Float64}, p2::Tuple{Float64, Float64}, p3::Tuple{Float64, Float64})
    e1_dist = distance(p1, p2)
    e2_dist = distance(p2, p3)
    e3_dist = distance(p3, p1)
    area = 1/2 * abs(p1[1] * (p2[2] - p3[2]) + p2[1] * (p3[2] - p1[2]) + p3[1] * (p1[2] - p2[2]))
    return e1_dist*e2_dist*e3_dist / (4*area)
end

function createSplines(parametrizations::Vector{Function}; num_points::Int64 = 64)::Vector{Int32}
    first_point_tag = -1
    splines = Vector{Int32}()
    for (j, func) in enumerate(parametrizations)
        points = getBoundaryPoints(func, num_points)
        point_tags = Vector{Int32}()
        for point in points
            # println(point)
            point_tag = gmsh.model.occ.addPoint(point..., 0.0)
            if first_point_tag == -1
                first_point_tag = point_tag
            end
            push!(point_tags, point_tag)
            # println(point)
        end
        if j == length(parametrizations)
            push!(point_tags, first_point_tag)
            # println(point_tags)
        else
            push!(point_tags, gmsh.model.occ.addPoint(func(1.0)..., 0.0))
        end
            # there are parametrizations that follow this one, should be connected in series
            
        spline = gmsh.model.occ.addSpline(point_tags)
        
        push!(splines, spline)
    end
    return splines
end

function initializeMeshVariables(parametrizations::Vector{Vector{Function}}, forcing_func::Function)
    tmp_bounds = getMeshBounds(parametrizations[1])
    global forcing[] = forcing_func
    tmp_mesh_len = [tmp_bounds[2][1] - tmp_bounds[1][1], tmp_bounds[2][2] - tmp_bounds[1][2]]
    global bounds = Tuple(Tuple(tmp_bounds[i][j] + (-1)^i * tmp_mesh_len[j] * 0.05 for j in 1:2) for i in 1:2)
    global mesh_len = max(bounds[2][1] - bounds[1][1], bounds[2][2] - bounds[1][2])
    # println(bounds)
    # println(mesh_len)
    global bndry_mesh_points = Vector{Vector{Tuple{Float64, Float64}}}()
    global bndry_mesh_orientations = Vector{Int32}()
    global num_boundaries = length(parametrizations)

    # first enforcement on boundary points: must be divisible by four for symmetry
    for i in 1:num_boundaries
        # first check if parametrization is closed
        if distance(parametrizations[i][1](0.0), parametrizations[i][end](1.0)) > 1e-8
            error("Parametrization ", i, " is not closed! start value ", parametrizations[i][1](0.0), " is not equal to end value ", parametrizations[i][end](1.0), "!")
        end
        boundary_points = Vector{Tuple{Float64, Float64}}()
        splines = createSplines(parametrizations[i])
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", curvature_size)
        gmsh.model.mesh.generate(1)
        for spline in splines
            _, node_coords, _ = gmsh.model.mesh.getNodes(1, spline, true, true)
            new_points = [(node_coords[3*k+1], node_coords[3*k+2]) for k in 0:Int64(length(node_coords)/3)-1]
            insert!(new_points, 1, popat!(new_points, length(new_points)-1))   
            new_points = new_points[1:end-1]

            append!(boundary_points, new_points)
        end
        push!(bndry_mesh_points, boundary_points)
        push!(bndry_mesh_splines, splines)
        # println(boundary_points)
        # println("Boundary ", i, " Length: ", length(bndry_mesh_points[i]))
        
        if i == 1
            push!(bndry_mesh_orientations, 1)
        else
            push!(bndry_mesh_orientations, -1)
        end
    end

    if num_boundaries > max_num_boundaries
        error("Too many boundaries! Should be < ", max_num_boundaries)
    end
end

function createQuadtreeMesh(forest)
    gmsh.clear()
    surface_tag = gmsh.model.addDiscreteEntity(2)
    element_type = 3
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

    gmsh.model.occ.synchronize()
    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.mesh.generate(2)
    msh = Inti.import_mesh(; dim = 2)

    # showMesh(msh)

    return msh
end

function getQuadFromPoint(tree::P4estTypes.Tree, point::Tuple{Float64, Float64})::P4estTypes.QuadrantWrapper
    for quad in tree
        coords = QuadToCoords(quad)
        if PolygonAlgorithms.contains(coords, point, atol = 0.0)
            return quad
        end
    end
end

function createQuadtreeMesh(parametrizations::Vector{Vector{Function}}, forcing_func::Function)
    # Creates the meshes describing the boundary region and the quadtree
    # important note: the first function in parametrizations is assumed to be the boundary, the others are assumed to be holes
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 3)
    Inti.clear_entities!()
    println("Creating Quadtree")

    initializeMeshVariables(parametrizations, forcing_func)
    println("Initialized Mesh Variables")

    # println(bndry_mesh_points)
    
    forest = createQuadtree()
    println("Created Quadtree")
    # point = (0.10795, -0.907795)
    # point = (-0.21563770872503343, 0.9951847266721969)
    # quad = getQuadFromPoint(forest[1], point)
    # GEO, DOM = get_GEO_DOM(quad)
    # println(quad)
    # println(GEO)
    # println(DOM)
    # showGeoQuadtreeMesh(forest)
    # error("stop")

    

    quad_mesh = createQuadtreeMesh(forest)
    println("Meshed Quadtree")

    # showMesh(quad_mesh)
    
    point_associations, internal_points = getInternalBoundaryPoints(forest[1])
    println("Created Internal Boundary")

    boundary_strip_meshes = createBoundaryStripMeshes(point_associations, internal_points, parametrizations)
    println("Meshed Boundary Strips")

    meshes = [quad_mesh, boundary_strip_meshes...]

    # showGeoQuadtreeMesh(meshes, forest)

    gmsh.clear()
    gmsh.finalize()

    return meshes
end

# maybe have two functions? one for domain integrals, one for boundary integrals
function getDomainQuadrature(mesh, order)
    domain = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, mesh)
    domain_mesh = Inti.view(mesh, domain)
    domain_quadrature = Inti.Quadrature(domain_mesh; qorder = order)
    return domain_quadrature
end

function getBoundaryQuadrature(mesh, order)
    boundary = get_boundary(mesh)
    boundary_mesh = Inti.view(mesh, boundary)
    boundary_quadrature = Inti.Quadrature(boundary_mesh; qorder = order)
    return boundary_quadrature
end

function calculateIntegrals(meshes, r0; dom_func=nothing, bndry_func=nothing, dom_order=4, bndry_order=6)
    # calculates based on specified functions; if a dom_func is given, calculate domain integral. If a bndry_func is given, calculate boundary integral.
    # first do quadtree mesh
    dom_int = nothing
    bndry_int = nothing
    if !isnothing(dom_func)
        dom_int = 0
        for (i, mesh) in enumerate(meshes)
            domain_quadrature = getDomainQuadrature(mesh, dom_order)
            if i == 1
                dom_int += integrate_quads(dom_func, domain_quadrature, mesh, dom_order, r0)
            else
                dom_int += Inti.integrate(dom_func, domain_quadrature)
            end
        end
    end
    if !isnothing(bndry_func)
        bndry_int = 0
        i = 1
        for (i, mesh) in enumerate(meshes[2:end])
            boundary_quadrature = getBoundaryQuadrature(mesh, bndry_order)
            # check this. based on orientation of boundary quadrature points, assuming now that it comes from bndry_mesh but could be wrong
            bndry_int += Inti.integrate(bndry_func, boundary_quadrature)
        end
    end
    if isnothing(dom_int) && isnothing(bndry_int)
        error("No input functions!")
    elseif isnothing(dom_int) && !isnothing(bndry_int)
        return bndry_int
    elseif !isnothing(dom_int) && isnothing(bndry_int)
        return dom_int
    else
        return dom_int, bndry_int
    end
end

function integrate_quads(func, quadrature, mesh, order, r0)
    sum = 0
    # hard coding order 17 point counts in here temporarily, will change later
    quad_count = 81
    # now separate quadrature points into lists based on distance from r0
    singular_points, near_singular_points, intermediate_points, far_points = [], [], [], []

    # temporary values. dont know how to get them
    singular_threshold = 0.05
    near_singular_threshold = 0.1
    intermediate_threshold = 0.15

    for q in quadrature
        dist = distance(q.coords, r0)
        if dist < singular_threshold
            push!(singular_points, q)
        elseif dist < near_singular_threshold
            push!(near_singular_points, q)
        elseif dist < intermediate_threshold
            push!(intermediate_points, q)
        else
            push!(far_points, q)
        end
    end
    for q in singular_points
        sum += func(q) * q.weight
    end
    for q in near_singular_points
        sum += func(q) * q.weight
    end
    for q in intermediate_points
        sum += func(q) * q.weight
    end
    for q in far_points
        sum += func(q) * q.weight
    end
    return sum
end

function generatePoints(deg, u, x, bounds)
    # x is a vector of deg points between 1 and -1, transform into quad dimensions
    x_normal = [xi / 2 + 1/2 for xi in x]
    xlen, ylen = bounds[1][2] - bounds[1][1], bounds[2][2] - bounds[2][1]
    x_xdim = [bounds[1][1] + xi*xlen for xi in x_normal]
    x_ydim = [bounds[2][1] + xi*ylen for xi in x_normal]
    # println(bounds)
    # println(xlen, ylen)
    # println(x_normal)
    # println(x_xdim)
    # println(x_ydim)
    return [[u((xi, yi)) for xi in x_xdim] for yi in x_ydim]
end

function calculateQuadVolumePotential(quad_mesh, u, target_points)
    # algorithm here:
    #  1 - perform algorithm setup steps (create C matrix)
    #  2 - for each target point: classify quads as singular, near singular, intermediate, and far based on distance from singularity
    #  3 - calculate potentials for each quad based on distance (paper algorithm for singular + near singular, quadrature for intermediate + far)

    # do I want this to have the mesh as input or the quadtree? really the quadtree is the play I think, its a larger object but info can be used 
    # from it for optimization (maybe) and I can construct the quadratures individually
    potentials = [0.0 for point in target_points]
    deg = 5;
    F_map = Dict{P4estTypes.QuadrantWrapper, Matrix{Float64}}()
    L_map = Dict{Tuple{Float64, Float64}, Matrix{Float64}}()
    # first set up legendre polynomials for singular evaluations
    P = ClassicalOrthogonalPolynomials.Legendre()
    
    x = ClassicalOrthogonalPolynomials.grid(P, deg)
    # ux = u.((x, x'))
    # C = ClassicalOrthogonalPolynomials.plan_transform(P, (deg, deg)) * ux

    singular_threshold = 0.05
    near_singular_threshold = 0.1
    intermediate_threshold = 0.15
    # println(forest)
    # separate into four quadratures for singular, near, intermediate, and far points
    for (i, point) in enumerate(target_points)
        singular_quads = Vector{P4estTypes.QuadrantWrapper}()
        near_quads = Vector{P4estTypes.QuadrantWrapper}()
        intermed_quads = Vector{P4estTypes.QuadrantWrapper}()
        far_quads = Vector{P4estTypes.QuadrantWrapper}()
        for quad in tree
            coords = QuadToCoords(quad)
            dist = minimum(distance(coord, point) for coord in coords)
            if dist < singular_threshold
                push!(singular_quads, quad)
            elseif dist < near_singular_threshold
                push!(near_quads, quad)
            elseif dist < intermediate_threshold
                push!(intermed_quads, quad)
            else
                push!(far_quads, quad)
            end
        end
    end
        # println(point)
    for quad in singular_quads
        if haskey(F_map, quad)
            
        else

        end
        coords = QuadToCoords(quad)
        bounds = ((coords[1][1], coords[3][1]), (coords[1][2], coords[3][2]))
        F = generatePoints(deg, u, x, bounds)
        println(typeof(ux))
        C = ClassicalOrthogonalPolynomials.plan_transform(P, (deg, deg)) * ux
        N = Float64.(MultivariateSingularIntegrals.newtoniansquare(big.(point), deg))
        potentials[i] += dot(N, C)
    end
        # do the same for near singular quads
    for quad in intermed_quads
        # get quadrature points for quad by generating grid in similar way to coords
        # then do manual quadrature? two options: implement my own quadrature scheme here, or use the inti quadrature 
    end
    # for now, use the paper technique on singular and near quadratures
    # calculate normally for intermediate and far quadratures (use inti volume potential??? I could also just calculate manually lol)
end

function calculateTriangleVolumePotential(mesh, density, target_points)
    quadrature = getDomainQuadrature(mesh, 17)
    volume_potential = Inti.volume_potential(; 
        op, 
        target = target_points, 
        source = quadrature,
        compression = (method = :none, ),
        correction = (method = :none, ),
    )
    laplacian_points = [density(q.coords) for q in quadrature]
    potentials = volume_potential * laplacian_points
    return potentials
end

function calculateVolumePotential(meshes, u, target_points)
    # calculates the volume potential over the quadratures. Assumes that the first quadrature passed represents the quadtree.
    target_potentials = calculateQuadVolumePotential(meshes[1], u, target_points)
    for mesh in meshes[2:end]
        tmp_target_potentials = calculateTriangleVolumePotential(mesh, u, target_points)
        for (i, potential) in enumerate(tmp_target_potentials)
            target_potentials[i] += potential
        end
    end
    return target_potentials
end

# TODO
# - fix issues when boundaries are close
# - SPEED UP PLEASE
# - change small triangles system
#    - right now we just cap the size of triangle formed for every boundary to meshsize. While the external boundary should have
#      meshsize-sized triangles, we shouldn't require every boundary to be like this, especially inner boundaries
# - fix issues with boundaries sometimes ignoring quads that are closer than other quads that get cut (prob issue with quad classification)

# sizing mesh
# set scalar field to determine local h value (determined by boundary h value)
# https://integralequations.github.io/Inti.jl/stable/tutorials/geo_and_meshes/#Curving-a-given-mesh
# test green's theorem