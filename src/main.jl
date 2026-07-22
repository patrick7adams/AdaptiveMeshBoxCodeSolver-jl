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
curvature_size = 200

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
            # if (GEO == regular && length(boundary) > 0) || (GEO == cut && length(boundary) > 0)
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
        # showEdgeList(good_edges)
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
            edge_connections = [getFirstPoints(good_edges)...]
            if Tuple(edge_connections) in good_edges
                pop!(good_edges, Tuple(edge_connections))
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
            if edge_list_associations[i] in association_set && length(association_set) > 1
                found_before = true
            end
        end
        if found_before
            continue
        end
        edge_list_association = [i]
        j = 1
        repeat = true
        count = 0
        while repeat && count < length(edge_list)
            repeat = false
            for (j, other_edge_set) in enumerate(edge_list)
                break_k = false
                for k in 1:5
                    # check first few elements of the next one for distance then cut off there
                    if distance(edge_set[end], other_edge_set[k]) <= 1e-10
                        if !(edge_list_associations[j] in edge_list_association)
                            push!(edge_list_association, edge_list_associations[j])
                        end
                        append!(edge_set, other_edge_set[k+1:end])
                        repeat = true
                        if distance(edge_set[1], edge_set[end]) <= 1e-10
                            repeat = false
                        end
                        break_k = true
                        break
                    end
                end
                if break_k
                    break
                end
            end
            count += 1
        end
        push!(final_edge_list_associations, edge_list_association)
        push!(final_edge_list, edge_set[1:end-1])
    end
    # for edge_set in final_edge_list
    #     showEdgeList(edge_set)
    # end
    # println(final_edge_list)
    # for i in 1:length(final_edge_list[1])
    #     if distance(final_edge_list[1][i], final_edge_list[1][mod1(i+1, length(final_edge_list[1]))]) < 1e-8
    #         println(final_edge_list[1][i])
    #         println(final_edge_list[1][mod1(i+1, length(final_edge_list[1]))])
    #         println(i)
    #     end
    # end
    # println(final_edge_list_associations)
    # error("hi")
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
                end
                if point[j] > ub[j]
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

function refinementTriangleStage(triangle, quad_points, func, Qw, R, degs, Iorder)
    # takes in a triangle, returns true if triangle should refine, false if not
    c, r = Inti.translation_and_scaling(triangle)
    d = func.(Inti.coords.(quad_points))
    # w = Inti.weight.(quad_points)
    mat = zeros(length(quad_points), length(quad_points))
    for (i, point) in enumerate(quad_points)
        pnt = 1/r * (point.coords .- c)
        for (j, deg) in enumerate(degs)
            mat[i, j] = (1/factorial(deg[1]))*pnt[1]^deg[1] * (1/factorial(deg[2]))*pnt[2]^deg[2]
        end
    end
    coeffs = mat \ d
    # coeffs = R * (Qw' * (sqrt.(w) .* d))
    # coeffs = Qw' * d
    
    good_coeffs = []
    for (deg, coeff) in zip(degs, coeffs)
        if deg[1]+deg[2] == Iorder
            push!(good_coeffs, coeff)
        end
    end
    rel_sum = norm(good_coeffs) / max(norm(coeffs), 1e-12)
    abs_sum = norm(good_coeffs)

    return abs_sum > 1e-1 && rel_sum > 1e-1
end

function generateMesh(surf, boundary_splines, linetags, func, Iorder, maxiter = 5)
    qorder = Inti.TRIANGLE_VR_IORDER_TO_QORDER[Iorder]
    num_points_per_quad = Inti.TRIANGLE_VR_ORDER_TO_NPTS[qorder]
    qrule = Inti.TRIANGLE_VR_QRULES[num_points_per_quad]

    size_field = -1

    degs = sort([(i, j) for i in 0:Iorder, j in 0:Iorder if i+j <= Iorder], by = sum)
    len_vander = Int64((Iorder+1)*(Iorder+2) / 2)

    # create reference vandermonde matrix for the given qrule
    mat = zeros(num_points_per_quad, len_vander)
    for (i, node) in enumerate(qrule)
        pnt = node[1]
        for (j, deg) in enumerate(degs)
            mat[i, j] = (1/factorial(deg[1]))*pnt[1]^deg[1] * (1/factorial(deg[2]))*pnt[2]^deg[2]
        end
    end

    A = sqrt.([node[2] for node in qrule]) .* mat
    F = qr(A)
    Qw = Matrix(F.Q)[:, 1:size(mat, 2)]
    R = Matrix(F.R)

    gmsh.model.occ.addPlaneSurface(surf)
    gmsh.model.occ.synchronize()
    # for tag in linetags
    #     gmsh.model.mesh.setTransfiniteCurve(tag, 2)
    # end
    gmsh.model.addPhysicalGroup(1, boundary_splines, -1, "Boundary")
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", curvature_size)
    gmsh.model.mesh.removeDuplicateNodes()

    for iter in 1:maxiter
        gmsh.model.mesh.clear()
        gmsh.model.mesh.generate(2)
        mesh = Inti.import_mesh(; dim=2)

        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, mesh)
        dom_mesh = view(mesh, dom)
        quadrature = Inti.Quadrature(dom_mesh; qorder = qorder)
        
        c1, c2 = 0, 0
        vertex_dict = Dict()
        data = []

        triangles = Inti.elements(mesh, Inti.LagrangeElement{Inti.ReferenceSimplex{2}, 3, StaticArraysCore.SVector{2, Float64}})
        for (i, elem) in enumerate(triangles)
            v1, v2, v3 = ((round(v[1], digits=8), round(v[2], digits=8)) for v in Inti.vertices(elem))
            quad_points = quadrature[(i-1)*num_points_per_quad+1:i*num_points_per_quad]
            if iter > 1
                center = Inti.center(elem)
                target, _ = gmsh.view.probe(size_field, center[1], center[2], 0.0)
                size = target[1]
            else
                size = maximum((norm(v1.-v2), norm(v2.-v3), norm(v1.-v3)))
            end

            if refinementTriangleStage(elem, quad_points, func, Qw, R, degs, Iorder)
                size /= 2.0
                c2 += 1
            else
                c1 += 1
            end

            for v in (v1, v2, v3)
                if haskey(vertex_dict, v)
                    push!(vertex_dict[v], size)
                else
                    vertex_dict[v] = [size]
                end
            end
        end

        if c2 == 0
            break
        end

        filtered_size_dict = Dict()
        for (k, v) in vertex_dict
            newsum = minimum(v)
            filtered_size_dict[k] = newsum
        end

        for (i, elem) in enumerate(triangles)
            v1, v2, v3 = ((round(v[1], digits=8), round(v[2], digits=8)) for v in Inti.vertices(elem))
            size1 = filtered_size_dict[v1]
            size2 = filtered_size_dict[v2]
            size3 = filtered_size_dict[v3]
            push!(data, v1[1], v2[1], v3[1], v1[2], v2[2], v3[2], 0.0, 0.0, 0.0, size1, size2, size3)
        end

        size_field = gmsh.view.add("mesh size field")
        gmsh.view.addListData(size_field, "ST", length(triangles), data)

        bg_field = gmsh.model.mesh.field.add("PostView")
        gmsh.model.mesh.field.setNumber(bg_field, "ViewTag", size_field)
        gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)
    end

    mesh = Inti.import_mesh(; dim=2)
    return mesh
end

function createBoundaryStripMeshes(boundary_associations, internal_points, parametrizations, func, Iorder)
    function getCurve(i)
        splines = createSplines(parametrizations[i])
        return gmsh.model.occ.addCurveLoop(splines), splines
    end

    iter = 1
    meshes = []
    gmsh.clear()
    boundary_splines = Vector{Int32}()
    surfaces = Vector{Int32}()
    linetags = Vector{Int32}()

    external_associations = [association for association in boundary_associations if 1 in association]
    internal_associations = [association for association in boundary_associations if !(1 in association)]
    
    boundary_curve, splines = getCurve(1)
    push!(surfaces, boundary_curve)
    push!(boundary_splines, splines...)

    for (i, association) in enumerate(external_associations)
        internal_curve_loop, internal_linetags = createCurveLoop(internal_points[i])
        append!(linetags, internal_linetags)
        push!(surfaces, internal_curve_loop)
        for j in association
            j != 1 || continue
            curve, splines = getCurve(j)
            push!(surfaces, curve)
            push!(boundary_splines, splines...)
        end
    end
    msh = generateMesh(surfaces, boundary_splines, linetags, func, Iorder)
    push!(meshes, msh)
    gmsh.clear()
    for (i, association) in enumerate(internal_associations)
        boundary_splines = Vector{Int32}()
        surfaces = Vector{Int32}()
        linetags = Vector{Int32}()
        
        internal_curve_loop, internal_linetags = createCurveLoop(internal_points[i+length(external_associations)])
        append!(linetags, internal_linetags)

        for j in association
            curve, splines = getCurve(j)
            push!(surfaces, curve)
            push!(boundary_splines, splines...)
        end
        msh = generateMesh([internal_curve_loop, surfaces...], boundary_splines, linetags, func, Iorder, maxiter = 15)
        push!(meshes, msh)
        gmsh.clear()
    end
    # print(meshes)
    # for mesh in meshes
    #     showMesh(mesh)
    # end

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

function createSplines(parametrizations::Vector{Function}; num_points::Int64 = 16)::Vector{Int32}
    first_point_tag = -1
    splines = Vector{Int32}()
    for (j, func) in enumerate(parametrizations)
        points = getBoundaryPoints(func, num_points)
        # println("-----")
        # println(points)
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
        # println(point_tags)
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

function createQuadtreeMesh(parametrizations::Vector{Vector{Function}}, forcing_func::Function, verbose=false, Iorder=4)
    # Creates the meshes describing the boundary region and the quadtree
    # important note: the first function in parametrizations is assumed to be the boundary, the others are assumed to be holes
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 3)
    Inti.clear_entities!()
    if verbose
        println("Creating Quadtree")
    end

    initializeMeshVariables(parametrizations, forcing_func)
    if verbose
        println("Initialized Mesh Variables")
    end

    # println(bndry_mesh_points)
    
    forest = createQuadtree()
    if verbose
        println("Created Quadtree")
    end

    quad_mesh = createQuadtreeMesh(forest)
    if verbose
        println("Meshed Quadtree")
    end
    
    point_associations, internal_points = getInternalBoundaryPoints(forest[1])
    if verbose
        println("Created Internal Boundary")
    end

    boundary_strip_meshes = createBoundaryStripMeshes(point_associations, internal_points, parametrizations, forcing_func, Iorder)
    if verbose
        println("Meshed Boundary Strips")
    end

    meshes = [quad_mesh, boundary_strip_meshes...]

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
    # domain = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, mesh)
    # boundary = Inti.boundary(domain)
    boundary_mesh = Inti.view(mesh, boundary)
    boundary_quadrature = Inti.Quadrature(boundary_mesh; qorder = order)
    return boundary_quadrature
end

function generatePoints(u, x, bounds)
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

    return stack([[u((xi, yi)) for xi in x_xdim] for yi in x_ydim], dims=1)
end

function separateMesh(meshes, quadrature, point)
    function roundPoint(pnt, digits=8)
        return (round(pnt[1], digits=digits), round(pnt[2], digits=digits))
    end
    # separates the mesh into:
    mesh = meshes[1]
    # 1. a quadrature of the nonsingular points
    # 2. a list of elements of the mesh for the singular / near-singular quads

    # first get the singular elements
    good_elements = []
    for elem in Inti.elements(mesh)
        coords = [(x[1], x[2]) for x in Inti.vertices(elem)]
        source_h = coords[3][1] - coords[1][1]
        center = ((coords[3][1] + coords[1][1])/2, (coords[3][2] + coords[1][2])/2)
        center_dist = distance(center, point)
        relative_dist = center_dist / source_h
        if relative_dist < 1.5
            push!(good_elements, elem)
        end
    end

    # NOW from the separated mesh we need to get the quadrature of the rest of the mesh
    n = length([x for x in Inti.elements(mesh)])
    points_per_quad = Int64(length(quadrature) / n)
    culled_quadrature = Vector{Inti.QuadratureNode}()
    for i in 0:n-1
        is_singular = false
        # get some quadrature point

        test_point = quadrature[Int64(floor((i+0.5)*points_per_quad))+1].coords
        for elem in good_elements
            coords = [(x[1], x[2]) for x in Inti.vertices(elem)]
            if (coords[1][1] <= test_point[1] && test_point[1] <= coords[3][1]) && (coords[1][2] <= test_point[2] && test_point[2] <= coords[3][2])
                is_singular = true
                break
            end
        end
        if !is_singular
            append!(culled_quadrature, quadrature[i*points_per_quad+1:(i+1)*points_per_quad])
        end
    end
    # showSeparatedMesh(meshes, good_elements, culled_quadrature, point)
    return good_elements, culled_quadrature
end

function getCorrectionMap(quad_mesh, order, u)
    Correction_map = Dict{Inti.LagrangeElement{Inti.ReferenceHyperCube{2}, 4, StaticArraysCore.SVector{2, Float64}}, Float64}()

    quadrature = Inti.Quadrature(quad_mesh; qorder=order)
    num_points_per_quad = Int64(ceil((order + 1) / 2))^2
    num_quads = Int64(length(quadrature) / num_points_per_quad)
    mid_point_coord = Int64(ceil(num_points_per_quad/2))
    for quad in Inti.elements(quad_mesh)
        coords = [(x[1], x[2]) for x in Inti.vertices(quad)]
        h = coords[3][1] - coords[1][1]
        start_idx = -1
        for i in 1:num_quads
            point = quadrature[(i-1)*num_points_per_quad + mid_point_coord].coords
            if (coords[1][1] <= point[1] && point[1] <= coords[3][1]) && (coords[1][2] <= point[2] && point[2] <= coords[3][2])
                start_idx = (i-1)*num_points_per_quad
                break
            end
        end
        integral = 0
        for quad_point in quadrature[start_idx+1:start_idx+num_points_per_quad]
            integral += u(quad_point.coords) * quad_point.weight 
        end
        correction_term = integral
        Correction_map[quad] = correction_term
    end
    return Correction_map
end

function calculateQuadVolumePotential(meshes, quadrature, u, target_points)
    quad_mesh = meshes[1]
    greens_fn = (x, y) -> 1/(2pi) * log(distance(x, y))
    potentials = [0.0 for point in target_points]
    deg = 5;
    F_map = Dict{Inti.LagrangeElement{Inti.ReferenceHyperCube{2}, 4, StaticArraysCore.SVector{2, Float64}}, Matrix{Float64}}()
    L_map = Dict{Vector{Float64}, Matrix{Float64}}()
    Correction_map = getCorrectionMap(quad_mesh, 5, u)

    # first set up legendre polynomials for singular evaluations
    P = ClassicalOrthogonalPolynomials.Legendre()
    
    x = ClassicalOrthogonalPolynomials.grid(P, deg)
    F_cache_hits = 0
    L_cache_hits = 0
    total_calcs = 0
    F_calc_time = 0
    L_calc_time = 0
    separate_mesh_time = 0
    @show length(target_points)
    true_start_time = time()
    for (i, point) in enumerate(target_points)
        # println(Float64(i) / Float64(length(target_points)))
        start_time = time()
        singular_quads, far_quadrature = separateMesh(meshes, quadrature, point)
        separate_mesh_time += time() - start_time
        for elem in singular_quads
            total_calcs += 1
            # first create F
            coords = [(x[1], x[2]) for x in Inti.vertices(elem)]
            
            bounds = ((coords[1][1], coords[3][1]), (coords[1][2], coords[3][2]))
            h = coords[3][1] - coords[1][1]
            if haskey(F_map, elem)
                F = F_map[elem]
                F_cache_hits += 1
            else
                start_time = time()
                ux = generatePoints(u, x, bounds)'
                tmp_F = ClassicalOrthogonalPolynomials.plan_transform(P, (deg, deg))
                F = (tmp_F * ux)
                F_map[elem] = F
                F_calc_time += time() - start_time
            end
            # now create L by first mapping point to the [[-1, 1], [-1, 1]] grid
            x_dist, y_dist = point[1] - coords[1][1], point[2] - coords[1][2]
            
            mapped_target_point = [round(-1.0+2/h*x_dist, digits=12), round(-1.0+2/h*y_dist, digits=12)]
            mapped_target_point += [1e-16, 1e-16]
            if haskey(L_map, mapped_target_point)
                L = L_map[mapped_target_point]
                L_cache_hits += 1
            else
                start_time = time()
                L = Float64.(MultivariateSingularIntegrals.newtoniansquare(big.(mapped_target_point), deg))
                L_map[mapped_target_point] = L
                L_calc_time += time() - start_time
            end
            I = dot(L, F)
            term = I*h^2 / (8*pi) - log(2/h)/(2*pi)*Correction_map[elem]
            potentials[i] += term

            # if elem in singular_quads
            #     println("--------")
            #     println("Coords: ", coords)
            #     println("Quad Contribution: ", term)
            #     println("Density integral over quad: ", I)
            #     println("Relative target point position: ", mapped_target_point)
            # end
        end

     
        domain_func = (q) -> greens_fn(point, q.coords) * u(q.coords)
        potentials[i] += sum(domain_func(q)*q.weight for q in far_quadrature)

        # for elem in Inti.elements(quad_mesh)
        #     println(Inti.vertices(elem))
        # end
        # println(sum(domain_func(q)*q.weight for q in far_quadrature))
    end
    full_run_time = time() - true_start_time
    @show full_run_time
    @show F_calc_time
    @show L_calc_time
    @show total_calcs
    @show F_cache_hits
    @show L_cache_hits
    @show separate_mesh_time
    return potentials
end

function calculateTriangleVolumePotential(quadrature, density, target_points, multiplicative_terms)
    op = Inti.Laplace(; dim=2)
    inside_target_points = Vector{Tuple{Float64, Float64}}()
    boundary_target_points = Vector{Tuple{Float64, Float64}}()
    outside_target_points = Vector{Tuple{Float64, Float64}}()
    inside_indices = Vector{Int64}()
    boundary_indices = Vector{Int64}()
    outside_indices = Vector{Int64}()
    # println(multiplicative_terms)
    for (i, point) in enumerate(target_points)
        if multiplicative_terms[i] == :on
            push!(boundary_target_points, point)
            push!(boundary_indices, i)
        elseif multiplicative_terms[i] == :inside
            push!(inside_target_points, point)
            push!(inside_indices, i)
        else
            push!(outside_target_points, point)
            push!(outside_indices, i)
        end
    end
    max_dist = 0.2
    potentials = [0.0 for x in target_points]
    laplacian_points = [density(q.coords) for q in quadrature]

    compression_method = (method = :fmm, tol=1e-12)

    if length(inside_target_points) > 0
        inside_volume_potential = Inti.volume_potential(; 
            op, 
            target = inside_target_points, 
            source = quadrature,
            compression = compression_method,
            correction = (method = :dim, maxdist=max_dist, target_location=:inside), 
        )
        inside_potentials = inside_volume_potential * laplacian_points
        for i in 1:length(inside_potentials)
            potentials[inside_indices[i]] = -inside_potentials[i] 
        end
    end
    if length(boundary_target_points) > 0
        boundary_volume_potential = Inti.volume_potential(; 
            op, 
            target = boundary_target_points, 
            source = quadrature,
            compression = compression_method,
            correction = (method = :dim, maxdist=max_dist, target_location=:on), 
        )
        boundary_potentials = boundary_volume_potential * laplacian_points
        for i in 1:length(boundary_potentials)
            potentials[boundary_indices[i]] = -boundary_potentials[i] 
        end
    end
    if length(outside_target_points) > 0
        outside_volume_potential = Inti.volume_potential(; 
            op, 
            target = outside_target_points, 
            source = quadrature,
            compression = compression_method,
            correction = (method = :dim, maxdist=max_dist, target_location=:outside), 
        )
        outside_potentials = outside_volume_potential * laplacian_points
        for i in 1:length(outside_potentials)
            potentials[outside_indices[i]] = -outside_potentials[i] 
        end
    end
    # println(potentials)
    return potentials
end

function calculateVolumePotential(quadratures, meshes, u, target_points, multiplicative_terms, verbose=false)
    # calculates the volume potential over the quadratures. Assumes that the first quadrature passed represents the quadtree.
    target_potentials = calculateQuadVolumePotential(meshes, quadratures[1], u, target_points)
    if verbose
        println("Finished calculating quad volume potential")
    end
    for (i, quadrature) in enumerate(quadratures[2:end])
        tmp_target_potentials = calculateTriangleVolumePotential(quadrature, u, target_points, multiplicative_terms[i+1])
        for (j, potential) in enumerate(tmp_target_potentials)
            target_potentials[j] += potential
        end
        if verbose
            println("Finished calculating triangle volume potential for boundary ", i)
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