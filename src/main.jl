# Project Requirements:
# - Look into p4est for managing the quadtree
# - Create a simple test for the solver with a known analytical solution 
#   to compare to on the region [0,1] x [0,1]
# - Create algorithm for seeing if quad should be refined, based on 
#   function values with chebyshev points
# - Create heatmap display for integral error at different levels of 
#   quadtree depth, generally create multiple methods of display
include("mesh_tools.jl")

using P4estTypes;
using Polynomials
using FastChebInterp
using Inti, Meshes, CairoMakie, StaticArrays, Gmsh
using PolygonAlgorithms
using GeometryBasics
using Colorfy
using FMM2D
using IterativeSolvers, LinearAlgebra
# using BenchmarkTools

@enum GEOS null=0 cut=1 regular=2 contains=3 outside=4

p = 12 # deg of interpolating polys for function refinement
len = 2^30# size of default quadtree
contain_tol = 1e-15

mesh_len = -1
num_boundaries = -1
bndry_mesh_points = []
bndry_mesh_orientations = []
global const forcing = Ref{Function}((x) -> 0)


const max_num_boundaries = 16

struct geo_dom_data
    GEO::GEOS
    DOM::NTuple{max_num_boundaries, Int32}
    PARENT_GEO::GEOS
    PARENT_DOM::NTuple{max_num_boundaries, Int32}
end

r0 = (-0.5, 0)
# forcing(x) = -2*pi^2 * sin(pi*x[1]) * sin(pi*x[2]) + exp(-600*((x[1]- r0[1])^2 + (x[2] - r0[2])^2))
# forcing(x) = exp(-100 * (x[1]^2 + x[2]^2))
# forcing(x) = exp(x[1]+x[2])

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
            geo_dom_val_tmp = geo_dom_data(null, null_list, parent_GEO, parent_DOM_NTuple)
            P4estTypes.unsafe_storeuserdata!(quadrant, geo_dom_val_tmp)
            GEO, DOM = classifyQuad(quadrant, bndry_mesh_points, bndry_mesh_orientations, 0.1)
            DOM_list = [x for x in DOM] 
            DOM_list = vcat(DOM_list, [-1 for i in 1:max_num_boundaries-length(DOM_list)])
            DOM_NTuple = NTuple{max_num_boundaries, Int32}(DOM_list)
            geo_dom_val = geo_dom_data(GEO, DOM_NTuple, parent_GEO, parent_DOM_NTuple)
            P4estTypes.unsafe_storeuserdata!(quadrant, geo_dom_val)
        end
    end
end

function getInternalBoundaryPoints(tree)
    # new algorithm: we want every point that is on the "outside" of the quadtree; essentially every cut

    # first find lowest quad level
    max_level = 0
    for quad in tree
        GEO, DOM = get_GEO_DOM(quad)
        lvl = P4estTypes.level(quad)
        if (GEO == regular || GEO == cut) && lvl > max_level 
            max_level = lvl
        end
    end
    edge_list = []
    for boundary in 1:num_boundaries
        point_connection_dict = Dict{Tuple{Float64, Float64}, Set{Tuple{Float64, Float64}}}()
        edge_counts = Dict{Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}, Int32}()
        for quadrant in tree
            GEO, DOM = get_GEO_DOM(quadrant)
            if (GEO == regular && boundary in DOM) || (GEO == cut && boundary in DOM) && boundary == minimum(DOM)
                level_diff = max_level - P4estTypes.level(quadrant)
                coords = QuadToCoords(quadrant)
                for i in 1:4
                    # add edge for each point to connection dict
                    point, right_point = coords[i], coords[mod1(i+1, 4)]                    
                    # if !haskey(point_connection_dict, point)
                    #     point_connection_dict[point] = Set{Tuple{Float64, Float64}}()
                    # end
                    # push!(point_connection_dict[point], middle_point, right_point)
                    # push!(point_connection_dict[middle_point], point, right_point)
                    # push!(point_connection_dict[right_point], point, middle_point)
                    intermediary_points = [((point[1]+right_point[1])/2*i, (point[2]+right_point[2])/2*i) for i in 1:2^level_diff-1]
                    all_edge_points = [point; intermediary_points; right_point]
                    for i in 1:length(all_edge_points)-1
                        edge = (all_edge_points[i], all_edge_points[i+1])
                        reverse_edge = (all_edge_points[i+1], all_edge_points[i])
                        found_edge = false
                        for other_edge in keys(edge_counts)
                            if distance(edge[1], other_edge[1]) < tol && distance(edge[2], other_edge[2]) < tol
                                edge_counts[other_edge] += 1
                                found_edge = true
                                break
                            elseif distance(reverse_edge[1], other_edge[1]) < tol && distance(reverse_edge[2], other_edge[2]) < tol
                                edge_counts[other_edge] += 1
                                found_edge = true
                                break
                            end
                        end
                        if !found_edge
                            edge_counts[edge] = 1
                        end
                    end
                end
            end
        end
        good_edges = []
        for edge in keys(edge_counts)
            count = edge_counts[edge]
            if count == 1 && PolygonAlgorithms.contains(bndry_mesh_points[boundary], edge[1], atol=1e-15)
                push!(good_edges, edge)
            end
        end
        start_edge = good_edges[1]
        edge_connections = [start_edge[1], start_edge[2]]
        k = 2
        while k < length(good_edges)
            for edge in good_edges
                if distance(edge[1], edge_connections[end]) < tol && distance(edge[2], edge_connections[end-1]) >= tol
                    push!(edge_connections, edge[2])
                    break
                elseif distance(edge[2], edge_connections[end]) < tol && distance(edge[1], edge_connections[end-1]) >= tol
                    push!(edge_connections, edge[1])
                    break
                end
            end
            k += 1
        end
        push!(edge_list, edge_connections)
    end
    return edge_list
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
    # viz!(Γ_msh; color = :red, segmentsize = 1)
    color_vals = [0.8, 0.5, 0, 0.2, 1]
    for i = 1:num_boundaries
        push!(color_vals, 0.8 - i*(0.2)/num_boundaries)
    end
    colorfier = Colorfy.Colorfier(color_vals)
    colors = Colorfy.colors(colorfier)
    for i in 1:length(geometries)
        viz!(geometries[i]; color=colors[i], showsegments = true, segmentsize = 1)
    end
    viz!(Γ_msh; color = :red, segmentsize = 1)
    display(fig)
end

function createMeshFromQuadtree(forest)
    box_list = []
    for quad in forest[1]
        coords = QuadToCoords(quad)
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

function getBoundaryPoints(parametrization::Function, num_points::Int64)::Vector{Tuple{Float64, Float64}}
    boundary = Vector{Tuple{Float64, Float64}}()
    for t in range(0, stop=1, length=num_points)
        point = parametrization(t)
        push!(boundary, point)
    end
    return boundary
end

function getMeshBounds(boundary_parametrization::Function)::Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}
    lb = [Inf, Inf]
    ub = [-Inf, -Inf]
    # using 1000 sample points, this doesn't have to be perfect:
    for i in range(0, stop = 1, length = 1000)
        val = boundary_parametrization(i)
        for j in 1:2
            if val[j] < lb[j]
                lb[j] = val[j]
            elseif val[j] > ub[j]
                ub[j] = val[j]
            end
        end
    end
    return ((lb[1], lb[2]), (ub[1], ub[2]))
end

function createMeshWithInternalBoundary(boundary_points, internal_points, meshsize, max_triangle_size, parametrizations)
    surface_tags = []
    all_boundary_linetags = []
    all_internal_linetags = []
    for i in 1:num_boundaries
        # boundary_linetag = Inti.gmsh_curve(0, 1; npts = 100, meshsize = max_triangle_size) do s
        #     return Inti.Point2D(parametrizations[i](s)[1], parametrizations[i](s)[2])
        # end
        points = []
        for t in range(0, stop = 1, length=Int32(ceil(100 / meshsize)))[1:end-1]
            point = (parametrizations[i](t)[1], parametrizations[i](t)[2])
            push!(points, gmsh.model.geo.addPoint(point..., 0.0, meshsize))
        end
        push!(points, points[1])
        polyloop = gmsh.model.geo.addPolyline(points)
        boundary_curve_loop = gmsh.model.geo.addCurveLoop([polyloop])
        internal_curve_loop, internal_linetags = createInternalCurveLoop(internal_points[i], max_triangle_size)
        if bndry_mesh_orientations[i] == 1
            surf = gmsh.model.geo.addPlaneSurface([boundary_curve_loop, internal_curve_loop])
        else
            surf = gmsh.model.geo.addPlaneSurface([internal_curve_loop, boundary_curve_loop])
        end
        push!(surface_tags, surf)
        push!(all_boundary_linetags, polyloop)
        push!(all_internal_linetags, internal_linetags...)
    end
    return all_boundary_linetags, all_internal_linetags, surface_tags
end

function createQuadtree(maxlevel = -1)
    conn = P4estTypes.Connectivity{4}(:unitsquare)
    forest = pxest(conn; min_level=0, init_function = initializeQuadrant, data_type = geo_dom_data)
    refine!(forest; refine = refinementStage, replace = replaceQuadrant, recursive = true, maxlevel = maxlevel)
    coarsen!(forest; coarsen = coarsenStage, init = initializeQuadrant, recursive = true)
    balance!(forest; init = initializeQuadrant)
    return forest
end

function initializeMeshVariables(parametrizations::Vector{Function}, forcing_func::Function)
    tmp_bounds = getMeshBounds(parametrizations[1])
    global forcing[] = forcing_func
    tmp_mesh_len = [tmp_bounds[2][1] - tmp_bounds[1][1], tmp_bounds[2][2] - tmp_bounds[1][2]]
    global bounds = Tuple(Tuple(tmp_bounds[i][j] + (-1)^i * tmp_mesh_len[j] * 0.05 for j in 1:2) for i in 1:2)
    global mesh_len = [bounds[2][1] - bounds[1][1], bounds[2][2] - bounds[1][2]]
    # println(bounds)
    # println(mesh_len)
    global bndry_mesh_points = Vector{Vector{Tuple{Float64, Float64}}}()
    global bndry_mesh_orientations = Vector{Int32}()
    global num_boundaries = length(parametrizations)

    for i in 1:num_boundaries
        boundary_points = getBoundaryPoints(parametrizations[i], 100)
        # println(boundary_points)
        push!(bndry_mesh_points, boundary_points)
        if i == 1
            push!(bndry_mesh_orientations, 1)
        else
            push!(bndry_mesh_orientations, -1)
        end
    end    
    # println(bndry_mesh_points)
    # println(bndry_mesh_orientations)
    # println(num_boundaries)

    if num_boundaries > max_num_boundaries
        error("Too many boundaries! Should be < ", max_num_boundaries)
    end
end

function createCombinedMesh(forest, surface_tag)
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
end

function createQuadtreeMesh(parametrizations::Vector{Function}, forcing_func::Function, meshsize::Float64, max_triangle_size::Float64)
    # Creates the meshes describing the boundary region and the quadtree
    # important note: the first function in parametrizations is assumed to be the boundary, the others are assumed to be holes
    
    gmsh.option.setNumber("General.Verbosity", 3)

    initializeMeshVariables(parametrizations, forcing_func)
    println("wassup")
    forest = createQuadtree(-1)
    println(forest[1])
    # println("hey")
    internalBoundaryPoints = getInternalBoundaryPoints(forest[1])
    
    gmsh.model.add("BoundaryRegions")

    boundary_linetags, internal_linetags, surface_tags = createMeshWithInternalBoundary(bndry_mesh_points, internalBoundaryPoints, meshsize, max_triangle_size, parametrizations)
    println("hi")
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.mesh.generate(2)
    println("no")
    for i in 1:num_boundaries
        boundary_str = string("Boundary ", i)
        gmsh.model.addPhysicalGroup(2, [surface_tags[i]], -1, boundary_str)
        gmsh.model.addPhysicalGroup(1, [boundary_linetags[i]], -1, boundary_str)
    end

    surface_tag = gmsh.model.addDiscreteEntity(2, -1, internal_linetags)

    createCombinedMesh(forest, surface_tag)
    
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(2, [surface_tag], -1, "Quadtree")
    
    gmsh.model.mesh.removeDuplicateNodes()
    
    gmsh.write("combined.msh")
    
    Inti.clear_entities!()
    combined_msh = Inti.import_mesh(; dim = 2)

    gmsh.clear()
    
    gmsh.finalize()

    return combined_msh
end

function calculateVolumePotential(dom_mesh, f)
    dom_quad = Inti.Quadrature(dom_mesh; qorder = 4)

    greens_fn = (r, x) -> 1/(2pi) * log(distance(x, r))

    volume_potential = (x) -> Inti.integrate((q) -> f(q.coords) * greens_fn(x, q.coords), dom_quad)
    return volume_potential
end

# driver code
# gmsh.initialize()
# gmsh.option.setNumber("General.Verbosity", 3) # turn to 4/5 for info or debug, 3 is all that is necessary I think though
# gmsh.open("testing.msh")
# msh = Inti.import_mesh(; dim = 2)
# calculateMeshvals(msh)
# mesh = createMesh(0.25, 0.25)
# showMesh(mesh)

# showMesh(mesh)
# gmsh.initialize()
# parametrizations::Vector{Function} = [(x) -> [cos(2*pi*x), sin(2*pi*x)]]
# quadtree_mesh = createQuadtreeMesh(mesh, (x) -> 1, 0.25, 0.25, parametrizations)
# showMesh(quadtree_mesh)
# dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, quadtree_mesh)


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