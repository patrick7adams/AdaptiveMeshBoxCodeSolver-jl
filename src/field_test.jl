using AdaptiveMeshSolver
using Inti
using LinearAlgebra
using Printf
using Gmsh
using CairoMakie
using Meshes
using StaticArrays
using Statistics

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
    # if length(points) > 0
    #     pointset = Meshes.PointSet([Meshes.Point((point[1], point[2])) for point in points])
    #     viz!(pointset; color = :red, pointsize = 3)
    # end
    display(fig)
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

function getBoundaryPoints(parametrization::Function, num_points::Int64)::Vector{Tuple{Float64, Float64}}
    boundary = Vector{Tuple{Float64, Float64}}()
    for t in range(0, stop=1, length=num_points+1)[1:end-1]
        point = parametrization(t)
        push!(boundary, point)
    end
    return boundary
end

function refinementTriangleStage(triangle, quad_points, func, degree)
    # takes in a triangle, returns true if triangle should refine, false if not
    # c, r = Inti.translation_and_scaling(triangle)
    c, r = Inti.translation_and_scaling(triangle)

    d = func.(Inti.coords.(quad_points))
    mat = zeros(num_points_per_quad, len_vander)
    for (i, point) in enumerate(quad_points)
        pnt = 1/r * (point.coords .- c)
        for (j, deg) in enumerate(degs)
            mat[i, j] = (1/factorial(deg[1]))*pnt[1]^deg[1] * (1/factorial(deg[2]))*pnt[2]^deg[2]
        end
    end
    coeffs = mat \ d
    
    good_coeffs = []
    for (deg, coeff) in zip(degs, coeffs)
        if deg[1]+deg[2] == degree
            push!(good_coeffs, coeff)
        end
    end
    rel_sum = norm(good_coeffs) / max(norm(coeffs), 1e-12)
    abs_sum = norm(good_coeffs)
    if triangle == Inti.LagrangeElement{Inti.ReferenceSimplex{2}, 3, SVector{2, Float64}}(SVector{2, Float64}[[-0.08570188488096438, 0.016586422460625295], [-0.08973315785630873, 0.012554927151873826], [-0.08442082395404617, 0.008774070863537256]])
        println(coeffs)
        println(mat)
        for (i, point) in enumerate(quad_points)
            println(1/r * (point.coords .- c))
            println(point.coords .- c)
        end
        @show r
        @show c
    end

    return abs_sum > 1e-1 && rel_sum > 1e-1, cond(mat)
end

# get coeffs of one triangle, refine

# function getCoeffs(d, w)
#     # return Qw' * (sqrt.(w) .* d)
#     return R \ (Qw' * (sqrt.(w) .* d))
# end

σ=0.05
x_0 = 0.0
y_0 = 0.0
u = (x, y) -> exp(-((x-x_0)^2 + (y-y_0)^2)/(2*σ^2))

P = (x, y) -> -(x-x_0) / σ^2 * u(x, y)
Q = (x, y) -> -(y-y_0) / σ^2 * u(x, y)
Δu = (x, y) -> (((x-x_0)^2 + (y-y_0)^2)/σ^4 - 2 / σ^2) * u(x, y)
function func(x)
    return Δu(x...)
end

# START CODE #

degree = 2 # can be 0 through 10
qorder = Inti.TRIANGLE_VR_IORDER_TO_QORDER[degree]
num_points_per_quad = Inti.TRIANGLE_VR_ORDER_TO_NPTS[qorder]
@show num_points_per_quad
qrule = Inti.TRIANGLE_VR_QRULES[num_points_per_quad]
weights = []

degs = sort([(i, j) for i in 0:degree, j in 0:degree if i+j <= degree], by = sum)
len_vander = Int64((degree+1)*(degree+2) / 2)

# create reference vandermonde matrix
# mat = zeros(num_points_per_quad, len_vander)
# for (i, node) in enumerate(qrule)
#     pnt = node[1]
#     for (j, deg) in enumerate(degs)
#         mat[i, j] = (1/factorial(deg[1]))*pnt[1]^deg[1] * (1/factorial(deg[2]))*pnt[2]^deg[2]
#     end
#     weight = node[2]
#     push!(weights, weight)
# end

# A = mat
# F = qr(Diagonal(sqrt.(w))A)
# Qw = Matrix(F.Q)[:, 1:size(mat, 2)]
# R = UpperTriangular(F.R)

gmsh.initialize()
gmsh.option.setNumber("General.Verbosity", 3)
gmsh.clear()
parametrizations::Vector{Function} = [(x) -> (cos(2*x*pi), sin(2*x*pi))]
splines = createSplines(parametrizations)
curve = gmsh.model.occ.addCurveLoop(splines)
surface = gmsh.model.occ.addPlaneSurface([curve])
gmsh.model.occ.synchronize()
gmsh.model.addPhysicalGroup(1, splines, -1, "Boundary")
gmsh.model.mesh.removeDuplicateNodes()
size_field = -1
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 30)

for iter in 1:15
    gmsh.model.mesh.clear()
    gmsh.model.mesh.generate(2)
    mesh = Inti.import_mesh(; dim=2)
    # showMesh(mesh)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    data = []
    dom = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, mesh)
    dom_mesh = view(mesh, dom)
    
    quadrature = Inti.Quadrature(dom_mesh; qorder = qorder)
    
    c1, c2 = 0, 0
    conds = []
    vertex_dict = Dict()
    triangles = Inti.elements(mesh, Inti.LagrangeElement{Inti.ReferenceSimplex{2}, 3, StaticArraysCore.SVector{2, Float64}})
    for (i, elem) in enumerate(triangles)
        v1, v2, v3 = ((round(v[1], digits=8), round(v[2], digits=8)) for v in Inti.vertices(elem))
        quad_points = quadrature[(i-1)*num_points_per_quad+1:i*num_points_per_quad]
        if iter > 1
            center = Inti.center(elem)
            target, hi = gmsh.view.probe(size_field, center[1], center[2], 0.0)
            size = target[1]
        else
            size = maximum((norm(v1.-v2), norm(v2.-v3), norm(v1.-v3)))
        end
        refine, cond = refinementTriangleStage(elem, quad_points, func, degree)
        push!(conds, cond)
        if refine
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

    @show triangles[argmax(conds)]
    @show maximum(conds)

    global size_field = gmsh.view.add("mesh size field")
    gmsh.view.addListData(size_field, "ST", length(triangles), data)

    bg_field = gmsh.model.mesh.field.add("PostView")
    gmsh.model.mesh.field.setNumber(bg_field, "ViewTag", size_field)
    gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)

    cfs = zeros(len_vander, 1)
    point = SVector{2, Float64}(0.0, -0.18)
    println("-----------------------------------------------------------------")
    counter = 0
    for (i, elem) in enumerate(triangles)
        verts = Inti.vertices(elem)
        if any(norm(vert - point) < 0.025 for vert in verts)
            counter += 1
            if counter > 10
                break
            end
            quad_points = quadrature[(i-1)*num_points_per_quad+1:i*num_points_per_quad]
            d = func.(Inti.coords.(quad_points))
            w = Inti.weight.(quad_points)
            c, r = Inti.translation_and_scaling(elem)

            d = func.(Inti.coords.(quad_points))
            mat = zeros(num_points_per_quad, len_vander)
            for (i, point) in enumerate(quad_points)
                pnt = 1/r * (point.coords .- c)
                for (j, deg) in enumerate(degs)
                    mat[i, j] = (1/factorial(deg[1]))*pnt[1]^deg[1] * (1/factorial(deg[2]))*pnt[2]^deg[2]
                end
            end
            cfs += mat \ d
            # @show verts
            # @show r
            for deg in 0:degree
                sum = 0
                for i in 1:length(degs)
                    if degs[i][1]+degs[i][2] == deg
                        sum += cfs[i]
                    end
                end
                sum /= deg+1
                # println(deg, " : ", sum)
            end
            # println("------------------")
        end
    end
    cfs = cfs ./ length(triangles)
    
    # @show cfs
    # @show degs
    # @show maximum(conds)
    # @show median(conds)
end
mesh = Inti.import_mesh(; dim=2)
showMesh(mesh)
println(mesh)

gmsh.finalize()