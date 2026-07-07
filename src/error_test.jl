using Inti
using Gmsh
using LinearAlgebra
using Meshes
using CairoMakie

function distance(a, b)
    return sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)
end

function getBoundaryPoints(parametrization::Function, num_points::Int64)::Vector{Tuple{Float64, Float64}}
    boundary = Vector{Tuple{Float64, Float64}}()
    for t in range(0, stop=1, length=num_points+1)[1:end-1]
        point = parametrization(t)
        push!(boundary, point)
    end
    return boundary
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

function showSeparatedMesh(mesh, quad, target, other_points)
    Ω = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    Γ = Inti.boundary(Ω)
    Ω_msh = @views mesh[Ω]
    Γ_msh = @views mesh[Γ]
    fig = viz(
        Ω_msh;
        segmentsize = 1,
        showsegments = true,
        axis = (aspect = DataAspect(),),
        figure = (; size = (1000, 1000)),
    )
    viz!(Γ_msh; color = :red, segmentsize = 2)

    domain_quad_coords = Meshes.PointSet(map((q) -> Meshes.Point(q.coords...), quad))
    viz!(domain_quad_coords; color = :blue, pointsize = 5)
    viz!(Meshes.Point(target); color = :green, pointsize = 20)

    for pnt in other_points
        viz!(Meshes.Point(pnt[1], pnt[2]); color = :red, pointsize = 10)
    end
    display(fig)
end

op = Inti.Laplace(; dim=2)

# compute exact solution
u = (x) -> x[1]^2 + x[2]^2
du = (x, normal) -> dot((2*x[1], 2*x[2]), normal)

laplacian_u = (x) -> 4.0

# now, create the mesh         
parametrizations::Vector{Function} = [(x) -> (cos(2*x*pi), sin(2*x*pi))]

# for offset in [0.1*t for t in -4:4]
for offset in [0.0]
    parametrizations2::Vector{Function} = [(t) -> (0.25*t, 0.0), 
        (t) -> (0.25+offset*t, 0.25*t), 
        (t) -> ((0.25+offset) + (0.4-offset)*t, 0.25-0.25*t), 
        (t) -> (0.65-0.4*t, -0.25*t), 
        (t) -> (0.25-0.25*t, -0.25+0.25*t)
    ]
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 1) 
    gmsh.clear()
    splines = createSplines(parametrizations)
    curve = gmsh.model.occ.addCurveLoop(splines)
    splines2 = createSplines(parametrizations2)
    curve2 = gmsh.model.occ.addCurveLoop(splines2)
    surface = gmsh.model.occ.addPlaneSurface([curve, curve2])
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(1, splines, -1, "Boundary")
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 30)
    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.mesh.generate(2)
    mesh = Inti.import_mesh(; dim=2)
    gmsh.finalize()

    # create quadratures + get target points
    domain = Inti.Domain((e) -> Inti.geometric_dimension(e) == 2, mesh)
    domain_mesh = Inti.view(mesh, domain)
    domain_quadrature = Inti.Quadrature(domain_mesh; qorder = 4)

    boundary = Inti.boundary(domain)
    boundary_mesh = Inti.view(mesh, boundary)
    boundary_quadrature = Inti.Quadrature(boundary_mesh; qorder = 6)

    # now compute volume potential with Inti routine
    # target = [(q.coords[1], q.coords[2]) for q in domain_quadrature]
    target = [(-Inf, -Inf)]
    goal_point = (0.25, 0.0)
    for q in domain_quadrature
        if distance(q.coords, goal_point) < distance(target[1], goal_point)
            target[1] = (q.coords[1], q.coords[2])
        end
    end

    # target = [(0.23017316817953776, 0.008379805131740959)]
    multiplicative_terms = [:inside for i in 1:length(target)]

    inside_target_points = Vector{Tuple{Float64, Float64}}()
    inside_indices = Vector{Int64}()
    for (i, point) in enumerate(target)
        if multiplicative_terms[i] == :inside
            push!(inside_target_points, point)
            push!(inside_indices, i)
        end
    end

    # set max dist to Inf, it doesn't matter (always gives back one triangle)
    max_dist = Inf
    potentials = [0.0 for x in target]
    laplacian_points = [laplacian_u(q.coords) for q in domain_quadrature]

    inside_volume_potential = Inti.volume_potential(; 
        op, 
        target = inside_target_points, 
        source = domain_quadrature,
        compression = (method = :none, ),
        correction = (method = :dim, maxdist=max_dist, target_location=:inside), 
    )
    inside_volume_potential_2 = Inti.volume_potential(; 
        op, 
        target = inside_target_points, 
        source = domain_quadrature,
        compression = (method = :none, ),
        correction = (method = :none, ),
    )

    inside_volume_potential_mat = Matrix(inside_volume_potential - inside_volume_potential_2)
    diff_coords = []
    for (i, q) in enumerate(domain_quadrature)
        if inside_volume_potential_mat[i] != 0.0
            push!(diff_coords, Vector(q.coords))
        end
    end
    println(diff_coords)

    inside_potentials = inside_volume_potential * laplacian_points
    for i in 1:length(inside_potentials)
        potentials[inside_indices[i]] = -inside_potentials[i] 
    end

    S, D = Inti.single_double_layer(;
        op,
        target,
        source = boundary_quadrature,
        compression = (method = :none, ), 
        correction = (method = :dim, target_location = :inside, maxdist = 0.4)
    )
    γ₀u = map(q -> u(q.coords), boundary_quadrature)
    γ₁u = map(q -> du(q.coords, q.normal), boundary_quadrature)
    potentials += S*γ₁u - D*γ₀u

    errors = abs.(potentials - u.(target))

    if offset == 0.0
        println("THIS IS THE IMPORTANT CASE")
    end
    println("Error: ", errors[1])
    println("---------")

    showSeparatedMesh(mesh, domain_quadrature, target[1], diff_coords)
end