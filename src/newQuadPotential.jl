function calculateQuadVolumePotential(meshes, quadrature, u, target_points)
    quad_mesh = meshes[1]
    greens_fn = (x, y) -> 1/(2pi) * log(distance(x, y))
    potentials = [0.0 for point in target_points]
    deg = 5;
    F_map = Dict{Inti.LagrangeElement{Inti.ReferenceHyperCube{2}, 4, StaticArraysCore.SVector{2, Float64}}, Matrix{Float64}}()
    L_map = Dict{Vector{Float64}, Matrix{Float64}}()
    Correction_map = getCorrectionMap(quad_mesh, quadrature, u)

    # first set up legendre polynomials for singular evaluations
    P = ClassicalOrthogonalPolynomials.Legendre()
    
    x = ClassicalOrthogonalPolynomials.grid(P, deg)

    # construct kdtree
    tree = NearestNeighbors.KDTree(target_points)

    # map quad mesh element -> iterable of quadrature points over that element
    elem_to_quad_points = Dict()
    num_points_per_quad = Int64(ceil((order + 1) / 2))^2
    for (i, elem) in enumerate(Inti.elements(quad_mesh))
        elem_to_quad_points[elem] = quadrature[(i-1)*num_points_per_quad+1:i*num_points_per_quad]
    end

    for elem in Inti.elements(quad_mesh)
        # loop through elements, get target points within range of center
        coords = [(x[1], x[2]) for x in Inti.vertices(elem)]
        bounds = ((coords[1][1], coords[3][1]), (coords[1][2], coords[3][2]))
        h = coords[3][1] - coords[1][1]
        center = ((coords[3][1] + coords[1][1])/2, (coords[3][2] + coords[1][2])/2)
        near_target_points = NearestNeighbors.inrange(tree, center, h*1.25)
        elem_quadrature_points = elem_to_quad_points[elem]
        if haskey(F_map, elem)
            F = F_map[elem]
        else
            ux = generatePoints(u, x, bounds)'
            tmp_F = ClassicalOrthogonalPolynomials.plan_transform(P, (deg, deg))
            F = (tmp_F * ux)
            F_map[elem] = F
        end

        # add in all potentials
        for (i, point) in enumerate(target_points)
            domain_func = (q) -> greens_fn(point, q.coords) * u(q.coords)
            potentials[i] += Inti.integrate(domain_func, quadrature)
        end
        
        correction_term = log(2/h)/(2*pi)*Correction_map[elem]

        # now calculate both singular method and nonsingular method quad contributions for each of these points
        for point_idx in near_target_points
            point = target_points[point_idx]
            # first singular method
            x_dist, y_dist = point[1] - coords[1][1], point[2] - coords[1][2]
            
            mapped_target_point = [round(-1.0+2/h*x_dist, digits=12), round(-1.0+2/h*y_dist, digits=12)]
            mapped_target_point += [1e-16, 1e-16]
            # write this statically instead
            if haskey(L_map, mapped_target_point)
                L = L_map[mapped_target_point]
            else
                L = Float64.(MultivariateSingularIntegrals.newtoniansquare(big.(mapped_target_point), deg))
                L_map[mapped_target_point] = L
            end
            I = dot(L, F)
            term = I*h^2 / (8*pi) - correction_term
            potentials[point_idx] += term

            # now just subtract this quad's normal quadrature contribution from the estimated potential
            domain_func = (q) -> greens_fn(point, q.coords) * u(q.coords)
            potentials[point_idx] -= sum(domain_func(q)*q.weight for q in elem_quadrature_points)
        end
    end
    return potentials
end