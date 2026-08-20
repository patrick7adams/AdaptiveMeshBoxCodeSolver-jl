using Revise
using AdaptiveMeshSolver, Inti, StaticArrays, LinearAlgebra

function run_stuff()
    σ = 0.06
    op = Inti.Laplace(; dim=2)
    x_0 = 0.5
    y_0 = 0.0

    # compute exact solution
    u = (x) -> exp(-((x[1]-x_0)^2 + (x[2]-y_0)^2)/(2*σ^2))
    du = (x, normal) -> dot((-(x[1]-x_0)/σ^2 * u(x), -(x[2]-y_0)/σ^2 * u(x)), normal)

    laplacian_u = (x) -> (((x[1]-x_0)^2 + (x[2]-y_0)^2)/σ^4 - 2 / σ^2) * u(x)

    # now create meshes
    parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
    meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, laplacian_u, false)
    println("Generated mesh!")
    domain_quadratures = AdaptiveMeshSolver.AdaptiveQuadrature(meshes, 8, 17)
    boundary_quadratures = [AdaptiveMeshSolver.getBoundaryQuadrature(mesh, 18) for mesh in meshes[2:end]]
    f_vals = stack(laplacian_u.(AdaptiveMeshSolver.target_points(domain_quadratures)))

    vol_pot = AdaptiveMeshSolver.adaptive_volume_potential(; 
        op=op, 
        source=domain_quadratures, 
        compression=(method = :fmm, tol=1e-14)
    )            
    potentials = vol_pot * f_vals

    target = Vector{SVector{2, Float64}}()
    multiplicative_terms = []
    for (i, quadrature) in enumerate(domain_quadratures)
        points = [SVector{2, Float64}(q.coords[1], q.coords[2]) for q in quadrature]
        append!(target, points)
    end
    
    for quad in boundary_quadratures
        S, D = Inti.single_double_layer(;
            op,
            target,
            source = quad,
            compression = (method = :none, ), 
            correction = (method = :dim, target_location = :inside, maxdist = Inf)
        )

        γ₀u = map(q -> u(q.coords), quad)
        γ₁u = map(q -> du(q.coords, q.normal), quad)
        
        contribution = S*γ₁u - D*γ₀u
        potentials += contribution
    end
    errors = abs.(potentials - u.(target))
    n = length(errors)
    L2_error = sqrt(1/n * sum(errors.^2))
    max_error = maximum(errors)
    @show L2_error
    @show max_error
    AdaptiveMeshSolver.showErrorMesh(meshes, target, errors)
end
run_stuff()