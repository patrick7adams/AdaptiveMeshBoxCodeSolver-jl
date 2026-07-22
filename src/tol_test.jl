using AdaptiveMeshSolver
using Inti
using LinearAlgebra
using Printf
degs = 4:8
tols = [10.0^(-x) for x in 4:8]

function greens_thm_stuff()
    for σ in [0.5, 0.25, 0.125]
        vals = []
        for degree in degs
            println(degree)
            tol_vals = []
            for tol in tols
                x_0 = 0.85
                y_0 = 0.0
                u = (x, y) -> exp(-((x-x_0)^2 + (y-y_0)^2)/(2*σ^2))

                P = (x, y) -> -(x-x_0) / σ^2 * u(x, y)
                Q = (x, y) -> -(y-y_0) / σ^2 * u(x, y)
                Δu = (x, y) -> (((x-x_0)^2 + (y-y_0)^2)/σ^4 - 2 / σ^2) * u(x, y)

                boundary_fn = (q) -> P(q.coords...)*q.normal[1] + Q(q.coords...)*q.normal[2]
                domain_fn = (q) -> Δu(q.coords...)

                forcing_function = (x) -> Δu(x...)

                parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
                meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, forcing_function, false, degree, tol)
                # AdaptiveMeshSolver.showMeshes(meshes)

                domain_quads = []
                boundary_quads = []
                for (i, mesh) in enumerate(meshes)
                    push!(domain_quads, AdaptiveMeshSolver.getDomainQuadrature(mesh, 12))
                    if i != 1
                        push!(boundary_quads, AdaptiveMeshSolver.getBoundaryQuadrature(mesh, 6)) # order 6 is fine
                    end
                end
                domain_integral = 0.0
                for quad in domain_quads
                    domain_integral += Inti.integrate(domain_fn, quad)
                end
                boundary_integral = 0.0
                for quad in boundary_quads
                    boundary_integral += Inti.integrate(boundary_fn, quad)
                end
                calculated_u_val = abs(domain_integral - boundary_integral)
                push!(tol_vals, calculated_u_val)
            end
            push!(vals, tol_vals)
        end
        println("Accuracy matrix for σ=",σ, ":\n")
        @printf("Tolerances ")
        for j in 1:length(vals[1])
            @printf("| %1.3e ", tols[j])
        end
        st = "-"^(12+length(vals[1])*12)
        @printf("|\n")
        println(st)
        for i in 1:length(vals)
            @printf("Degree %2d  ", degs[i])
            for j in 1:length(vals[1])
                @printf("| %1.3e ", vals[i][j])
            end
            @printf("|\n")
        end
        @printf("\n\n\n")
    end
end

function greens_identity_test()
    for σ in [0.5, 0.25, 0.125]
        vals = []
        for degree in degs
            println(degree)
            tol_vals = []
            for tol in tols
                op = Inti.Laplace(; dim=2)
                x_0 = 0.85
                y_0 = 0.0

                # compute exact solution
                u = (x) -> exp(-((x[1]-x_0)^2 + (x[2]-y_0)^2)/(2*σ^2))
                du = (x, normal) -> dot((-(x[1]-x_0)/σ^2 * u(x), -(x[2]-y_0)/σ^2 * u(x)), normal)

                laplacian_u = (x) -> (((x[1]-x_0)^2 + (x[2]-y_0)^2)/σ^4 - 2 / σ^2) * u(x)

                # now create meshes
                parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]
                meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, laplacian_u, false, degree, tol)
                domain_quadratures = [AdaptiveMeshSolver.getDomainQuadrature(mesh, 4) for mesh in meshes]
                boundary_quadratures = [AdaptiveMeshSolver.getBoundaryQuadrature(mesh, 6) for mesh in meshes[2:end]]
                AdaptiveMeshSolver.showMeshes(meshes)
                
                target = []
                multiplicative_terms = []
                for (i, quadrature) in enumerate(domain_quadratures)
                    points = [(q.coords[1], q.coords[2]) for q in quadrature]
                    append!(target, points)
                end
                println(length(target))
                
                count = 1
                for (i, quadrature) in enumerate(domain_quadratures)
                    l = length(quadrature)
                    terms = [i >= count && i < count+l ? :inside : :outside for i in 1:length(target)]
                    push!(multiplicative_terms, terms)
                    count += l
                end
                
                potentials = AdaptiveMeshSolver.calculateVolumePotential(domain_quadratures, meshes, laplacian_u, target, multiplicative_terms)

                for quad in boundary_quadratures
                    S, D = Inti.single_double_layer(;
                        op,
                        target,
                        source = quad,
                        compression = (method = :none, ), 
                        correction = (method = :dim, target_location = :inside, maxdist = 0.4)
                    )

                    γ₀u = map(q -> u(q.coords), quad)
                    γ₁u = map(q -> du(q.coords, q.normal), quad)
                    
                    contribution = S*γ₁u - D*γ₀u
                    potentials += contribution
                end
                errors = abs.(potentials - u.(target))
                println(maximum(errors))
                
                push!(tol_vals, maximum(errors))
            end
            push!(vals, tol_vals)
        end
        println("Accuracy matrix for σ=",σ, ":\n")
        @printf("Tolerances ")
        for j in 1:length(vals[1])
            @printf("| %1.3e ", tols[j])
        end
        st = "-"^(12+length(vals[1])*12)
        @printf("|\n")
        println(st)
        for i in 1:length(vals)
            @printf("Degree %2d  ", degs[i])
            for j in 1:length(vals[1])
                @printf("| %1.3e ", vals[i][j])
            end
            @printf("|\n")
        end
        @printf("\n\n\n")
    end
end

greens_identity_test()