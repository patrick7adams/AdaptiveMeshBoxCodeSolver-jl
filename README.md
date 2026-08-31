# AdaptiveMeshBoxCodeSolver-jl

A package for efficient adaptive mesh construction and volume potential application. This package will be combined with Inti.jl soon.

## Usage
First, we perform some setup. We will show the usage of this package by utilizing Green's Third Identity:
$$ u(x) = \int_{\partial V} \left( G(x, y) \frac{\partial u}{\partial n}(y) - u(y) \frac{\partial G}{\partial n}(x, y)\right)dS(y) + \int_V G(x, y) \Delta u(y) dV(y)$$

We will let $u$ be a strong shifted gaussian, and through Green's Third Identity, we will reconstruct $u$ from $\frac{\partial u}{\partial n}$ and $\Delta u$.
```julia
op = Inti.Laplace(; dim=2)

σ = 0.06
x_0 = 0.5
y_0 = 0.0

# create function
u = (x) -> exp(-((x[1]-x_0)^2 + (x[2]-y_0)^2)/(2*σ^2))
du = (x, normal) -> dot((-(x[1]-x_0)/σ^2 * u(x), -(x[2]-y_0)/σ^2 * u(x)), normal)
laplacian_u = (x) -> (((x[1]-x_0)^2 + (x[2]-y_0)^2)/σ^4 - 2 / σ^2) * u(x)

# parametrize function boundary
parametrizations::Vector{Vector{Function}} = [[(x) -> (cos(x*2*pi), sin(x*2*pi))]]

# create the adaptive mesh
meshes = AdaptiveMeshSolver.createQuadtreeMesh(parametrizations, laplacian_u, false)

# now create quadratures and function vals
domain_quadratures = AdaptiveMeshSolver.AdaptiveQuadrature(meshes, 8, 12)
boundary_quadratures = [AdaptiveMeshSolver.getBoundaryQuadrature(mesh, 18) for mesh in meshes[2:end]]
f_vals = stack(laplacian_u.(AdaptiveMeshSolver.target_points(domain_quadratures)))

# finally create volume potential operator
vol_pot = AdaptiveMeshSolver.adaptive_volume_potential(; 
    op=op, 
    source=domain_quadratures, 
    compression=(method = :fmm, tol=1e-14)
)       
# apply operator   
@potentials = vol_pot * f_vals


# we can create target point vector
target = Vector{SVector{2, Float64}}()
for (i, quadrature) in enumerate(domain_quadratures)
    points = [SVector{2, Float64}(q.coords[1], q.coords[2]) for q in quadrature]
    append!(target, points)
end
# then use this to compute boundary contributions
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
# finally, compare reconstructed solution with u
errors = abs.(potentials - u.(target))
```

Now, we can look at the error and application time:
```julia
julia> @show maximum(errors):
> maximum(errors) = 4.935651540249886e-11

julia> @time potentials = vol_pot * f_vals
> 0.308969 seconds (221 allocations: 4.627 MiB)
```