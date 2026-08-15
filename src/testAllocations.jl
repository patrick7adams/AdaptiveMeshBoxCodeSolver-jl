using AdaptiveMeshSolver, Inti

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
# error("bruh")
domain_quadratures = AdaptiveMeshSolver.AdaptiveQuadrature(meshes, 8, 17)
boundary_quadratures = [AdaptiveMeshSolver.getBoundaryQuadrature(mesh, 18) for mesh in meshes[2:end]]
# AdaptiveMeshSolver.showMeshes(meshes)
# println(meshes[2])
# AdaptiveMeshSolver.showMesh(meshes[1], [[-0.18906250000000005, -0.30550088025258804], [-0.20625000000000004, -0.48125000000000007]])
# error("HI")
f_vals = stack(laplacian_u.(AdaptiveMeshSolver.target_points(domain_quadratures)))

vol_pot = AdaptiveMeshSolver.adaptive_volume_potential(; 
    op=op, 
    source=domain_quadratures, 
    compression=(method = :fmm, tol=1e-14)
)
potentials = vol_pot * f_vals