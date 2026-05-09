using P4estTypes;
using Polynomials
using FastChebInterp
using Inti, Meshes, CairoMakie, StaticArrays, Gmsh
using PolygonAlgorithms
using GeometryBasics
using Colorfy
using FMM2D
using IterativeSolvers, LinearAlgebra

meshsize = 0.25
max_triangle_size = 0.25

function showMesh(msh)
    println("Showing mesh!!!")
    Ω = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, msh)
    Γ = Inti.boundary(Ω)
    Ω_msh = @views msh[Ω]
    Γ_msh = @views msh[Γ]
    fig = viz(
        Ω_msh;
        segmentsize = 1,
        showsegments = true,
        axis = (aspect = DataAspect(),),
        figure = (; size = (500, 400)),
    )
    viz!(Γ_msh; color = :red, segmentsize = 1)
    display(fig)
end

function calculate_area(mesh, expected)
    dom = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, mesh)
    quad = Inti.Quadrature(mesh[dom]; qorder = 10)
    area = Inti.integrate(x -> 1, quad)
    println("Error in area: ", area - expected)
end

gmsh.initialize()
gmsh.option.setNumber("General.Verbosity", 3) # turn to 4/5 for info or debug, 3 is all that is necessary I think though

circle = Inti.gmsh_curve(0, 2pi;meshsize) do s
    return Inti.Point2D(cos(s), sin(s))
end
curve_loop = gmsh.model.occ.addCurveLoop([circle])

circle2 = Inti.gmsh_curve(0, 2pi;meshsize) do s
    return Inti.Point2D(0.2*cos(s)+0.5, 0.1*sin(s)+0.5)
end
inner_loop2 = gmsh.model.occ.addCurveLoop([circle2])

inscribedPoints = [(0, 0), (0.5, 0.05), (0, 0.1), (-0.05, 0.4), (-0.1, 0.1), (-0.4, 0.05), (-0.1, 0), (-0.05, -0.4), (0, 0)]

point_tags = []
for point in inscribedPoints
    tag = gmsh.model.occ.addPoint(point[1], point[2], 0, max_triangle_size)
    push!(point_tags, tag)
end
curve_tags = []
for i in 1:length(point_tags)-1
    tag = gmsh.model.occ.addLine(point_tags[i], point_tags[i+1])
    push!(curve_tags, tag)
end
inner_loop = gmsh.model.occ.addCurveLoop(curve_tags)
surf = gmsh.model.occ.addPlaneSurface([curve_loop, inner_loop2])
gmsh.model.occ.synchronize()



gmsh.model.addPhysicalGroup(1, [circle], -1, "curve_loop")
gmsh.model.addPhysicalGroup(1, [circle2], -1, "inner_loop")
gmsh.model.addPhysicalGroup(1, curve_tags, -1, "inner_loop2")

gmsh.option.setNumber("Mesh.MeshSizeMax", max_triangle_size)
gmsh.model.mesh.generate(2)
gmsh.model.mesh.setOrder(2)
# gmsh.write("testing.msh")

mesh = Inti.import_mesh(; dim=2)
# showMesh(mesh)
star_area = 0.085
oval_area = 0.1*0.2*pi
area = pi - star_area - oval_area

dom = Inti.Domain(Inti.entities(mesh)) do ent
    return Inti.geometric_dimension(ent) == 2
end

for element in Inti.elements(view(mesh, dom))
    cent = Inti.center(element)
    # now associate center with the closest boundary
end

gmsh.finalize()

parameterization_1 = (t) -> [cos(2*pi*t), sin(2*pi*t)]
parameterization_2 = (t) -> [0.2*cos(2*pi*t)+0.5, 0.1*sin(2*pi*t)+0.5]
entity_parameterizations = Dict{Inti.EntityKey, Function}()
for e in Inti.all_keys(dom)
    # println(e)
    l = Inti.labels(e)
    if "curve_loop" in l
        entity_parameterizations[e] = parameterization_1
    elseif "inner_loop2" in l
        entity_parameterizations[e] = parameterization_2
    end
    # no need to look at inner_loop, its already flat
end
println(entity_parameterizations)
theta = 6
curve_mesh = Inti.curve_mesh(msh, parameterization_1, theta)
showMesh(curve_mesh)



# gmsh.initialize()
# meshsize = 0.075
# gmsh.option.setNumber("Mesh.MeshSizeMax", meshsize)
# gmsh.option.setNumber("Mesh.MeshSizeMin", meshsize)

# # Three circles
# c1 = gmsh.model.occ.addDisk(0, 0, 0, 1, 1)
# c2 = gmsh.model.occ.addDisk(0, 3.0, 0, 1, 1)
# c3 = gmsh.model.occ.addDisk(0, 8.0, 0, 2, 2)
# gmsh.model.occ.synchronize()

# # Add tags for stable identification of the entities
# gmsh.model.addPhysicalGroup(2, [c1], -1, "c1")
# gmsh.model.addPhysicalGroup(2, [c2], -1, "c2")
# gmsh.model.addPhysicalGroup(2, [c3], -1, "c3")

# gmsh.model.mesh.generate(2)
# msh = Inti.import_mesh(; dim = 2)

