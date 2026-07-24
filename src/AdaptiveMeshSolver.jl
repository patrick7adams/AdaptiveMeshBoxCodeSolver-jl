module AdaptiveMeshSolver

using P4estTypes;
using NearestNeighbors
using Polynomials
using FastChebInterp
using Inti, Meshes, CairoMakie, StaticArrays, Gmsh
using PolygonAlgorithms
using GeometryBasics
using Colorfy
using FMM2D
using IterativeSolvers, LinearAlgebra
using Statistics
using ClassicalOrthogonalPolynomials, MultivariateSingularIntegrals

include("mesh_tools.jl")

include("main.jl")

end