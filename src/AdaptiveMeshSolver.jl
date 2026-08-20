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
using Serialization
using ClassicalOrthogonalPolynomials, MultivariateSingularIntegrals
using NearestNeighbors
using PolynomialBases
using LinearMaps
using SparseArrays
using InteractiveUtils
using BenchmarkTools
using Cthulhu
using Infiltrator
using Profile, PProf

include("mesh_tools.jl")

include("main.jl")

end