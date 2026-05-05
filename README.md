# AdaptiveMeshBoxCodeSolver-jl
Julia adaptive mesh box code solver

## Current Unfixed Bugs
Issues when boundaries are too close to each other without sufficient refinement. Not intending to fix for a second
Current issue: certain functions cause infinite refinement
- for some reason, the coefficient matrix of the chebyshev polynomial for some of these points is 1x1 containing only Inf
