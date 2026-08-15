With a few improvements, I have gotten the time for quad_contribution down from ~7 seconds to ~2.5 seconds. The memory allocation is also down significantly, from 7-8 gb down to ~1.3 gb.

However, the runtime and the memory allocation are still both entirely dominated by the calculation of the far-field corrections for the near-singular terms:
- For the rest of the function: 50-150 allocations, 0.000158 seconds per iteration
- For the correction calculations: 95117 allocations, 0.002442 seconds per iteration

so this needs to be improved more. I've spent a while optimizing this already but I'll work on it more

Some other @timev stuff:
Nothing in the loop (so all other applications):
0.356918 seconds (286.38 k allocations: 26.838 MiB, 55.12% compilation time)
elapsed time (ns):  3.56918177e8
gc time (ns):       0
bytes allocated:    28141352
pool allocs:        286010
non-pool GC allocs: 2
malloc() calls:     367
free() calls:       0
minor collections:  0
full collections:   0

coords = [(x[1], x[2]) for x in Inti.vertices(elem)] adds 0.05 seconds to the runtime??? tho memory addition is fine

near_target_points = NearestNeighbors.inrange(tree, center, h*1.5)[1] adds nearly 1600 malloc calls
0.458701 seconds (431.33 k allocations: 47.733 MiB, 64.17% compilation time)
elapsed time (ns):  4.58700605e8
gc time (ns):       0
bytes allocated:    50051432
pool allocs:        429306
non-pool GC allocs: 2
malloc() calls:     2023
free() calls:       0
minor collections:  0
full collections:   0

filter_target_points(filtered_target_points, near_target_points, target, coords, h, L_map) also nearly adds 0.05 seconds, another 1600 mallocs
0.492556 seconds (447.16 k allocations: 60.702 MiB, 63.68% compilation time)
elapsed time (ns):  4.92555618e8
gc time (ns):       0
bytes allocated:    63650304
pool allocs:        443522
non-pool GC allocs: 2
malloc() calls:     3632
free() calls:       0
minor collections:  0
full collections:   0

this is with everything except for the correction
0.612254 seconds (494.94 k allocations: 67.429 MiB, 2.86% gc time, 57.03% compilation time)
elapsed time (ns):  6.12253869e8
gc time (ns):       17494455
bytes allocated:    70704568
pool allocs:        490388
non-pool GC allocs: 912
malloc() calls:     3637
free() calls:       1715
minor collections:  1
full collections:   0

this is with the correction
3.021611 seconds (87.53 M allocations: 1.363 GiB, 5.35% gc time, 12.38% compilation time)
elapsed time (ns):  3.021611447e9
gc time (ns):       161803827
bytes allocated:    1463479608
pool allocs:        87529657
non-pool GC allocs: 912
malloc() calls:     3639
free() calls:       3812
minor collections:  16
full collections:   0

key findings currently: when looking at second_corrections, the internal loop allocated on average 2-2.5 times per iter

second_correction, kernel, target, p, elem_q_points, elem_q_weights, f_vals, k, j
correction_contribution!(
    StaticArraysCore.MArray{Tuple{224}, Float64, 1, 224},  - second_correction
    AdaptiveMeshSolver.var"#62#63",  - kernel
    Array{StaticArraysCore.SArray{Tuple{2}, Float64, 1, 2}, 1},  - target
    Int64,  - p
    Array{StaticArraysCore.SArray{Tuple{2}, Float64, 1, 2}, 1},  - elem_q_points
    Array{Float64, 1}, - elem_q_weights
    Array{Float64, 1}, - f_vals
    Int64, - k
    Int64  - j
)

target uses gc_small_alloc
elem_q_points does too
elem_q_weights does three
f_vals does four

so running getindex

FOUND ISSUE: it was the kernel function not being defined locally??? very strange, temporarily just using the kernel func here

Average runtime / allocations for certain parts of the program:
Application, no quad near singular correction: ~0.13 seconds, 216 heap allocations
Application, no problematic correction: ~0.24 seconds, 76.80k heap allocations
Full application: ~0.91 seconds, 78.62k heap allocations

Each map time:
FMM map - ~0.08 seconds, 76 allocations
quad -> triangle map - ~0.002 seconds, 3 allocations
quad -> quad map - ~0.81 seconds, 78.41k allocations
triangle -> triangle map - 0.047 seconds, 140 allocations

SO quad -> quad is STILL BAD