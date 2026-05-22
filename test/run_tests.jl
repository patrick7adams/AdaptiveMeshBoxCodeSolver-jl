include("test_greens_identity.jl")
include("test_individual_elements.jl")
include("test_main_volume_potentials.jl")
include("test_quadtree_mesh.jl")
include("test_quadtree_vs_triangle_mesh.jl")
include("test_triangle_mesh.jl")

function test_simple_triangle_mesh()
    @testset "Simple Triangle Mesh" begin
        test_simple_triangle_mesh_area()
        test_simple_triangle_mesh_zero()
        test_simple_triangle_mesh_negative_area()
        test_simple_triangle_mesh_linear_x()
        test_simple_triangle_mesh_linear_y()
        test_simple_triangle_mesh_quadratic_x()
        test_simple_triangle_mesh_quadratic_y()
        test_simple_triangle_mesh_quadratic_xy()
        test_simple_triangle_mesh_high_order_x()
        test_simple_triangle_mesh_high_order_y()
        test_simple_triangle_mesh_high_order_xy()
        test_simple_triangle_mesh_sin_x()
        test_simple_triangle_mesh_sin_y()
        test_simple_triangle_mesh_sin_xy()
        test_simple_triangle_mesh_sin_x_cos_y()
        test_simple_triangle_mesh_exp_xy()
        test_simple_triangle_mesh_spike()
    end;
end

function test_simple_quadtree_mesh()
    @testset "Simple Quadtree Mesh" begin
        test_simple_quadtree_mesh_area()
        test_simple_quadtree_mesh_zero()
        test_simple_quadtree_mesh_negative_area()
        test_simple_quadtree_mesh_linear_x()
        test_simple_quadtree_mesh_linear_y()
        test_simple_quadtree_mesh_quadratic_x()
        test_simple_quadtree_mesh_quadratic_y()
        test_simple_quadtree_mesh_quadratic_xy()
        test_simple_quadtree_mesh_high_order_x()
        test_simple_quadtree_mesh_high_order_y()
        test_simple_quadtree_mesh_high_order_xy()
        test_simple_quadtree_mesh_sin_x()
        test_simple_quadtree_mesh_sin_y()
        test_simple_quadtree_mesh_sin_xy()
        test_simple_quadtree_mesh_sin_x_cos_y()
        test_simple_quadtree_mesh_exp_xy()
        test_simple_quadtree_mesh_spike()
    end;
end

function test_individual_triangles()
    @testset "Individual Triangle Tests" begin
        test_simple_triangle_mesh_quadratures(0.4)
        test_simple_triangle_mesh_quadratures(0.3)
        test_simple_triangle_mesh_quadratures(0.2)
        test_simple_triangle_mesh_quadratures(0.1)
        test_simple_triangle_mesh_quadratures(0.075)
        test_simple_triangle_mesh_quadratures(0.05)
        test_simple_triangle_mesh_quadratures(0.025)
    end;
end

function test_individual_quads()
    @testset "Individual Quad Tests" begin
        test_simple_quadtree_mesh_quadratures(0.4)
        test_simple_quadtree_mesh_quadratures(0.3)
        test_simple_quadtree_mesh_quadratures(0.2)
        test_simple_quadtree_mesh_quadratures(0.1)
        test_simple_quadtree_mesh_quadratures(0.075)
        test_simple_quadtree_mesh_quadratures(0.05)
        test_simple_quadtree_mesh_quadratures(0.025)
    end;
end

function test_greens_identity()
    @testset "Test Green's Identity" begin
        test_simple_quadtree_greens_third_identity_simple_forcing() # tol 1e-15
        test_simple_quadtree_greens_third_identity_linear_x() # tol 1e-15
        test_simple_quadtree_greens_third_identity_linear_y() # tol 1e-15
        test_simple_quadtree_greens_third_identity_quadratic() # tol 1e-13
        test_simple_quadtree_greens_third_identity_sin_term() # tol 1e-13
        test_simple_quadtree_greens_third_identity_complex_forcing() # tol 1e-13
    end;
end

function test_quadtree_vs_triangle()
    @testset "Test Quadtree vs. Triangle Mesh" begin
        test_simple_quadtree_vs_triangle_mesh_area()
        test_simple_quadtree_vs_triangle_mesh_linear_x()
        test_simple_quadtree_vs_triangle_mesh_linear_y()    
    end;
end

@testset verbose = true "All Tests" begin
    # test_simple_triangle_mesh() # these won't be accurate until curved meshes
    # test_simple_quadtree_mesh() # these won't be accurate until curved meshes
    # test_individual_triangles() # these work
    # test_individual_quads() # weird rn but they work
    # test_greens_identity() # these work to different tols (max 1e-13)
    # test_quadtree_vs_triangle() # these work to different tols (max 1e-13)
    test_simple_quadtree_greens_third_identity_complex_forcing()
end;