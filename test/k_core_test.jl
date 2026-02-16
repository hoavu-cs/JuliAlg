using Test
using Graphs
using JuliAlg

@testset "K-Core Decomposition Tests" begin

    @testset "empty graph" begin
        g = SimpleGraph(0)
        core = k_core_decomposition(g)
        @test isempty(core)
    end

    @testset "single vertex, no edges" begin
        g = SimpleGraph(1)
        core = k_core_decomposition(g)
        @test core[1] == 0
    end

    @testset "isolated vertices" begin
        g = SimpleGraph(4)
        core = k_core_decomposition(g)
        for v in 1:4
            @test core[v] == 0
        end
    end

    @testset "single edge" begin
        g = SimpleGraph(2)
        add_edge!(g, 1, 2)
        core = k_core_decomposition(g)
        @test core[1] == 1
        @test core[2] == 1
    end

    @testset "triangle K3" begin
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 3)
        core = k_core_decomposition(g)
        for v in 1:3
            @test core[v] == 2
        end
    end

    @testset "complete graph K4" begin
        g = complete_graph(4)
        core = k_core_decomposition(g)
        for v in 1:4
            @test core[v] == 3
        end
    end

    @testset "complete graph K5" begin  
        g = complete_graph(5)       
        core = k_core_decomposition(g)
        for v in 1:5
            @test core[v] == 4
        end
    end

    @testset "path graph" begin
        # 1 - 2 - 3 - 4 - 5
        # All vertices have core number 1 (endpoints have degree 1)
        g = path_graph(5)
        core = k_core_decomposition(g)
        for v in 1:5
            @test core[v] == 1
        end
    end

    @testset "star graph" begin
        # Center vertex 1 connected to 2,3,4,5
        # All have core number 1 (leaves have degree 1)
        g = star_graph(5)
        core = k_core_decomposition(g)
        for v in 1:5
            @test core[v] == 1
        end
    end

    @testset "cycle graph C6" begin
        g = cycle_graph(6)
        core = k_core_decomposition(g)
        for v in 1:6
            @test core[v] == 2
        end
    end

    @testset "K4 with pendant vertex" begin
        # K4 on {1,2,3,4}, vertex 5 connected only to vertex 1
        # K4 vertices: core 3, pendant: core 1
        g = complete_graph(4)
        add_vertex!(g)
        add_edge!(g, 1, 5)
        core = k_core_decomposition(g)
        for v in 1:4
            @test core[v] == 3
        end
        @test core[5] == 1
    end

    @testset "K4 + K3 connected by bridge" begin
        # K4 on {1,2,3,4} (core 3), K3 on {5,6,7} (core 2), bridge (4,5)
        g = complete_graph(4)
        for v in 5:7
            add_vertex!(g)
        end
        add_edge!(g, 5, 6)
        add_edge!(g, 6, 7)
        add_edge!(g, 5, 7)
        add_edge!(g, 4, 5)
        core = k_core_decomposition(g)
        for v in 1:4
            @test core[v] == 3
        end
        for v in 5:7
            @test core[v] == 2
        end
    end

    @testset "tree graph" begin
        # Trees have all core numbers = 1 (except isolated), since removing
        # leaves iteratively peels the whole tree at k=1
        g = SimpleGraph(7)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        add_edge!(g, 2, 4)
        add_edge!(g, 2, 5)
        add_edge!(g, 3, 6)
        add_edge!(g, 3, 7)
        core = k_core_decomposition(g)
        for v in 1:7
            @test core[v] == 1
        end
    end

    @testset "Petersen graph" begin
        # Petersen graph is 3-regular, and its k-core decomposition gives all vertices core 3
        g = smallgraph(:petersen)
        core = k_core_decomposition(g)
        for v in 1:10
            @test core[v] == 3
        end
    end

    @testset "two triangles sharing an edge" begin
        # Vertices 1,2,3 form triangle, vertices 2,3,4 form triangle
        # Shared edge (2,3). All vertices have core 2.
        g = SimpleGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        add_edge!(g, 2, 3)
        add_edge!(g, 2, 4)
        add_edge!(g, 3, 4)
        core = k_core_decomposition(g)
        for v in 1:4
            @test core[v] == 2
        end
    end

    @testset "core number ≤ degree for all vertices" begin
        # Property: core(v) ≤ degree(v) always holds
        g = SimpleGraph(8)
        for (u, v) in [(1,2),(1,3),(1,5),(2,3),(2,4),(3,6),(4,5),(4,7),(5,8),(6,7),(6,8),(7,8)]
            add_edge!(g, u, v)
        end
        core = k_core_decomposition(g)
        for v in vertices(g)
            @test core[v] ≤ degree(g, v)
        end
    end

    @testset "k-core subgraph property" begin
        # For each core number k, the subgraph induced by vertices with core ≥ k
        # should have minimum degree ≥ k
        g = complete_graph(5)
        for v in 6:8
            add_vertex!(g)
        end
        add_edge!(g, 5, 6)
        add_edge!(g, 6, 7)
        add_edge!(g, 7, 8)
        add_edge!(g, 6, 8)

        core = k_core_decomposition(g)
        max_core = maximum(values(core))

        for k in 0:max_core
            verts_k = [v for v in vertices(g) if core[v] ≥ k]
            if !isempty(verts_k)
                sg, vmap = induced_subgraph(g, verts_k)
                for v in vertices(sg)
                    @test degree(sg, v) ≥ k
                end
            end
        end
    end

    @testset "mixed components" begin
        # Isolated vertex (core 0) + edge (core 1) + triangle (core 2)
        g = SimpleGraph(6)
        # vertex 1: isolated
        add_edge!(g, 2, 3)           # edge
        add_edge!(g, 4, 5)           # triangle
        add_edge!(g, 5, 6)
        add_edge!(g, 4, 6)
        core = k_core_decomposition(g)
        @test core[1] == 0
        @test core[2] == 1
        @test core[3] == 1
        for v in 4:6
            @test core[v] == 2
        end
    end

    @testset "complex graph with multiple cores" begin
        g = SimpleGraph(15)
        add_edge!(g, 1, 3)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 3, 5)
        add_edge!(g, 3, 8)
        add_edge!(g, 6, 8)
        add_edge!(g, 6, 10)
        add_edge!(g, 6, 11)
        add_edge!(g, 6, 9)
        add_edge!(g, 7, 8)
        add_edge!(g, 7, 9)
        add_edge!(g, 8, 9)
        add_edge!(g, 8, 10)
        add_edge!(g, 9, 10)
        add_edge!(g, 9, 14)
        add_edge!(g, 10, 11)
        add_edge!(g, 10, 13)
        add_edge!(g, 10, 14)
        add_edge!(g, 11, 12)
        add_edge!(g, 13, 15)
        add_edge!(g, 14, 15)

        true_core = Dict(1 => 1, 
                        2 => 1, 
                        3 => 1,
                        4 => 1,
                        5 => 1,
                        12 => 1,
                        11 => 2,
                        7 => 2,
                        14 => 2,
                        15 => 2,
                        13 => 2,
                        8 => 3,
                        6 => 3,
                        9 => 3,
                        10 => 3)

        core = k_core_decomposition(g)
        for v in vertices(g)
            @test core[v] == true_core[v]
        end
    end

end