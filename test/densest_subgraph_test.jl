using Test
using Graphs
using Combinatorics
using Random
using JuliAlg

const subgraph_density = JuliAlg.density

@testset "Densest Subgraph Tests" begin

    @testset "density helper" begin
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 3)

        @test subgraph_density(g, Int[]) == 0.0
        @test subgraph_density(g, [1]) == 0.0
        @test subgraph_density(g, [1, 2]) == 0.5
        @test subgraph_density(g, [1, 2, 3]) == 1.0
    end

    @testset "single edge" begin
        g = SimpleGraph(2)
        add_edge!(g, 1, 2)

        S, d = densest_subgraph(g)
        @test d ≈ 0.5 atol=1e-6
        @test Set(S) == Set([1, 2])
    end

    @testset "triangle K3" begin
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 3)

        S, d = densest_subgraph(g)
        @test Set(S) == Set([1, 2, 3])
        @test d ≈ 1.0 atol=1e-6
    end

    @testset "complete graph K4" begin
        g = complete_graph(4)

        S, d = densest_subgraph(g)
        @test Set(S) == Set(1:4)
        @test d ≈ 1.5 atol=1e-6
    end

    @testset "complete graph K5" begin
        g = complete_graph(5)

        S, d = densest_subgraph(g)
        @test Set(S) == Set(1:5)
        @test d ≈ 2.0 atol=1e-6
    end

    @testset "K4 with pendant vertex" begin
        # K4 on {1,2,3,4}, vertex 5 connected only to vertex 1
        # K4 density: 6/4 = 1.5, whole graph density: 7/5 = 1.4
        g = complete_graph(4)
        add_vertex!(g)
        add_edge!(g, 1, 5)

        S, d = densest_subgraph(g)
        @test Set(S) == Set([1, 2, 3, 4])
        @test d ≈ 1.5 atol=1e-6
    end

    @testset "star graph" begin
        # star_graph(5): vertex 1 center, edges to 2,3,4,5
        # density = 4/5 = 0.8
        g = star_graph(5)

        S, d = densest_subgraph(g)
        @test d ≈ 0.8 atol=1e-6
        @test Set(S) == Set(1:5)
    end

    @testset "cycle graph" begin
        g = cycle_graph(6)
        # 6 edges, 6 vertices → density = 1.0

        S, d = densest_subgraph(g)
        @test d ≈ 1.0 atol=1e-6
        @test Set(S) == Set(1:6)
    end

    @testset "no edges" begin
        g = SimpleGraph(4)

        S, d = densest_subgraph(g)
        @test d ≈ 0.0 atol=1e-6
    end

    @testset "directed graph throws error" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        @test_throws ArgumentError densest_subgraph(g)
    end

    @testset "K5 with two pendants" begin
        # K5 on {1..5} has density 2.0
        # Adding pendants 6,7 gives density 12/7 ≈ 1.71 for whole graph
        g = complete_graph(5)
        add_vertex!(g)
        add_vertex!(g)
        add_edge!(g, 1, 6)
        add_edge!(g, 2, 7)

        S, d = densest_subgraph(g)
        @test Set(S) == Set(1:5)
        @test d ≈ 2.0 atol=1e-6
    end

    @testset "dense core with sparse periphery" begin
        # K5 on {1..5} (density 2.0) + 5 pendant vertices
        g = complete_graph(5)
        for i in 6:10
            add_vertex!(g)
            add_edge!(g, i - 5, i)
        end

        S, d = densest_subgraph(g)
        @test Set(S) == Set(1:5)
        @test d ≈ 2.0 atol=1e-6
    end

    @testset "K5 and K3 connected by a bridge edge" begin
        # K5 on {1..5} (density 2.0), K3 on {6,7,8} (density 1.0)
        # connected by edge (5,6)
        # Whole graph: 14 edges / 8 vertices = 1.75
        # Densest subgraph should be K5
        g = complete_graph(5)
        for v in 6:8
            add_vertex!(g)
        end
        add_edge!(g, 6, 7)
        add_edge!(g, 7, 8)
        add_edge!(g, 6, 8)
        add_edge!(g, 5, 6)  # bridge

        S, d = densest_subgraph(g)
        @test Set(S) == Set(1:5)
        @test d ≈ 2.0 atol=1e-6
    end

    @testset "K_{3,4} -- K3 -- square chain" begin
        # K_{3,4} on {1..7}: parts {1,2,3} and {4,5,6,7}, 12 edges, density 12/7 ≈ 1.714
        # K3 on {8,9,10}: 3 edges, density 1.0
        # Square on {11,12,13,14}: 4 edges, density 1.0
        # Bridges: (7,8) and (10,11)
        # Whole graph: 21 edges / 14 vertices = 1.5
        # Densest subgraph should be K_{3,4}
        g = SimpleGraph(14)
        for a in 1:3, b in 4:7
            add_edge!(g, a, b)
        end
        add_edge!(g, 8, 9)
        add_edge!(g, 9, 10)
        add_edge!(g, 8, 10)
        add_edge!(g, 11, 12)
        add_edge!(g, 12, 13)
        add_edge!(g, 13, 14)
        add_edge!(g, 14, 11)
        add_edge!(g, 7, 8)   # bridge K_{3,4} -- K3
        add_edge!(g, 10, 11)  # bridge K3 -- square

        # Brute force verification
        best_d = 0.0
        for k in 1:nv(g)
            for subset in combinations(1:nv(g), k)
                d = subgraph_density(g, subset)
                best_d = max(best_d, d)
            end
        end

        S, d = densest_subgraph(g)
        @test Set(S) == Set(1:7)
        @test d ≈ 12.0 / 7.0 atol=1e-6
        @test d ≈ best_d atol=1e-6
    end

    @testset "random G(1000, 0.01) with planted K10" begin
        Random.seed!(42)
        n = 1000
        p = 0.01
        clique_vertices = 1:10

        g = SimpleGraph(n)
        for i in 1:n, j in (i+1):n
            if rand() < p
                add_edge!(g, i, j)
            end
        end

        # Plant K10 on vertices 1..10
        for i in clique_vertices, j in clique_vertices
            i < j && add_edge!(g, i, j)
        end

        S, d = densest_subgraph(g)

        # Returned density must match the density helper
        @test d ≈ subgraph_density(g, S) atol=1e-6
        # Must be at least as dense as the planted K10
        @test d ≥ subgraph_density(g, collect(clique_vertices)) - 1e-6
        # The planted clique should be fully contained in the densest subgraph
        @test issubset(Set(clique_vertices), Set(S))
    end

    @testset "brute force verification: small graph" begin
        # Triangle {1,2,3} + square {3,4,5,6} sharing vertex 3
        g = SimpleGraph(6)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 5)
        add_edge!(g, 5, 6)
        add_edge!(g, 4, 6)

        # Brute force over all non-empty subsets
        best_d = 0.0
        for k in 1:nv(g)
            for subset in combinations(1:nv(g), k)
                d = subgraph_density(g, subset)
                best_d = max(best_d, d)
            end
        end

        S, d = densest_subgraph(g)
        @test d ≈ best_d atol=1e-6
    end

    @testset "brute force verification: irregular graph" begin
        # Petersen-like graph with extra edges
        g = SimpleGraph(8)
        for (u, v) in [(1,2),(1,3),(1,5),(2,3),(2,4),(3,6),(4,5),(4,7),(5,8),(6,7),(6,8),(7,8)]
            add_edge!(g, u, v)
        end

        best_d = 0.0
        for k in 1:nv(g)
            for subset in combinations(1:nv(g), k)
                d = subgraph_density(g, subset)
                best_d = max(best_d, d)
            end
        end

        S, d = densest_subgraph(g)
        @test d ≈ best_d atol=1e-6
    end

end

@testset "Densest Subgraph Peeling Tests" begin
    densest_peeling = JuliAlg.densest_subgraph_peeling

    @testset "single edge" begin
        g = SimpleGraph(2)
        add_edge!(g, 1, 2)

        S, d = densest_peeling(g)
        @test d ≈ 0.5 atol=1e-6
        @test Set(S) == Set([1, 2])
    end

    @testset "triangle K3" begin
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 3)

        S, d = densest_peeling(g)
        optimal = 1.0
        @test d ≥ optimal / 2 - 1e-6
        @test d ≈ subgraph_density(g, collect(S)) atol=1e-6
    end

    @testset "complete graph K4" begin
        g = complete_graph(4)

        S, d = densest_peeling(g)
        optimal = 1.5
        @test d ≥ optimal / 2 - 1e-6
        @test d ≈ subgraph_density(g, collect(S)) atol=1e-6
    end

    @testset "complete graph K5" begin
        g = complete_graph(5)

        S, d = densest_peeling(g)
        optimal = 2.0
        @test d ≥ optimal / 2 - 1e-6
        @test d ≈ subgraph_density(g, collect(S)) atol=1e-6
    end

    @testset "K4 with pendant vertex" begin
        g = complete_graph(4)
        add_vertex!(g)
        add_edge!(g, 1, 5)

        S, d = densest_peeling(g)
        optimal = 1.5
        @test d ≥ optimal / 2 - 1e-6
        @test d ≈ subgraph_density(g, collect(S)) atol=1e-6
    end

    @testset "star graph" begin
        g = star_graph(5)

        S, d = densest_peeling(g)
        optimal = 0.8
        @test d ≥ optimal / 2 - 1e-6
        @test d ≈ subgraph_density(g, collect(S)) atol=1e-6
    end

    @testset "cycle graph" begin
        g = cycle_graph(6)

        S, d = densest_peeling(g)
        optimal = 1.0
        @test d ≥ optimal / 2 - 1e-6
        @test d ≈ subgraph_density(g, collect(S)) atol=1e-6
    end

    @testset "no edges" begin
        g = SimpleGraph(4)

        S, d = densest_peeling(g)
        @test d ≈ 0.0 atol=1e-6
    end

    @testset "directed graph throws error" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        @test_throws ArgumentError densest_peeling(g)
    end

    @testset "K5 with two pendants" begin
        g = complete_graph(5)
        add_vertex!(g)
        add_vertex!(g)
        add_edge!(g, 1, 6)
        add_edge!(g, 2, 7)

        S, d = densest_peeling(g)
        optimal = 2.0
        @test d ≥ optimal / 2 - 1e-6
        @test d ≈ subgraph_density(g, collect(S)) atol=1e-6
    end

    @testset "dense core with sparse periphery" begin
        g = complete_graph(5)
        for i in 6:10
            add_vertex!(g)
            add_edge!(g, i - 5, i)
        end

        S, d = densest_peeling(g)
        optimal = 2.0
        @test d ≥ optimal / 2 - 1e-6
        @test d ≈ subgraph_density(g, collect(S)) atol=1e-6
    end

    @testset "K5 and K3 connected by a bridge edge" begin
        g = complete_graph(5)
        for v in 6:8
            add_vertex!(g)
        end
        add_edge!(g, 6, 7)
        add_edge!(g, 7, 8)
        add_edge!(g, 6, 8)
        add_edge!(g, 5, 6)

        S, d = densest_peeling(g)
        optimal = 2.0
        @test d ≥ optimal / 2 - 1e-6
        @test d ≈ subgraph_density(g, collect(S)) atol=1e-6
    end

    @testset "K_{3,4} -- K3 -- square chain" begin
        g = SimpleGraph(14)
        for a in 1:3, b in 4:7
            add_edge!(g, a, b)
        end
        add_edge!(g, 8, 9)
        add_edge!(g, 9, 10)
        add_edge!(g, 8, 10)
        add_edge!(g, 11, 12)
        add_edge!(g, 12, 13)
        add_edge!(g, 13, 14)
        add_edge!(g, 14, 11)
        add_edge!(g, 7, 8)
        add_edge!(g, 10, 11)

        # Brute force optimal
        best_d = 0.0
        for k in 1:nv(g)
            for subset in combinations(1:nv(g), k)
                d = subgraph_density(g, subset)
                best_d = max(best_d, d)
            end
        end

        S, d = densest_peeling(g)
        @test d ≥ best_d / 2 - 1e-6
        @test d ≈ subgraph_density(g, collect(S)) atol=1e-6
    end

    @testset "random G(1000, 0.01) with planted K10" begin
        Random.seed!(42)
        n = 1000
        p = 0.01
        clique_vertices = 1:10

        g = SimpleGraph(n)
        for i in 1:n, j in (i+1):n
            if rand() < p
                add_edge!(g, i, j)
            end
        end

        for i in clique_vertices, j in clique_vertices
            i < j && add_edge!(g, i, j)
        end

        _, optimal = densest_subgraph(g)

        S, d = densest_peeling(g)
        @test d ≥ optimal / 2 - 1e-6
        @test d ≈ subgraph_density(g, collect(S)) atol=1e-6
    end

    @testset "brute force verification: small graph" begin
        g = SimpleGraph(6)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 5)
        add_edge!(g, 5, 6)
        add_edge!(g, 4, 6)

        best_d = 0.0
        for k in 1:nv(g)
            for subset in combinations(1:nv(g), k)
                d = subgraph_density(g, subset)
                best_d = max(best_d, d)
            end
        end

        S, d = densest_peeling(g)
        @test d ≥ best_d / 2 - 1e-6
        @test d ≈ subgraph_density(g, collect(S)) atol=1e-6
    end

    @testset "brute force verification: irregular graph" begin
        g = SimpleGraph(8)
        for (u, v) in [(1,2),(1,3),(1,5),(2,3),(2,4),(3,6),(4,5),(4,7),(5,8),(6,7),(6,8),(7,8)]
            add_edge!(g, u, v)
        end

        best_d = 0.0
        for k in 1:nv(g)
            for subset in combinations(1:nv(g), k)
                d = subgraph_density(g, subset)
                best_d = max(best_d, d)
            end
        end

        S, d = densest_peeling(g)
        @test d ≥ best_d / 2 - 1e-6
        @test d ≈ subgraph_density(g, collect(S)) atol=1e-6
    end
end

@testset "Densest At-Most-K Subgraph Tests" begin
    densest_at_most_k = JuliAlg.densest_at_most_k_subgraph

    # Brute force reference implementation
    function bf_densest_k(G, k)
        best_d = 0.0
        best_S = Int[]
        for size in 1:min(k, nv(G))
            for subset in combinations(1:nv(G), size)
                d = subgraph_density(G, subset)
                if d > best_d
                    best_d = d
                    best_S = subset
                end
            end
        end
        return best_S, best_d
    end

    @testset "k >= n delegates to densest_subgraph" begin
        g = complete_graph(4)
        S, d = densest_at_most_k(g, 5)
        @test Set(S) == Set(1:4)
        @test d ≈ 1.5 atol=1e-6
    end

    @testset "k == n delegates to densest_subgraph" begin
        g = complete_graph(3)
        S, d = densest_at_most_k(g, 3)
        @test Set(S) == Set(1:3)
        @test d ≈ 1.0 atol=1e-6
    end

    @testset "directed graph throws error" begin
        g = SimpleDiGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        @test_throws ArgumentError densest_at_most_k(g, 2)
    end

    @testset "triangle K3, k=2" begin
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 3)

        S, d = densest_at_most_k(g, 2)
        @test length(S) ≤ 2
        @test d ≈ 0.5 atol=1e-6
    end

    @testset "complete graph K5, k=3" begin
        g = complete_graph(5)
        S, d = densest_at_most_k(g, 3)
        @test length(S) ≤ 3
        @test d ≈ 1.0 atol=1e-6  # any 3 vertices form K3: 3 edges / 3 = 1.0
    end

    @testset "complete graph K5, k=4" begin
        g = complete_graph(5)
        S, d = densest_at_most_k(g, 4)
        @test length(S) ≤ 4
        @test d ≈ 1.5 atol=1e-6  # any 4 vertices form K4: 6 edges / 4 = 1.5
    end

    @testset "K4 with pendant, k=4" begin
        # K4 on {1,2,3,4}, vertex 5 connected only to 1
        # Best k=4 is K4 with density 6/4 = 1.5
        g = complete_graph(4)
        add_vertex!(g)
        add_edge!(g, 1, 5)

        S, d = densest_at_most_k(g, 4)
        @test length(S) ≤ 4
        @test Set(S) == Set(1:4)
        @test d ≈ 1.5 atol=1e-6
    end

    @testset "K4 with pendant, k=3" begin
        g = complete_graph(4)
        add_vertex!(g)
        add_edge!(g, 1, 5)

        S, d = densest_at_most_k(g, 3)
        @test length(S) ≤ 3
        @test d ≈ 1.0 atol=1e-6  # any 3 from K4 form K3
    end

    @testset "star graph, k=2" begin
        # star_graph(5): vertex 1 center, edges to 2,3,4,5
        # Best k=2: center + any leaf = 1 edge / 2 = 0.5
        g = star_graph(5)
        S, d = densest_at_most_k(g, 2)
        @test length(S) ≤ 2
        @test d ≈ 0.5 atol=1e-6
        @test 1 in S  # center must be included
    end

    @testset "no edges, k=2" begin
        g = SimpleGraph(4)
        S, d = densest_at_most_k(g, 2)
        @test d ≈ 0.0 atol=1e-6
    end

    @testset "single edge, k=1" begin
        # Single vertex has no edges, density = 0
        g = SimpleGraph(2)
        add_edge!(g, 1, 2)
        S, d = densest_at_most_k(g, 1)
        @test d ≈ 0.0 atol=1e-6
    end

    @testset "K5+K3 bridge, k=5" begin
        g = complete_graph(5)
        for v in 6:8
            add_vertex!(g)
        end
        add_edge!(g, 6, 7)
        add_edge!(g, 7, 8)
        add_edge!(g, 6, 8)
        add_edge!(g, 5, 6)  # bridge

        S, d = densest_at_most_k(g, 5)
        @test length(S) ≤ 5
        @test Set(S) == Set(1:5)
        @test d ≈ 2.0 atol=1e-6
    end

    @testset "K5+K3 bridge, k=3" begin
        g = complete_graph(5)
        for v in 6:8
            add_vertex!(g)
        end
        add_edge!(g, 6, 7)
        add_edge!(g, 7, 8)
        add_edge!(g, 6, 8)
        add_edge!(g, 5, 6)

        S, d = densest_at_most_k(g, 3)
        @test length(S) ≤ 3
        @test d ≈ 1.0 atol=1e-6  # any K3 gives density 1.0
    end

    @testset "cycle graph C6, k=4" begin
        g = cycle_graph(6)
        S, d = densest_at_most_k(g, 4)
        _, bf_d = bf_densest_k(g, 4)
        @test length(S) ≤ 4
        @test d ≈ bf_d atol=1e-6
    end

    @testset "brute force: triangle + square sharing vertex" begin
        # Triangle {1,2,3} + square {3,4,5,6} sharing vertex 3
        g = SimpleGraph(6)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 5)
        add_edge!(g, 5, 6)
        add_edge!(g, 4, 6)

        for k in 1:5
            S, d = densest_at_most_k(g, k)
            _, bf_d = bf_densest_k(g, k)
            @test length(S) ≤ k
            @test d ≈ bf_d atol=1e-6
        end
    end

    @testset "brute force: irregular graph" begin
        g = SimpleGraph(8)
        for (u, v) in [(1,2),(1,3),(1,5),(2,3),(2,4),(3,6),(4,5),(4,7),(5,8),(6,7),(6,8),(7,8)]
            add_edge!(g, u, v)
        end

        for k in 1:7
            S, d = densest_at_most_k(g, k)
            _, bf_d = bf_densest_k(g, k)
            @test length(S) ≤ k
            @test d ≈ bf_d atol=1e-6
        end
    end

    @testset "brute force: K_{3,4} chain graph" begin
        g = SimpleGraph(14)
        for a in 1:3, b in 4:7
            add_edge!(g, a, b)
        end
        add_edge!(g, 8, 9)
        add_edge!(g, 9, 10)
        add_edge!(g, 8, 10)
        add_edge!(g, 11, 12)
        add_edge!(g, 12, 13)
        add_edge!(g, 13, 14)
        add_edge!(g, 14, 11)
        add_edge!(g, 7, 8)
        add_edge!(g, 10, 11)

        for k in [3, 5, 7]
            S, d = densest_at_most_k(g, k)
            _, bf_d = bf_densest_k(g, k)
            @test length(S) ≤ k
            @test d ≈ bf_d atol=1e-6
        end
    end

    @testset "brute force: random G(12, 0.3)" begin
        Random.seed!(123)
        g = SimpleGraph(12)
        for i in 1:12, j in (i+1):12
            if rand() < 0.3
                add_edge!(g, i, j)
            end
        end

        for k in [2, 4, 6, 8]
            S, d = densest_at_most_k(g, k)
            _, bf_d = bf_densest_k(g, k)
            @test length(S) ≤ k
            @test d ≈ bf_d atol=1e-6
        end
    end

    @testset "planted K7 in random G(50, 0.1), k=5" begin
        Random.seed!(99)
        n = 50
        clique_vertices = 1:7

        g = SimpleGraph(n)
        for i in 1:n, j in (i+1):n
            if rand() < 0.1
                add_edge!(g, i, j)
            end
        end

        # Plant K7 on vertices 1..7
        for i in clique_vertices, j in clique_vertices
            i < j && add_edge!(g, i, j)
        end

        S, d = densest_at_most_k(g, 5)
        @test length(S) ≤ 5

        # Any 5 vertices from K7 form K5 with density C(5,2)/5 = 2.0
        # The result should be at least as dense as that
        k5_density = subgraph_density(g, collect(clique_vertices)[1:5])
        @test d ≥ k5_density - 1e-6

        # The selected vertices should come from the planted clique
        @test issubset(Set(S), Set(clique_vertices))
    end

    @testset "damaged K6 and K5 planted in random G(100, 0.1), k=5" begin
        Random.seed!(42)
        n = 100
        clique1 = collect(1:6)   # K6 with 5 edges removed
        clique2 = collect(7:11)  # K5 with 4 edges removed

        g = SimpleGraph(n)
        for i in 1:n, j in (i+1):n
            if rand() < 0.1
                add_edge!(g, i, j)
            end
        end

        # Plant K6 on vertices 1..6, then remove 5 random edges
        k6_edges = [(i, j) for i in clique1 for j in clique1 if i < j]
        for (i, j) in k6_edges
            add_edge!(g, i, j)
        end
        Random.seed!(7)
        for (i, j) in shuffle(k6_edges)[1:5]
            rem_edge!(g, i, j)
        end

        # Plant K5 on vertices 7..11, then remove 4 random edges
        k5_edges = [(i, j) for i in clique2 for j in clique2 if i < j]
        for (i, j) in k5_edges
            add_edge!(g, i, j)
        end
        Random.seed!(13)
        for (i, j) in shuffle(k5_edges)[1:4]
            rem_edge!(g, i, j)
        end

        S, d = densest_at_most_k(g, 5)
        @test length(S) ≤ 5

        # Brute force best 5-subset from the planted vertices (C(11,5) = 462)
        planted = vcat(clique1, clique2)
        best_planted_d = 0.0
        for subset in combinations(planted, 5)
            sd = subgraph_density(g, subset)
            best_planted_d = max(best_planted_d, sd)
        end

        # Algorithm should find density at least as good as the best planted subset
        @test d ≥ best_planted_d - 1e-6
    end
end