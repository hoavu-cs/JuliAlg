using Test
using Graphs
using Combinatorics
using Random
using JuliOpt

const subgraph_density = JuliOpt.density

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

        S, λ, d = densest_subgraph(g)
        @test d ≈ 0.5 atol=1e-6
        @test Set(S) == Set([1, 2])
    end

    @testset "triangle K3" begin
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 3)

        S, λ, d = densest_subgraph(g)
        @test Set(S) == Set([1, 2, 3])
        @test d ≈ 1.0 atol=1e-6
    end

    @testset "complete graph K4" begin
        g = complete_graph(4)

        S, λ, d = densest_subgraph(g)
        @test Set(S) == Set(1:4)
        @test d ≈ 1.5 atol=1e-6
    end

    @testset "complete graph K5" begin
        g = complete_graph(5)

        S, λ, d = densest_subgraph(g)
        @test Set(S) == Set(1:5)
        @test d ≈ 2.0 atol=1e-6
    end

    @testset "K4 with pendant vertex" begin
        # K4 on {1,2,3,4}, vertex 5 connected only to vertex 1
        # K4 density: 6/4 = 1.5, whole graph density: 7/5 = 1.4
        g = complete_graph(4)
        add_vertex!(g)
        add_edge!(g, 1, 5)

        S, λ, d = densest_subgraph(g)
        @test Set(S) == Set([1, 2, 3, 4])
        @test d ≈ 1.5 atol=1e-6
    end

    @testset "star graph" begin
        # star_graph(5): vertex 1 center, edges to 2,3,4,5
        # density = 4/5 = 0.8
        g = star_graph(5)

        S, λ, d = densest_subgraph(g)
        @test d ≈ 0.8 atol=1e-6
        @test Set(S) == Set(1:5)
    end

    @testset "cycle graph" begin
        g = cycle_graph(6)
        # 6 edges, 6 vertices → density = 1.0

        S, λ, d = densest_subgraph(g)
        @test d ≈ 1.0 atol=1e-6
        @test Set(S) == Set(1:6)
    end

    @testset "no edges" begin
        g = SimpleGraph(4)

        S, λ, d = densest_subgraph(g)
        @test d ≈ 0.0 atol=1e-6
    end

    @testset "K5 with two pendants" begin
        # K5 on {1..5} has density 2.0
        # Adding pendants 6,7 gives density 12/7 ≈ 1.71 for whole graph
        g = complete_graph(5)
        add_vertex!(g)
        add_vertex!(g)
        add_edge!(g, 1, 6)
        add_edge!(g, 2, 7)

        S, λ, d = densest_subgraph(g)
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

        S, λ, d = densest_subgraph(g)
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

        S, λ, d = densest_subgraph(g)
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

        S, λ, d = densest_subgraph(g)
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

        S, λ, d = densest_subgraph(g)

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

        S, λ, d = densest_subgraph(g)
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

        S, λ, d = densest_subgraph(g)
        @test d ≈ best_d atol=1e-6
    end

end