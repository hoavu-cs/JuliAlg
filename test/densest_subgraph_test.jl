using Test
using Graphs
using Combinatorics
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