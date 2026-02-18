using Test
using Graphs
using LinearAlgebra
using JuliAlg

"""
Compute ground-truth PageRank by building the Google matrix and solving
the linear system (I - M)r = 0 with sum(r) = 1.
"""
function pagerank_ground_truth(
    G::AbstractGraph;
    α::Float64=0.85,
    weights::Union{Dict{Tuple{Int, Int}, Float64}, Nothing}=nothing
)
    n = nv(G)
    n == 0 && return Float64[]

    # Build column-stochastic transition matrix H
    H = zeros(n, n)
    for v in 1:n
        nbrs = outneighbors(G, v)
        if !isempty(nbrs)
            if weights !== nothing
                total_w = sum(get(weights, (v, u), 0.0) for u in nbrs)
                if total_w > 0
                    for u in nbrs
                        w = get(weights, (v, u), 0.0)
                        H[u, v] = w / total_w
                    end
                else
                    # All weights zero → treat as dangling
                    H[:, v] .= 1.0 / n
                end
            else
                for u in nbrs
                    H[u, v] = 1.0 / length(nbrs)
                end
            end
        else
            # Dangling node: redistribute uniformly
            H[:, v] .= 1.0 / n
        end
    end

    # Google matrix: M = α * H + (1 - α) * (1/n) * 11ᵀ
    e = ones(n)
    M = α * H + (1.0 - α) / n * (e * e')

    # Stationary distribution: eigenvector of M for eigenvalue 1
    # Solve (M - I)r = 0 with constraint sum(r) = 1
    # Replace last row with sum constraint
    A = M - I
    A[end, :] .= 1.0
    b = zeros(n)
    b[end] = 1.0
    r = A \ b
    return r
end

@testset "PageRank Tests" begin

    @testset "Empty graph" begin
        g = SimpleDiGraph(0)
        r = JuliAlg.pagerank(g)
        @test isempty(r)
    end

    @testset "Single node" begin
        g = SimpleDiGraph(1)
        r = JuliAlg.pagerank(g)
        expected = pagerank_ground_truth(g)
        @test r ≈ expected atol=1e-6
    end

    @testset "Two-node directed edge" begin
        g = SimpleDiGraph(2)
        add_edge!(g, 1, 2)
        r = JuliAlg.pagerank(g, nothing, 0.85, 1000, 1e-10)
        expected = pagerank_ground_truth(g)
        @test r ≈ expected atol=1e-6
    end

    @testset "Complete graph K3" begin
        g = SimpleDiGraph(3)
        for i in 1:3, j in 1:3
            i != j && add_edge!(g, i, j)
        end
        r = JuliAlg.pagerank(g, nothing, 0.85, 1000, 1e-10)
        expected = pagerank_ground_truth(g)
        @test r ≈ expected atol=1e-6
    end

    @testset "Cycle graph C3" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 1)
        r = JuliAlg.pagerank(g, nothing, 0.85, 1000, 1e-10)
        expected = pagerank_ground_truth(g)
        @test r ≈ expected atol=1e-6
    end

    @testset "Star graph directed outward" begin
        g = SimpleDiGraph(5)
        for i in 2:5
            add_edge!(g, 1, i)
        end
        r = JuliAlg.pagerank(g, nothing, 0.85, 1000, 1e-10)
        expected = pagerank_ground_truth(g)
        @test r ≈ expected atol=1e-6
    end

    @testset "Star graph directed inward" begin
        g = SimpleDiGraph(5)
        for i in 2:5
            add_edge!(g, i, 1)
        end
        r = JuliAlg.pagerank(g, nothing, 0.85, 1000, 1e-10)
        expected = pagerank_ground_truth(g)
        @test r ≈ expected atol=1e-6
        # Center should have highest rank
        @test argmax(r) == 1
    end

    @testset "Dangling node handling" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        # Node 2 points nowhere, node 3 is isolated
        r = JuliAlg.pagerank(g, nothing, 0.85, 1000, 1e-10)
        expected = pagerank_ground_truth(g)
        @test r ≈ expected atol=1e-6
    end

    @testset "All nodes dangling (no edges)" begin
        g = SimpleDiGraph(5)
        r = JuliAlg.pagerank(g)
        expected = pagerank_ground_truth(g)
        @test r ≈ expected atol=1e-6
        @test all(isapprox.(r, 0.2; atol=1e-6))
    end

    @testset "Single self-loop" begin
        g = SimpleDiGraph(1)
        add_edge!(g, 1, 1)
        r = JuliAlg.pagerank(g)
        expected = pagerank_ground_truth(g)
        @test r ≈ expected atol=1e-6
    end

    @testset "Damping factor 0 gives uniform" begin
        g = SimpleDiGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        r = JuliAlg.pagerank(g, nothing, 0.0)
        expected = pagerank_ground_truth(g, α=0.0)
        @test r ≈ expected atol=1e-6
        @test all(isapprox.(r, 0.25; atol=1e-6))
    end

    @testset "Different α factors" begin
        g = SimpleDiGraph(5)
        for i in 2:5
            add_edge!(g, 1, i)
        end
        add_edge!(g, 3, 1)

        for α in [0.5, 0.85, 0.99]
            r = JuliAlg.pagerank(g, nothing, α, 1000, 1e-10)
            expected = pagerank_ground_truth(g, α=α)
            @test r ≈ expected atol=1e-5
        end
    end

    @testset "Weighted PageRank" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        weights = Dict((1, 2) => 2.0, (1, 3) => 1.0)
        r = JuliAlg.pagerank(g, weights, 0.85, 1000, 1e-10)
        expected = pagerank_ground_truth(g, weights=weights)
        @test r ≈ expected atol=1e-6
        @test r[2] > r[3]
    end

    @testset "Equal weights match unweighted" begin
        g = SimpleDiGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        add_edge!(g, 1, 4)
        weights = Dict((1, 2) => 1.0, (1, 3) => 1.0, (1, 4) => 1.0)
        r_weighted = JuliAlg.pagerank(g, weights)
        r_unweighted = JuliAlg.pagerank(g)
        @test r_weighted ≈ r_unweighted atol=1e-6
    end

    @testset "Missing weight defaults to zero" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        weights = Dict((1, 2) => 1.0)  # missing (1,3)
        r = JuliAlg.pagerank(g, weights, 0.85, 1000, 1e-10)
        expected = pagerank_ground_truth(g, weights=weights)
        @test r ≈ expected atol=1e-6
    end

    @testset "Zero weight edges" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        weights = Dict((1, 2) => 0.0, (1, 3) => 1.0)
        r = JuliAlg.pagerank(g, weights, 0.85, 1000, 1e-10)
        expected = pagerank_ground_truth(g, weights=weights)
        @test r ≈ expected atol=1e-6
    end

    @testset "Undirected graph" begin
        g = SimpleGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 1, 4)
        r = JuliAlg.pagerank(g, nothing, 0.85, 1000, 1e-10)
        expected = pagerank_ground_truth(g)
        @test r ≈ expected atol=1e-6
    end


    @testset "Larger random graph" begin
        import Random
        Random.seed!(42)
        g = SimpleDiGraph(50)
        for i in 1:50, j in 1:50
            if i != j && rand() < 0.08
                add_edge!(g, i, j)
            end
        end
        r = JuliAlg.pagerank(g, nothing, 0.85, 1000, 1e-10)
        expected = pagerank_ground_truth(g)
        @test r ≈ expected atol=1e-5
    end

    @testset "Random graph with high-indegree hub node" begin
        import Random
        Random.seed!(77)
        g = SimpleDiGraph(200)
        # Erdos-Renyi edges with p=0.01
        for i in 1:200, j in 1:200
            if i != j && rand() < 0.01
                add_edge!(g, i, j)
            end
        end
        # 100 random nodes point to node 1
        extra_sources = Random.shuffle(2:200)[1:100]
        for u in extra_sources
            add_edge!(g, u, 1)
        end
        r = JuliAlg.pagerank(g, nothing, 0.85, 1000, 1e-10)
        expected = pagerank_ground_truth(g)
        @test r ≈ expected atol=1e-5
        # Node 1 should have the highest rank due to extra incoming edges
        @test argmax(r) == 1
    end

    @testset "Larger weighted random graph" begin
        import Random
        Random.seed!(99)
        g = SimpleDiGraph(30)
        weights = Dict{Tuple{Int, Int}, Float64}()
        for i in 1:30, j in 1:30
            if i != j && rand() < 0.1
                add_edge!(g, i, j)
                weights[(i, j)] = rand() * 10.0
            end
        end
        r = JuliAlg.pagerank(g, weights, 0.85, 1000, 1e-10)
        expected = pagerank_ground_truth(g, weights=weights)
        @test r ≈ expected atol=1e-5
    end

    @testset "Large sparse graph (2000 nodes, 1000 edges)" begin
        import Random
        Random.seed!(42)
        g = SimpleDiGraph(10000)
        edges_added = 0
        while edges_added < 100000
            u, v = rand(1:10000), rand(1:10000)
            if u != v && !has_edge(g, u, v)
                add_edge!(g, u, v)
                edges_added += 1
            end
        end

        r = JuliAlg.pagerank(g, nothing, 0.85, 200, 1e-8)
        @test length(r) == 10000
        @test sum(r) ≈ 1.0 atol=1e-4
        @test all(r .>= 0.0)
    end
end