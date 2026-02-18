using Test
using Graphs
using Combinatorics
using JuliAlg

n_simulations = 10000

@testset "Influence Maximization Tests" begin

    @testset "small graph 1: star topology" begin
        g = SimpleDiGraph(5)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        add_edge!(g, 1, 4)
        add_edge!(g, 2, 1)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 2)
        
        weights = Dict(
            (1, 2) => 0.8,
            (1, 3) => 0.5,
            (1, 4) => 0.3,
            (2, 1) => 0.6,
            (2, 3) => 0.9,
            (3, 4) => 0.7,
            (4, 2) => 0.4
        )
        k = 2
        
        best_spread = 0.0
        for subset in combinations(1:5, k)
            spread = simulate_ic(g, weights, collect(subset); n_simulations=n_simulations)
            best_spread = max(best_spread, spread)
        end
        
        solution, spread = influence_maximization_ic(g, weights, k; n_simulations_small=n_simulations ÷ 10, n_simulations=n_simulations)
        @test spread >= 0.63 * best_spread
    end

    @testset "small graph 2: linear chain" begin
        g = SimpleDiGraph(6)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 5)
        add_edge!(g, 5, 6)
        
        weights = Dict(
            (1, 2) => 0.7,
            (2, 3) => 0.6,
            (3, 4) => 0.8,
            (4, 5) => 0.5,
            (5, 6) => 0.9
        )
        k = 2
        
        best_spread = 0.0
        for subset in combinations(1:6, k)
            spread = simulate_ic(g, weights, collect(subset); n_simulations=n_simulations)
            best_spread = max(best_spread, spread)
        end
        
        solution, spread = influence_maximization_ic(g, weights, k; n_simulations_small=n_simulations ÷ 10, n_simulations=n_simulations)
        @test spread >= 0.63 * best_spread
    end

    @testset "small graph 3: bidirectional triangle" begin
        g = SimpleDiGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 1)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 2)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 3)
        add_edge!(g, 1, 4)
        add_edge!(g, 4, 1)
        
        weights = Dict(
            (1, 2) => 0.9,
            (2, 1) => 0.8,
            (2, 3) => 0.7,
            (3, 2) => 0.6,
            (3, 4) => 0.8,
            (4, 3) => 0.7,
            (1, 4) => 0.5,
            (4, 1) => 0.4
        )
        k = 1
        
        best_spread = 0.0
        for subset in combinations(1:4, k)
            spread = simulate_ic(g, weights, collect(subset); n_simulations=n_simulations)
            best_spread = max(best_spread, spread)
        end
        
        solution, spread = influence_maximization_ic(g, weights, k; n_simulations_small=n_simulations ÷ 10, n_simulations=n_simulations)
        @test spread >= 0.63 * best_spread
    end

    @testset "small graph 4: complete directed graph" begin
        g = SimpleDiGraph(5)
        for i in 1:5
            for j in 1:5
                if i != j
                    add_edge!(g, i, j)
                end
            end
        end
        
        weights = Dict{Tuple{Int, Int}, Float64}()
        for i in 1:5
            for j in 1:5
                if i != j
                    weights[(i, j)] = 0.3 + 0.1 * ((i + j) % 5)
                end
            end
        end
        k = 2
        
        best_spread = 0.0
        for subset in combinations(1:5, k)
            spread = simulate_ic(g, weights, collect(subset); n_simulations=n_simulations)
            best_spread = max(best_spread, spread)
        end
        
        solution, spread = influence_maximization_ic(g, weights, k; n_simulations_small=n_simulations ÷ 10, n_simulations=n_simulations)
        @test spread >= 0.63 * best_spread
    end

    @testset "small graph 5: hub and spoke" begin
        g = SimpleDiGraph(7)
        # Central hub (node 1) connects to all spokes
        for i in 2:7
            add_edge!(g, 1, i)
        end
        # Some spokes connect back to hub
        add_edge!(g, 3, 1)
        add_edge!(g, 5, 1)
        # Some inter-spoke connections
        add_edge!(g, 2, 4)
        add_edge!(g, 4, 6)
        
        weights = Dict(
            (1, 2) => 0.7,
            (1, 3) => 0.8,
            (1, 4) => 0.6,
            (1, 5) => 0.9,
            (1, 6) => 0.5,
            (1, 7) => 0.7,
            (3, 1) => 0.4,
            (5, 1) => 0.3,
            (2, 4) => 0.6,
            (4, 6) => 0.5
        )
        k = 2
        
        best_spread = 0.0
        for subset in combinations(1:7, k)
            spread = simulate_ic(g, weights, collect(subset); n_simulations=n_simulations)
            best_spread = max(best_spread, spread)
        end
        
        solution, spread = influence_maximization_ic(g, weights, k; n_simulations_small=n_simulations ÷ 10, n_simulations=n_simulations)
        @test spread >= 0.63 * best_spread
    end

    @testset "small graph 6: cycle with shortcuts" begin
        g = SimpleDiGraph(6)
        # Main cycle
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 5)
        add_edge!(g, 5, 6)
        add_edge!(g, 6, 1)
        # Shortcuts
        add_edge!(g, 1, 4)
        add_edge!(g, 3, 6)
        
        weights = Dict(
            (1, 2) => 0.7,
            (2, 3) => 0.8,
            (3, 4) => 0.6,
            (4, 5) => 0.7,
            (5, 6) => 0.8,
            (6, 1) => 0.5,
            (1, 4) => 0.9,
            (3, 6) => 0.8
        )
        k = 3
        
        best_spread = 0.0
        for subset in combinations(1:6, k)
            spread = simulate_ic(g, weights, collect(subset); n_simulations=n_simulations)
            best_spread = max(best_spread, spread)
        end
        
        solution, spread = influence_maximization_ic(g, weights, k; n_simulations_small=n_simulations ÷ 10, n_simulations=n_simulations)
        @test spread >= 0.63 * best_spread
    end

    @testset "small graph 7: binary tree" begin
        g = SimpleDiGraph(7)
        # Level 0 -> Level 1
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        # Level 1 -> Level 2
        add_edge!(g, 2, 4)
        add_edge!(g, 2, 5)
        add_edge!(g, 3, 6)
        add_edge!(g, 3, 7)
        
        weights = Dict(
            (1, 2) => 0.8,
            (1, 3) => 0.8,
            (2, 4) => 0.7,
            (2, 5) => 0.7,
            (3, 6) => 0.6,
            (3, 7) => 0.6
        )
        k = 1
        
        best_spread = 0.0
        for subset in combinations(1:7, k)
            spread = simulate_ic(g, weights, collect(subset); n_simulations=n_simulations)
            best_spread = max(best_spread, spread)
        end
        
        solution, spread = influence_maximization_ic(g, weights, k; n_simulations_small=n_simulations ÷ 10, n_simulations=n_simulations)
        @test spread >= 0.63 * best_spread
    end

    @testset "small graph 8: dense core with periphery" begin
        g = SimpleDiGraph(8)
        # Dense core (nodes 1-4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 1)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 2)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 3)
        add_edge!(g, 4, 1)
        add_edge!(g, 1, 4)
        # Periphery connections
        add_edge!(g, 1, 5)
        add_edge!(g, 2, 6)
        add_edge!(g, 3, 7)
        add_edge!(g, 4, 8)
        
        weights = Dict(
            (1, 2) => 0.9,
            (2, 1) => 0.9,
            (2, 3) => 0.9,
            (3, 2) => 0.9,
            (3, 4) => 0.9,
            (4, 3) => 0.9,
            (4, 1) => 0.9,
            (1, 4) => 0.9,
            (1, 5) => 0.4,
            (2, 6) => 0.5,
            (3, 7) => 0.4,
            (4, 8) => 0.3
        )
        k = 2
        
        best_spread = 0.0
        for subset in combinations(1:8, k)
            spread = simulate_ic(g, weights, collect(subset); n_simulations=n_simulations)
            best_spread = max(best_spread, spread)
        end
        
        solution, spread = influence_maximization_ic(g, weights, k; n_simulations_small=n_simulations ÷ 10, n_simulations=n_simulations)
        @test spread >= 0.63 * best_spread
    end

    @testset "small graph 9: two connected components" begin
        g = SimpleDiGraph(8)
        # Component 1 (nodes 1-4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 1)
        # Component 2 (nodes 5-8)
        add_edge!(g, 5, 6)
        add_edge!(g, 6, 7)
        add_edge!(g, 7, 8)
        add_edge!(g, 8, 5)
        # Bridge
        add_edge!(g, 2, 6)
        
        weights = Dict(
            (1, 2) => 0.8,
            (2, 3) => 0.7,
            (3, 4) => 0.8,
            (4, 1) => 0.6,
            (5, 6) => 0.7,
            (6, 7) => 0.8,
            (7, 8) => 0.7,
            (8, 5) => 0.6,
            (2, 6) => 0.5
        )
        k = 3

        
        best_spread = 0.0
        for subset in combinations(1:8, k)
            spread = simulate_ic(g, weights, collect(subset); n_simulations=n_simulations)
            best_spread = max(best_spread, spread)
        end
        
        solution, spread = influence_maximization_ic(g, weights, k; n_simulations_small=n_simulations ÷ 10, n_simulations=n_simulations)
        @test spread >= 0.63 * best_spread
    end

    @testset "small graph 10: high influence probabilities" begin
        g = SimpleDiGraph(5)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        add_edge!(g, 2, 4)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 5)
        add_edge!(g, 2, 5)
        add_edge!(g, 3, 5)
        
        weights = Dict(
            (1, 2) => 0.95,
            (1, 3) => 0.95,
            (2, 4) => 0.90,
            (3, 4) => 0.90,
            (4, 5) => 0.95,
            (2, 5) => 0.85,
            (3, 5) => 0.85
        )
        k = 1
        
        best_spread = 0.0
        for subset in combinations(1:5, k)
            spread = simulate_ic(g, weights, collect(subset); n_simulations=n_simulations)
            best_spread = max(best_spread, spread)
        end
        
        solution, spread = influence_maximization_ic(g, weights, k; n_simulations_small=n_simulations ÷ 10, n_simulations=n_simulations)
        @test spread >= 0.63 * best_spread
    end

    @testset "random graph: Erdős-Rényi with n=20, p=0.2" begin
        using Random
        Random.seed!(42)  # For reproducibility
        
        n = 20
        edge_prob = 0.2
        k = 5
        
        # Generate random directed graph
        g = SimpleDiGraph(n)
        weights = Dict{Tuple{Int, Int}, Float64}()
        
        for u in 1:n
            for v in 1:n
                if u != v && rand() < edge_prob
                    add_edge!(g, u, v)
                    weights[(u, v)] = rand()  # Random weight in [0, 1]
                end
            end
        end
        
        # For random graphs of this size, brute force is infeasible (C(20,5) = 15,504)
        # So we just test that the algorithm runs and produces reasonable output
        solution, spread = influence_maximization_ic(g, weights, k; n_simulations_small=n_simulations ÷ 10, n_simulations=n_simulations)
        
        @test length(solution) <= k
        @test length(solution) == length(unique(solution))  # No duplicates
        @test all(1 <= v <= n for v in solution)  # All nodes in valid range
        @test spread >= k  # At minimum, k seed nodes are activated
        @test spread <= n  # Cannot exceed total number of nodes
        
        # Additional sanity check: single best node should give less spread than k nodes
        best_single_spread = 0.0
        for v in 1:n
            single_spread = simulate_ic(g, weights, [v]; n_simulations=n_simulations)
            best_single_spread = max(best_single_spread, single_spread)
        end
        @test spread >= best_single_spread  # k nodes should beat 1 node
    end


    @testset "random graph: Erdős-Rényi with n=40, p=0.2" begin
        using Random
        Random.seed!(42)  # For reproducibility
        
        n = 40
        edge_prob = 0.2
        k = 5
        
        # Generate random directed graph
        g = SimpleDiGraph(n)
        weights = Dict{Tuple{Int, Int}, Float64}()
        
        for u in 1:n
            for v in 1:n
                if u != v && rand() < edge_prob
                    add_edge!(g, u, v)
                    weights[(u, v)] = rand()  # Random weight in [0, 1]
                end
            end
        end
        
        # For random graphs of this size, brute force is infeasible (C(40,5) = 658,008)
        # So we just test that the algorithm runs and produces reasonable output
        solution, spread = influence_maximization_ic(g, weights, k; n_simulations_small=n_simulations ÷ 10, n_simulations=n_simulations)
        
        @test length(solution) <= k
        @test length(solution) == length(unique(solution))  # No duplicates
        @test all(1 <= v <= n for v in solution)  # All nodes in valid range
        @test spread >= k  # At minimum, k seed nodes are activated
        @test spread <= n  # Cannot exceed total number of nodes
        
        # Additional sanity check: single best node should give less spread than k nodes
        best_single_spread = 0.0
        for v in 1:n
            single_spread = simulate_ic(g, weights, [v]; n_simulations=n_simulations)
            best_single_spread = max(best_single_spread, single_spread)
        end
        @test spread >= best_single_spread  # k nodes should beat 1 node
    end

    @testset "undirected graph: w(u,v) direction takes precedence" begin
        g = SimpleGraph(2)
        add_edge!(g, 1, 2)

        # w(1,2) = 1.0 takes precedence over w(2,1) = 0.0 since 1 < 2
        weights = Dict((1, 2) => 1.0, (2, 1) => 0.0)
        spread = simulate_ic(g, weights, [1]; n_simulations=5000)
        @test spread ≈ 2.0 atol=0.1  # node 1 always activates node 2

        # k=0 returns empty solution regardless of weights
        solution, spread = influence_maximization_ic(g, weights, 0; n_simulations_small=50, n_simulations=200)
        @test solution == Int[]
        @test spread == 0.0
    end

    @testset "undirected graphs: unweighted (uniform probabilities)" begin
        @testset "path graph with prob=1: seed always reaches all nodes" begin
            g = SimpleGraph(5)
            for i in 1:4; add_edge!(g, i, i+1); end
            weights = Dict{Tuple{Int,Int}, Float64}()
            for i in 1:4
                weights[(i, i+1)] = 1.0
                weights[(i+1, i)] = 1.0
            end
            # With prob=1, any seed reaches all 5 nodes
            for seed in 1:5
                spread = simulate_ic(g, weights, [seed]; n_simulations=100)
                @test spread ≈ 5.0 atol=0.01
            end
        end

        @testset "star graph: hub seed maximizes spread" begin
            g = SimpleGraph(5)
            for i in 2:5; add_edge!(g, 1, i); end
            weights = Dict{Tuple{Int,Int}, Float64}()
            for i in 2:5
                weights[(1, i)] = 0.9
                weights[(i, 1)] = 0.9
            end
            solution, _ = influence_maximization_ic(g, weights, 1;
                n_simulations_small=500, n_simulations=3000)
            @test 1 in solution  # hub should be selected
        end

        @testset "symmetric cycle: all seeds have equal expected spread" begin
            g = SimpleGraph(4)
            for i in 1:4; add_edge!(g, i, mod1(i+1, 4)); end
            weights = Dict{Tuple{Int,Int}, Float64}()
            for i in 1:4
                weights[(i, mod1(i+1, 4))] = 0.5
                weights[(mod1(i+1, 4), i)] = 0.5
            end
            spreads = [simulate_ic(g, weights, [v]; n_simulations=5000) for v in 1:4]
            # Symmetric graph → all seeds give roughly equal spread
            @test maximum(spreads) - minimum(spreads) < 0.3
        end
    end

    @testset "undirected graphs: weighted" begin
        @testset "weighted path: high-weight edges increase spread" begin
            g = SimpleGraph(4)
            add_edge!(g, 1, 2); add_edge!(g, 2, 3); add_edge!(g, 3, 4)
            # High weight on (1,2), low on (2,3) and (3,4)
            weights = Dict(
                (1,2)=>1.0, (2,1)=>1.0,
                (2,3)=>0.1, (3,2)=>0.1,
                (3,4)=>0.1, (4,3)=>0.1
            )
            spread_from_1 = simulate_ic(g, weights, [1]; n_simulations=5000)
            spread_from_4 = simulate_ic(g, weights, [4]; n_simulations=5000)
            # Seeding node 1 reaches node 2 reliably; seeding node 4 is isolated
            @test spread_from_1 > spread_from_4
        end

        @testset "random undirected weighted graph" begin
            import Random
            Random.seed!(88)
            n = 15
            g = SimpleGraph(n)
            weights = Dict{Tuple{Int,Int}, Float64}()
            for i in 1:n, j in i+1:n
                if rand() < 0.3
                    add_edge!(g, i, j)
                    w = rand() * 0.5 + 0.1
                    weights[(i,j)] = w
                    weights[(j,i)] = w
                end
            end
            k = 3
            solution, spread = influence_maximization_ic(g, weights, k;
                n_simulations_small=500, n_simulations=2000)
            @test length(solution) <= k
            @test length(solution) == length(unique(solution))
            @test all(1 <= v <= n for v in solution)
            @test spread >= k
            @test spread <= n
        end
    end

end