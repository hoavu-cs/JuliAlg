using Test
using Graphs
using Random
using JuliAlg

# ---------------------------------------------------------------------------
# Brute-force reference implementation
# Enumerate every subset of cross-edges and return the max-weight valid matching.
# ---------------------------------------------------------------------------
function bf_max_weight_matching(
    G::AbstractGraph,
    L::Vector{Int},
    R::Vector{Int},
    weights::Dict{Tuple{Int,Int},Float64}
)
    Lset = Set(L)
    Rset = Set(R)
    cross_edges = [(u, v) for u in L for v in neighbors(G, u) if v in Rset]
    n_edges = length(cross_edges)

    best_val = 0.0
    best_matching = Tuple{Int,Int}[]

    for mask in 0:(1 << n_edges - 1)
        subset = [cross_edges[i] for i in 1:n_edges if (mask >> (i - 1)) & 1 == 1]
        matched_L = Set{Int}()
        matched_R = Set{Int}()
        valid = true
        for (u, v) in subset
            if u in matched_L || v in matched_R
                valid = false; break
            end
            push!(matched_L, u); push!(matched_R, v)
        end
        if valid
            w = sum(
                get(weights, (u, v), get(weights, (v, u), 1.0))
                for (u, v) in subset;
                init=0.0
            )
            if w > best_val
                best_val = w
                best_matching = copy(subset)
            end
        end
    end
    return best_val, best_matching
end

# Validity checker: no vertex appears twice in the matching.
function is_valid_matching(matching, Lset, Rset)
    matched_L = Set{Int}()
    matched_R = Set{Int}()
    for (u, v) in matching
        u in Lset && v in Rset || return false
        u in matched_L && return false
        v in matched_R && return false
        push!(matched_L, u); push!(matched_R, v)
    end
    return true
end

# ---------------------------------------------------------------------------

@testset "Weighted Bipartite Matching Tests" begin

    @testset "Input validation" begin
        @testset "L and R overlap" begin
            g = SimpleGraph(3)
            add_edge!(g, 1, 3)
            @test_throws ErrorException weighted_bipartite_matching(
                g, [1, 2], [2, 3]
            )
        end

        @testset "L ∪ R does not cover all vertices" begin
            g = SimpleGraph(4)
            add_edge!(g, 1, 3)
            @test_throws ErrorException weighted_bipartite_matching(
                g, [1, 2], [3]  # vertex 4 missing
            )
        end
    end

    @testset "No cross edges (empty matching)" begin
        # L and R exist but there are no edges between them
        g = SimpleGraph(4)
        add_edge!(g, 1, 2)   # intra-L edge, ignored
        add_edge!(g, 3, 4)   # intra-R edge, ignored
        L = [1, 2]; R = [3, 4]
        weights = Dict{Tuple{Int,Int},Float64}()

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        @test val ≈ 0.0 atol=1e-6
        @test isempty(matching)
    end

    @testset "Single cross edge" begin
        g = SimpleGraph(2)
        add_edge!(g, 1, 2)
        L = [1]; R = [2]
        weights = Dict((1, 2) => 7.5)

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        @test val ≈ 7.5 atol=1e-6
        @test length(matching) == 1
        @test (1, 2) in matching
    end

    @testset "Weight competition for a single R-vertex" begin
        # L={1}, R={2,3}: both L→2 and L→3 available; pick higher weight
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        L = [1]; R = [2, 3]
        weights = Dict((1, 2) => 3.0, (1, 3) => 7.0)

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        @test val ≈ 7.0 atol=1e-6
        @test length(matching) == 1
        @test (1, 3) in matching
    end

    @testset "Both L-vertices matched to distinct R-vertices" begin
        g = SimpleGraph(4)
        add_edge!(g, 1, 3)
        add_edge!(g, 2, 4)
        L = [1, 2]; R = [3, 4]
        weights = Dict((1, 3) => 5.0, (2, 4) => 3.0)

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        @test val ≈ 8.0 atol=1e-6
        @test Set(matching) == Set([(1, 3), (2, 4)])
    end

    @testset "Classic conflict: shared R-vertex forces trade-off" begin
        # L={1,2}, R={3,4}: edges (1,3)=10, (2,3)=8, (2,4)=3
        # Best: match (1,3)=10 and (2,4)=3 = 13, better than (2,3)=8 alone
        g = SimpleGraph(4)
        add_edge!(g, 1, 3); add_edge!(g, 2, 3); add_edge!(g, 2, 4)
        L = [1, 2]; R = [3, 4]
        weights = Dict((1, 3) => 10.0, (2, 3) => 8.0, (2, 4) => 3.0)

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        @test val ≈ 13.0 atol=1e-6
        @test Set(matching) == Set([(1, 3), (2, 4)])
    end

    @testset "K_{2,2}: complete bipartite, optimal assignment" begin
        # L={1,2}, R={3,4}
        # weights: (1,3)=4, (1,4)=1, (2,3)=2, (2,4)=6
        # Options: {(1,3),(2,4)} = 10  vs  {(1,4),(2,3)} = 3
        g = SimpleGraph(4)
        add_edge!(g, 1, 3); add_edge!(g, 1, 4)
        add_edge!(g, 2, 3); add_edge!(g, 2, 4)
        L = [1, 2]; R = [3, 4]
        weights = Dict((1, 3) => 4.0, (1, 4) => 1.0, (2, 3) => 2.0, (2, 4) => 6.0)

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        @test val ≈ 10.0 atol=1e-6
        @test Set(matching) == Set([(1, 3), (2, 4)])
        @test is_valid_matching(matching, Set(L), Set(R))
    end

    @testset "K_{3,3}: classic 3×3 assignment problem" begin
        # L={1,2,3}, R={4,5,6}
        # Weight matrix:
        #       4    5    6
        #  1 [  3    2    1 ]
        #  2 [  6    4    5 ]
        #  3 [  1    7    2 ]
        # Perfect matchings and their total weights:
        #  (1→4,2→5,3→6) = 3+4+2 = 9
        #  (1→4,2→6,3→5) = 3+5+7 = 15  ← optimal
        #  (1→5,2→4,3→6) = 2+6+2 = 10
        #  (1→5,2→6,3→4) = 2+5+1 = 8
        #  (1→6,2→4,3→5) = 1+6+7 = 14
        #  (1→6,2→5,3→4) = 1+4+1 = 6
        g = SimpleGraph(6)
        for u in 1:3, v in 4:6; add_edge!(g, u, v); end
        L = [1, 2, 3]; R = [4, 5, 6]
        weights = Dict(
            (1, 4) => 3.0, (1, 5) => 2.0, (1, 6) => 1.0,
            (2, 4) => 6.0, (2, 5) => 4.0, (2, 6) => 5.0,
            (3, 4) => 1.0, (3, 5) => 7.0, (3, 6) => 2.0,
        )

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        @test val ≈ 15.0 atol=1e-6
        @test Set(matching) == Set([(1, 4), (2, 6), (3, 5)])
        @test is_valid_matching(matching, Set(L), Set(R))
    end

    @testset "Unbalanced bipartite (|L| > |R|)" begin
        # L={1,2,3}, R={4,5}: at most 2 can be matched
        g = SimpleGraph(5)
        add_edge!(g, 1, 4); add_edge!(g, 2, 4); add_edge!(g, 3, 5)
        L = [1, 2, 3]; R = [4, 5]
        weights = Dict((1, 4) => 1.0, (2, 4) => 9.0, (3, 5) => 5.0)

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        # Optimal: (2,4)=9 + (3,5)=5 = 14
        @test val ≈ 14.0 atol=1e-6
        @test Set(matching) == Set([(2, 4), (3, 5)])
        @test is_valid_matching(matching, Set(L), Set(R))
    end

    @testset "Empty weights dict — default weight 1.0 (max cardinality)" begin
        # K_{2,3}: L={1,2}, R={3,4,5}; all edges weight 1 by default
        g = SimpleGraph(5)
        for u in 1:2, v in 3:5; add_edge!(g, u, v); end
        L = [1, 2]; R = [3, 4, 5]
        weights = Dict{Tuple{Int,Int},Float64}()

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        # Maximum matching has 2 edges (limited by |L|)
        @test val ≈ 2.0 atol=1e-6
        @test length(matching) == 2
        @test is_valid_matching(matching, Set(L), Set(R))
    end

    @testset "Partial weights — unspecified edges get default 1.0" begin
        g = SimpleGraph(4)
        add_edge!(g, 1, 3); add_edge!(g, 1, 4); add_edge!(g, 2, 3)
        L = [1, 2]; R = [3, 4]
        # Only specify (1,4)=5.0; others default to 1.0
        weights = Dict((1, 4) => 5.0)

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        # Options: (1,4)=5 + (2,3)=1 = 6  vs  (1,3)=1 alone = 1
        @test val ≈ 6.0 atol=1e-6
        @test Set(matching) == Set([(1, 4), (2, 3)])
    end

    @testset "Weight given as (R,L) key — symmetrization" begin
        # User provides weight under the reverse-direction key (v, u)
        g = SimpleGraph(2)
        add_edge!(g, 1, 2)
        L = [1]; R = [2]
        # Weight stored with R-vertex first
        weights = Dict((2, 1) => 9.0)

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        @test val ≈ 9.0 atol=1e-6
        @test length(matching) == 1
    end

    @testset "Intra-partition edges are ignored" begin
        # Edge (1,2) is within L, edge (3,4) is within R; only (1,3) is cross
        g = SimpleGraph(4)
        add_edge!(g, 1, 2)   # L-L edge, ignored
        add_edge!(g, 3, 4)   # R-R edge, ignored
        add_edge!(g, 1, 3)   # L-R edge
        L = [1, 2]; R = [3, 4]
        weights = Dict((1, 3) => 4.0)

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        @test val ≈ 4.0 atol=1e-6
        @test Set(matching) == Set([(1, 3)])
    end

    @testset "Returned matching validity" begin
        Random.seed!(7)
        n = 8
        g = SimpleGraph(n)
        L = collect(1:4); R = collect(5:8)
        weights = Dict{Tuple{Int,Int},Float64}()
        for u in L, v in R
            if rand() < 0.6
                add_edge!(g, u, v)
                weights[(u, v)] = rand() * 10.0
            end
        end

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        # Matching must be valid
        @test is_valid_matching(matching, Set(L), Set(R))
        # Every returned edge must exist in G
        for (u, v) in matching
            @test has_edge(g, u, v)
        end
        # Returned value must equal the sum of matched edge weights
        computed = sum(get(weights, (u, v), get(weights, (v, u), 1.0)) for (u, v) in matching; init=0.0)
        @test val ≈ computed atol=1e-6
    end

    @testset "Brute force comparison: random K_{3,4}" begin
        Random.seed!(42)
        g = SimpleGraph(7)
        L = [1, 2, 3]; R = [4, 5, 6, 7]
        weights = Dict{Tuple{Int,Int},Float64}()
        for u in L, v in R
            add_edge!(g, u, v)
            weights[(u, v)] = round(rand() * 10.0, digits=2)
        end

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        bf_val, _ = bf_max_weight_matching(g, L, R, weights)

        @test val ≈ bf_val atol=1e-6
        @test is_valid_matching(matching, Set(L), Set(R))
    end

    @testset "Brute force comparison: sparse random bipartite" begin
        Random.seed!(99)
        n_L = 4; n_R = 4
        g = SimpleGraph(n_L + n_R)
        L = collect(1:n_L); R = collect(n_L+1:n_L+n_R)
        weights = Dict{Tuple{Int,Int},Float64}()
        for u in L, v in R
            if rand() < 0.5
                add_edge!(g, u, v)
                weights[(u, v)] = rand() * 10.0
            end
        end

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        bf_val, _ = bf_max_weight_matching(g, L, R, weights)

        @test val ≈ bf_val atol=1e-6
        @test is_valid_matching(matching, Set(L), Set(R))
    end

    @testset "Larger balanced K_{5,5}" begin
        g = SimpleGraph(10)
        L = collect(1:5); R = collect(6:10)
        for u in L, v in R; add_edge!(g, u, v); end
        # Diagonal assignment: (i, 5+i) has weight i*10, all others weight 1
        weights = Dict{Tuple{Int,Int},Float64}()
        for u in L, v in R
            weights[(u, v)] = (v == u + 5) ? Float64(u * 10) : 1.0
        end
        # Optimal: match diagonal (1,6)=10,(2,7)=20,(3,8)=30,(4,9)=40,(5,10)=50 = 150
        # vs any other permutation that loses the diagonal gets far less

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        @test val ≈ 150.0 atol=1e-6
        @test Set(matching) == Set([(i, i + 5) for i in 1:5])
        @test is_valid_matching(matching, Set(L), Set(R))
    end

    @testset "Custom graph 1" begin
        g = SimpleGraph(7)
        add_edge!(g, 1, 4)
        add_edge!(g, 1, 6)
        add_edge!(g, 2, 4)
        add_edge!(g, 2, 7)
        add_edge!(g, 3, 5)
        add_edge!(g, 3, 7)
        L = [1, 2, 3]; 
        R = [4, 5, 6, 7]
        weights = Dict((1, 4) => 0.5, (1, 6) => 0.3, (2, 4) => 0.4, (2, 7) => 0.1, (3, 5) => 0.7, (3, 7) => 0.2)

        val, matching = weighted_bipartite_matching(g, L, R; weights)
        @test val ≈ 0.4 + 0.3 + 0.7 atol=1e-6
        @test Set(matching) == Set([(1, 6), (2, 4), (3, 5)])
    end

end