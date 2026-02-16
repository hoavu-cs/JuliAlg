using Test
using JuliAlg
using Combinatorics

# helper to validate a set cover solution
function validate_set_cover(subsets, costs, selected, total_cost)
    # indices valid
    @test all(1 .<= selected .<= length(subsets))
    # no duplicate sets
    @test length(selected) == length(unique(selected))
    # total cost matches
    @test total_cost ≈ sum(costs[i] for i in selected)
    # all elements covered
    all_elements = Set{Int}()
    for s in subsets
        union!(all_elements, s)
    end
    covered = Set{Int}()
    for i in selected
        union!(covered, subsets[i])
    end
    @test all_elements ⊆ covered
end

@testset "Set Cover" begin
    @testset "basic example" begin
        subsets = [Int[1, 2, 3], Int[2, 4], Int[3, 4, 5], Int[5]]
        costs = Float64[3.0, 2.0, 3.0, 1.0]
        total_cost, selected = set_cover(subsets, costs)
        validate_set_cover(subsets, costs, selected, total_cost)
    end

    @testset "single set covers all" begin
        subsets = [Int[1, 2, 3], Int[1], Int[2, 3]]
        costs = Float64[1.0, 5.0, 5.0]
        total_cost, selected = set_cover(subsets, costs)
        @test total_cost == 1.0
        @test selected == Int[1]
        validate_set_cover(subsets, costs, selected, total_cost)
    end

    @testset "disjoint sets" begin
        subsets = [Int[1, 2], Int[3, 4], Int[5, 6]]
        costs = Float64[1.0, 1.0, 1.0]
        total_cost, selected = set_cover(subsets, costs)
        @test total_cost == 3.0
        @test sort(selected) == Int[1, 2, 3]
        validate_set_cover(subsets, costs, selected, total_cost)
    end

    @testset "greedy picks cost-effective sets" begin
        # set 1 covers all 4 elements for cost 4 (ratio 1.0)
        # sets 2,3 each cover 2 elements for cost 1 (ratio 2.0) and together cover all
        subsets = [Int[1, 2, 3, 4], Int[1, 2], Int[3, 4]]
        costs = Float64[4.0, 1.0, 1.0]
        total_cost, selected = set_cover(subsets, costs)
        @test total_cost == 2.0
        @test sort(selected) == Int[2, 3]
        validate_set_cover(subsets, costs, selected, total_cost)
    end

    @testset "unit costs" begin
        subsets = [Int[1, 2, 3], Int[4, 5], Int[1, 4], Int[2, 3, 5]]
        costs = Float64[1.0, 1.0, 1.0, 1.0]
        total_cost, selected = set_cover(subsets, costs)
        @test length(selected) <= 3
        validate_set_cover(subsets, costs, selected, total_cost)
    end

    @testset "identical sets picks cheapest" begin
        subsets = [Int[1, 2], Int[1, 2], Int[1, 2]]
        costs = Float64[5.0, 2.0, 3.0]
        total_cost, selected = set_cover(subsets, costs)
        @test total_cost == 2.0
        @test selected == Int[2]
    end

    @testset "large overlapping sets" begin
        subsets = [
            Int[1, 2, 3, 4, 5],
            Int[4, 5, 6, 7, 8],
            Int[7, 8, 9, 10],
            Int[1, 6, 9],
            Int[2, 3, 10]
        ]
        costs = Float64[5.0, 5.0, 4.0, 3.0, 3.0]
        total_cost, selected = set_cover(subsets, costs)
        validate_set_cover(subsets, costs, selected, total_cost)
    end

    @testset "approximation quality vs brute force" begin
        subsets = [
            Int[1, 2],
            Int[3, 4],
            Int[5, 6],
            Int[1, 3, 5],
            Int[2, 4, 6],
            Int[1, 2, 3, 4, 5, 6]
        ]
        costs = Float64[2.0, 2.0, 2.0, 3.0, 3.0, 7.0]
        m = length(subsets)

        # derive universe from subsets
        univ_set = Set{Int}()
        for s in subsets
            union!(univ_set, s)
        end

        # brute force optimal
        best_cost = Inf
        for r in 1:m
            for combo in combinations(1:m, r)
                covered = Set{Int}()
                for i in combo
                    union!(covered, subsets[i])
                end
                if univ_set ⊆ covered
                    c = sum(costs[i] for i in combo)
                    best_cost = min(best_cost, c)
                end
            end
        end

        greedy_cost, selected = set_cover(subsets, costs)
        # greedy should achieve within O(ln n) of optimal
        n = length(univ_set)
        harmonic = sum(1.0 / k for k in 1:n)
        @test greedy_cost <= harmonic * best_cost + 1e-9
        validate_set_cover(subsets, costs, selected, greedy_cost)
    end


    @testset "large instance with known optimal" begin
        using Random
        Random.seed!(123)
        n = 1_000_000

        # First 5 sets partition the universe with some overlap
        s1 = collect(1:250_000)
        s2 = collect(200_001:450_000)
        s3 = collect(400_001:650_000)
        s4 = collect(600_001:850_000)
        s5 = collect(800_001:1_000_000)

        # The first 5 sets cover everything, cost 1.0 each (ratio 20k-25k per unit)
        # Random sets cover only 20000 elements, cost 10.0 each (ratio 2000 per unit)
        # Greedy should always prefer the first 5 sets
        subsets = Vector{Vector{Int}}(undef, 1000)
        subsets[1] = s1
        subsets[2] = s2
        subsets[3] = s3
        subsets[4] = s4
        subsets[5] = s5
        for i in 6:1000
            subsets[i] = sort!(unique(rand(1:n, 20000)))
        end

        costs = fill(10.0, 1000)
        costs[1:5] .= 1.0

        total_cost, selected = set_cover(subsets, costs)
        validate_set_cover(subsets, costs, selected, total_cost)
        # Optimal is the first 5 sets at cost 5.0; greedy should find this
        @test total_cost == 5.0
        @test sort(selected) == [1, 2, 3, 4, 5]
    end

    @testset "single element per set" begin
        subsets = [Int[1], Int[2], Int[3]]
        costs = Float64[10.0, 20.0, 30.0]
        total_cost, selected = set_cover(subsets, costs)
        @test total_cost == 60.0
        @test sort(selected) == Int[1, 2, 3]
        validate_set_cover(subsets, costs, selected, total_cost)
    end
end
