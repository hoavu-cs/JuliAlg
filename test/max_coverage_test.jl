using Test
using JuliOpt
using Combinatorics

# helper to validate a max coverage solution
function validate_max_coverage(subsets, k, selected, num_covered)
    # indices valid
    @test all(1 .<= selected .<= length(subsets))
    # no duplicate sets
    @test length(selected) == length(unique(selected))
    # at most k sets selected
    @test length(selected) <= k
    # coverage count matches
    covered = Set{Int}()
    for i in selected
        union!(covered, subsets[i])
    end
    @test length(covered) == num_covered
end

@testset "Max Coverage" begin
    @testset "basic example" begin
        subsets = [Int[1, 2, 3], Int[2, 4], Int[3, 4, 5], Int[5]]
        k = Int(2)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 5
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "k equals number of subsets" begin
        subsets = [Int[1, 2], Int[3, 4], Int[5, 6]]
        k = Int(3)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 6
        @test sort(selected) == Int[1, 2, 3]
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "k exceeds number of subsets" begin
        subsets = [Int[1, 2], Int[3]]
        k = Int(10)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 3
        @test sort(selected) == Int[1, 2]
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "k equals one" begin
        subsets = [Int[1], Int[1, 2, 3], Int[2, 3]]
        k = Int(1)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 3
        @test selected == Int[2]
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "disjoint sets" begin
        subsets = [Int[1, 2], Int[3, 4], Int[5, 6], Int[7, 8]]
        k = Int(2)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 4
        @test length(selected) == 2
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "identical sets" begin
        subsets = [Int[1, 2, 3], Int[1, 2, 3], Int[1, 2, 3]]
        k = Int(2)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 3
        @test length(selected) == 1
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "greedy picks largest first" begin
        # set 1 covers 5 elements, sets 2-4 each cover 2 unique elements
        subsets = [Int[1, 2, 3, 4, 5], Int[6, 7], Int[8, 9], Int[10, 11]]
        k = Int(2)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 7
        @test 1 ∈ selected
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "overlapping sets" begin
        subsets = [Int[1, 2, 3, 4], Int[3, 4, 5, 6], Int[5, 6, 7, 8]]
        k = Int(2)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 8
        @test sort(selected) == Int[1, 3]
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "single element per set" begin
        subsets = [Int[1], Int[2], Int[3], Int[4], Int[5]]
        k = Int(3)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 3
        @test length(selected) == 3
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "approximation quality vs brute force" begin
        subsets = [
            Int[1, 2, 3],
            Int[3, 4, 5],
            Int[5, 6, 7],
            Int[7, 8, 9],
            Int[1, 4, 7, 10],
            Int[2, 5, 8, 11],
            Int[3, 6, 9, 12]
        ]
        k = Int(3)
        m = length(subsets)

        # brute force optimal
        best_coverage = 0
        for combo in combinations(1:m, k)
            covered = Set{Int}()
            for i in combo
                union!(covered, subsets[i])
            end
            best_coverage = max(best_coverage, length(covered))
        end

        greedy_covered, selected = max_coverage(subsets, k)
        # greedy should achieve (1 - 1/e) of optimal
        @test greedy_covered >= (1 - 1/ℯ) * best_coverage - 1e-9
        validate_max_coverage(subsets, k, selected, greedy_covered)
    end

    @testset "approximation quality large instance" begin
        subsets = [
            Int[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            Int[11, 12, 13, 14, 15],
            Int[6, 7, 8, 9, 10, 11, 12],
            Int[1, 3, 5, 7, 9, 13, 15],
            Int[2, 4, 6, 8, 10, 14],
            Int[1, 2, 3, 11, 12, 13],
            Int[4, 5, 14, 15],
            Int[16, 17, 18, 19, 20]
        ]
        k = Int(3)
        m = length(subsets)

        # brute force optimal
        best_coverage = 0
        for combo in combinations(1:m, k)
            covered = Set{Int}()
            for i in combo
                union!(covered, subsets[i])
            end
            best_coverage = max(best_coverage, length(covered))
        end

        greedy_covered, selected = max_coverage(subsets, k)
        @test greedy_covered >= (1 - 1/ℯ) * best_coverage - 1e-9
        validate_max_coverage(subsets, k, selected, greedy_covered)
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

        # The first 5 sets cover everything; random sets cover only 20000 elements
        # Greedy should pick the 5 large sets first
        subsets = Vector{Vector{Int}}(undef, 1000)
        subsets[1] = s1
        subsets[2] = s2
        subsets[3] = s3
        subsets[4] = s4
        subsets[5] = s5
        for i in 6:1000
            subsets[i] = sort!(unique(rand(1:n, 20000)))
        end

        k = Int64(5)
        num_covered, selected = max_coverage(subsets, k)
        validate_max_coverage(subsets, k, selected, num_covered)
        # The first 5 sets cover all 1M elements
        @test num_covered == n
        @test sort(selected) == [1, 2, 3, 4, 5]
    end
end