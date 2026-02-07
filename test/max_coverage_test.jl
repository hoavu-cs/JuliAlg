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
    covered = Set{Int64}()
    for i in selected
        union!(covered, subsets[i])
    end
    @test length(covered) == num_covered
end

@testset "Max Coverage" begin
    @testset "basic example" begin
        subsets = [Int64[1, 2, 3], Int64[2, 4], Int64[3, 4, 5], Int64[5]]
        k = Int64(2)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 5
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "k equals number of subsets" begin
        subsets = [Int64[1, 2], Int64[3, 4], Int64[5, 6]]
        k = Int64(3)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 6
        @test sort(selected) == Int64[1, 2, 3]
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "k exceeds number of subsets" begin
        subsets = [Int64[1, 2], Int64[3]]
        k = Int64(10)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 3
        @test sort(selected) == Int64[1, 2]
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "k equals one" begin
        subsets = [Int64[1], Int64[1, 2, 3], Int64[2, 3]]
        k = Int64(1)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 3
        @test selected == Int64[2]
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "disjoint sets" begin
        subsets = [Int64[1, 2], Int64[3, 4], Int64[5, 6], Int64[7, 8]]
        k = Int64(2)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 4
        @test length(selected) == 2
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "identical sets" begin
        subsets = [Int64[1, 2, 3], Int64[1, 2, 3], Int64[1, 2, 3]]
        k = Int64(2)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 3
        @test length(selected) == 1
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "greedy picks largest first" begin
        # set 1 covers 5 elements, sets 2-4 each cover 2 unique elements
        subsets = [Int64[1, 2, 3, 4, 5], Int64[6, 7], Int64[8, 9], Int64[10, 11]]
        k = Int64(2)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 7
        @test 1 ∈ selected
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "overlapping sets" begin
        subsets = [Int64[1, 2, 3, 4], Int64[3, 4, 5, 6], Int64[5, 6, 7, 8]]
        k = Int64(2)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 8
        @test sort(selected) == Int64[1, 3]
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "single element per set" begin
        subsets = [Int64[1], Int64[2], Int64[3], Int64[4], Int64[5]]
        k = Int64(3)
        num_covered, selected = max_coverage(subsets, k)
        @test num_covered == 3
        @test length(selected) == 3
        validate_max_coverage(subsets, k, selected, num_covered)
    end

    @testset "approximation quality vs brute force" begin
        subsets = [
            Int64[1, 2, 3],
            Int64[3, 4, 5],
            Int64[5, 6, 7],
            Int64[7, 8, 9],
            Int64[1, 4, 7, 10],
            Int64[2, 5, 8, 11],
            Int64[3, 6, 9, 12]
        ]
        k = Int64(3)
        m = length(subsets)

        # brute force optimal
        best_coverage = 0
        for combo in combinations(1:m, k)
            covered = Set{Int64}()
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
            Int64[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            Int64[11, 12, 13, 14, 15],
            Int64[6, 7, 8, 9, 10, 11, 12],
            Int64[1, 3, 5, 7, 9, 13, 15],
            Int64[2, 4, 6, 8, 10, 14],
            Int64[1, 2, 3, 11, 12, 13],
            Int64[4, 5, 14, 15],
            Int64[16, 17, 18, 19, 20]
        ]
        k = Int64(3)
        m = length(subsets)

        # brute force optimal
        best_coverage = 0
        for combo in combinations(1:m, k)
            covered = Set{Int64}()
            for i in combo
                union!(covered, subsets[i])
            end
            best_coverage = max(best_coverage, length(covered))
        end

        greedy_covered, selected = max_coverage(subsets, k)
        @test greedy_covered >= (1 - 1/ℯ) * best_coverage - 1e-9
        validate_max_coverage(subsets, k, selected, greedy_covered)
    end
end