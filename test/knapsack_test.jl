using Test
using OffsetArrays
using JuliAlg

# helper to validate correctness
function validate_knapsack_solution(W, weights, values, chosen_idxs, opt_value)
    # indices valid
    @test all(1 .<= chosen_idxs .<= length(weights))
    # unique indices (no duplicates)
    @test length(chosen_idxs) == length(unique(chosen_idxs))
    # total weight within capacity
    total_weight = sum(weights[i] for i in chosen_idxs; init=0)
    @test total_weight <= W
    # total value matches reported optimum
    @test sum(values[i] for i in chosen_idxs; init=0) == opt_value
end

@testset "Knapsack" begin
    @testset "exact_knapsack" begin
        @testset "classic example" begin
            W = 10
            weights = Int[2, 3, 4, 5]
            values  = Int[3, 4, 5, 6]
            opt, chosen = exact_knapsack(W, weights, values)
            @test opt == 13
            validate_knapsack_solution(W, weights, values, chosen, opt)
        end

        @testset "single item fits" begin
            W = 5
            weights = Int[3]
            values  = Int[10]
            opt, chosen = exact_knapsack(W, weights, values)
            @test opt == 10
            @test chosen == Int[1]
        end

        @testset "single item too heavy" begin
            W = 2
            weights = Int[5]
            values  = Int[10]
            opt, chosen = exact_knapsack(W, weights, values)
            @test opt == 0
            @test isempty(chosen)
        end

        @testset "all items fit" begin
            W = 100
            weights = Int[10, 20, 30]
            values  = Int[5, 10, 15]
            opt, chosen = exact_knapsack(W, weights, values)
            @test opt == 30
            @test sort(chosen) == Int[1, 2, 3]
        end

        @testset "no items fit" begin
            W = 5
            weights = Int[10, 20, 30]
            values  = Int[100, 200, 300]
            opt, chosen = exact_knapsack(W, weights, values)
            @test opt == 0
            @test isempty(chosen)
        end

        @testset "empty input" begin
            W = 10
            weights = Int[]
            values  = Int[]
            opt, chosen = exact_knapsack(W, weights, values)
            @test opt == 0
            @test isempty(chosen)
        end

        @testset "tie in value different weights" begin
            W = 10
            weights = Int[5, 8, 10]
            values  = Int[10, 10, 10]
            opt, chosen = exact_knapsack(W, weights, values)
            @test opt == 10
            validate_knapsack_solution(W, weights, values, chosen, opt)
        end

        @testset "optimal subset selection" begin
            W = 50
            weights = Int[10, 20, 30]
            values  = Int[60, 100, 120]
            opt, chosen = exact_knapsack(W, weights, values)
            @test opt == 220
            validate_knapsack_solution(W, weights, values, chosen, opt)
        end

        @testset "multiple items same efficiency" begin
            W = 15
            weights = Int[5, 5, 5, 5]
            values  = Int[10, 10, 10, 10]
            opt, chosen = exact_knapsack(W, weights, values)
            @test opt == 30
            validate_knapsack_solution(W, weights, values, chosen, opt)
        end

        @testset "greedy fails but DP succeeds" begin
            W = 10
            weights = Int[1, 2, 3, 4, 5]
            values  = Int[1, 6, 10, 16, 21]
            opt, chosen = exact_knapsack(W, weights, values)
            @test opt == 38  
            validate_knapsack_solution(W, weights, values, chosen, opt)
        end
    end

    @testset "ptas_knapsack" begin
        @testset "classic example with tight epsilon" begin
            W = 10
            weights = Int[2, 3, 4, 5]
            values  = Int[3, 4, 5, 6]
            epsilon = 0.1
            opt_approx, chosen = ptas_knapsack(W, epsilon, weights, values)
            opt_exact, _ = exact_knapsack(W, weights, values)
            @test opt_approx >= (1 - epsilon) * opt_exact
            validate_knapsack_solution(W, weights, values, chosen, opt_approx)
        end

        @testset "single item" begin
            W = 5
            weights = Int[3]
            values  = Int[10]
            epsilon = 0.2
            opt, chosen = ptas_knapsack(W, epsilon, weights, values)
            @test opt == 10
            @test chosen == Int[1]
        end

        @testset "empty input" begin
            W = 10
            weights = Int[]
            values  = Int[]
            epsilon = 0.1
            opt, chosen = ptas_knapsack(W, epsilon, weights, values)
            @test opt == 0
            @test isempty(chosen)
        end

        @testset "approximation quality various epsilon" begin
            W = 50
            weights = Int[10, 20, 30, 15, 25]
            values  = Int[60, 100, 120, 80, 110]
            opt_exact, _ = exact_knapsack(W, weights, values)
            
            for epsilon in [0.5, 0.3, 0.1, 0.05]
                opt_approx, chosen = ptas_knapsack(W, epsilon, weights, values)
                @test opt_approx >= (1 - epsilon) * opt_exact
                validate_knapsack_solution(W, weights, values, chosen, opt_approx)
            end
        end

        @testset "zero scaling fallback" begin
            W = 10
            weights = Int[1, 2, 3]
            values  = Int[1, 1, 1]
            epsilon = 10.0  # huge epsilon forces fallback
            opt, chosen = ptas_knapsack(W, epsilon, weights, values)
            opt_exact, _ = exact_knapsack(W, weights, values)
            @test opt == opt_exact
            validate_knapsack_solution(W, weights, values, chosen, opt)
        end

        @testset "negative epsilon fallback" begin
            W = 10
            weights = Int[2, 3, 4]
            values  = Int[3, 4, 5]
            epsilon = -0.1
            opt, chosen = ptas_knapsack(W, epsilon, weights, values)
            opt_exact, _ = exact_knapsack(W, weights, values)
            @test opt == opt_exact
            validate_knapsack_solution(W, weights, values, chosen, opt)
        end

        @testset "all scaled to zero fallback" begin
            W = 100
            weights = Int[10, 20, 30]
            values  = Int[1, 2, 3]
            epsilon = 0.9  # high epsilon, small values
            opt, chosen = ptas_knapsack(W, epsilon, weights, values)
            @test opt >= 0
            validate_knapsack_solution(W, weights, values, chosen, opt)
        end

        @testset "identical values" begin
            W = 20
            weights = Int[5, 10, 15]
            values  = Int[10, 10, 10]
            epsilon = 0.2
            opt_approx, chosen = ptas_knapsack(W, epsilon, weights, values)
            opt_exact, _ = exact_knapsack(W, weights, values)
            @test opt_approx >= (1 - epsilon) * opt_exact
            validate_knapsack_solution(W, weights, values, chosen, opt_approx)
        end

        @testset "large epsilon still reasonable" begin
            W = 50
            weights = Int[10, 20, 30]
            values  = Int[60, 100, 120]
            epsilon = 0.5  # 50% approximation
            opt_approx, chosen = ptas_knapsack(W, epsilon, weights, values)
            opt_exact, _ = exact_knapsack(W, weights, values)
            @test opt_approx >= (1 - epsilon) * opt_exact
            validate_knapsack_solution(W, weights, values, chosen, opt_approx)
        end
    end
end

