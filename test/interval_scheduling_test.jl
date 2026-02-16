using Test
using JuliAlg

# helper to validate correctness
function validate_solution(starts, ends, weights, chosen_idxs, opt)
    # indices valid
    @test all(1 .<= chosen_idxs .<= length(starts))

    # non-overlapping (sorted by end time)
    chosen = collect(chosen_idxs)
    sort!(chosen, by = i -> ends[i])
    for k in 2:length(chosen)
        @test ends[chosen[k-1]] <= starts[chosen[k]]
    end

    # total weight matches optimum
    @test sum(weights[i] for i in chosen_idxs) == opt
end

@testset "Weighted Interval Scheduling" begin

    @testset "classic example" begin
        starts = Int[1, 2, 4, 6, 5, 7]
        ends   = Int[3, 5, 6, 7, 8, 9]
        wts    = Int[5, 6, 5, 4, 11, 2]

        opt, chosen = weighted_interval_scheduling(starts, ends, wts)
        @test opt == 17
        validate_solution(starts, ends, wts, chosen, opt)
    end

    @testset "single job" begin
        starts = Int[1]
        ends   = Int[2]
        wts    = Int[10]

        opt, chosen = weighted_interval_scheduling(starts, ends, wts)
        @test opt == 10
        @test chosen == Int[1]
    end

    @testset "all overlapping" begin
        starts = Int[1, 1, 1]
        ends   = Int[5, 5, 5]
        wts    = Int[2, 10, 7]

        opt, chosen = weighted_interval_scheduling(starts, ends, wts)
        @test opt == 10
        validate_solution(starts, ends, wts, chosen, opt)
    end

    @testset "touching endpoints allowed" begin
        starts = Int[1, 3, 6]
        ends   = Int[3, 6, 9]
        wts    = Int[5, 5, 5]

        opt, chosen = weighted_interval_scheduling(starts, ends, wts)
        @test opt == 15
        validate_solution(starts, ends, wts, chosen, opt)
    end

    @testset "tie case" begin
        starts = Int[1, 3, 1]
        ends   = Int[3, 5, 5]
        wts    = Int[5, 5, 10]

        opt, chosen = weighted_interval_scheduling(starts, ends, wts)
        @test opt == 10
        validate_solution(starts, ends, wts, chosen, opt)
    end
end

