using Test
using JuliOpt

@testset "JuliOpt" begin
    include("knapsack_test.jl")
    include("bin_packing_test.jl")
    include("interval_scheduling_test.jl")
    include("influence_maximization_test.jl")
    include("set_cover_test.jl")
    include("max_coverage_test.jl")
    include("densest_subgraph_test.jl")
end
