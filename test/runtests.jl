using Test
using JuliAlg

@testset "JuliAlg" begin
    include("knapsack_test.jl")
    include("bin_packing_test.jl")
    include("interval_scheduling_test.jl")
    include("influence_maximization_test.jl")
    include("set_cover_test.jl")
    include("max_coverage_test.jl")
    include("densest_subgraph_test.jl")
    include("k_core_test.jl")
    include("pagerank_test.jl")
    include("centrality_test.jl")
end
