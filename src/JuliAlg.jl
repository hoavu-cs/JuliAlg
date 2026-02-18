module JuliAlg

export check_graph_integrity,
       bin_packing,
       exact_knapsack,
       ptas_knapsack,
       weighted_interval_scheduling,
       influence_maximization_ic,
       simulate_ic,
       set_cover,
       max_coverage,
       densest_subgraph,
       k_core_decomposition,
       pagerank,
       bw_centrality

include("algorithms/combinatorial/knapsack.jl")
include("algorithms/combinatorial/bin_packing.jl")
include("algorithms/combinatorial/interval_scheduling.jl")
include("algorithms/graphs/graph_utils.jl")
include("algorithms/graphs/influence_maximization.jl")
include("algorithms/combinatorial/set_cover.jl")
include("algorithms/combinatorial/max_coverage.jl")
include("algorithms/graphs/densest_subgraph.jl")
include("algorithms/graphs/k_core.jl")
include("algorithms/graphs/pagerank.jl")
include("algorithms/graphs/centrality.jl")

end
