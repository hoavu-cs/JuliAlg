module JuliOpt

export bin_packing,
       exact_knapsack,
       ptas_knapsack,
       weighted_interval_scheduling,
       influence_maximization_ic,
       simulate_ic

include("algorithms/combinatorial/knapsack.jl")
include("algorithms/combinatorial/bin_packing.jl")
include("algorithms/combinatorial/interval_scheduling.jl")
include("algorithms/graphs/influence_maximization.jl")

end
