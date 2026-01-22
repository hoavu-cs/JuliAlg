module JuliOpt

export bin_packing,
       exact_knapsack,
       ptas_knapsack,
       weighted_interval_scheduling

include("algorithms/knapsack.jl")
include("algorithms/bin_packing.jl")
include("algorithms/interval_scheduling.jl")

end
