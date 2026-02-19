using JuMP
using HiGHS


function weighted_bipartite_matching(L::AbstractVector{Int}, R::AbstractVector{Int}, weights::Dict{Tuple{Int, Int}, Float64})
    model = Model(HiGHS.Optimizer)
    @variable(model, x[L, R], Bin)

    # ensure weights are symmetric
    for (u, v) ∈ keys(weights)
        weights[(v, u)] = weights[(u, v)]
    end

    for u ∈ L
        for v ∈ R

        end
    end

end