using JuMP
using HiGHS
using Graphs

"""
    weighted_bipartite_matching(G::SimpleGraph, L::Vector{Int}, R::Vector{Int}; weights::Dict{Tuple{Int, Int}, Float64} = Dict{Tuple{Int, Int}, Float64}())

Find the maximum weight matching in a bipartite graph using a linear programming formulation.

# Arguments
- `G::SimpleGraph`: An undirected graph representing the bipartite graph.
- `L::Vector{Int}`: Vector of vertices in the left partition.
- `R::Vector{Int}`: Vector of vertices in the right partition.
- `weights::Dict{Tuple{Int, Int}, Float64}`: Optional dictionary of edge weights. If not provided, 
  all edges default to weight 1.0. Weights are symmetric - `(u, v)` and `(v, u)` have the same weight.

# Returns
- `Tuple{Float64, Vector{Tuple{Int, Int}}}`: A tuple containing:
  - The total weight of the matching (maximum total weight achieved).
  - A vector of tuples representing the edges in the matching.

# Throws
- Error if `L` and `R` are not disjoint or do not cover all vertices in `G`.

# Example
```julia
using Graphs

# Create a simple bipartite graph with 4 vertices
# L = [1, 2], R = [3, 4]
g = SimpleGraph(4)
add_edge!(g, 1, 3)
add_edge!(g, 1, 4)
add_edge!(g, 2, 3)
add_edge!(g, 2, 4)

# Find maximum weight matching with custom weights
weights = Dict((1, 3) => 1.0, (1, 4) => 2.0, (2, 3) => 2.0, (2, 4) => 1.0)
total_weight, matching = weighted_bipartite_matching(g, [1, 2], [3, 4]; weights=weights)
```
"""
function weighted_bipartite_matching(
    G::SimpleGraph, L::Vector{Int}, 
    R::Vector{Int}; 
    weights::Dict{Tuple{Int, Int}, Float64} = Dict{Tuple{Int, Int}, Float64}()) 

    # check L and R are disjoint and cover all vertices
    n = nv(G)
    Lset = Set(L)
    Rset = Set(R)

    if length(intersect(Lset, Rset)) > 0 || length(union(Lset, Rset)) != n
        error("L and R must be disjoint and cover all vertices")
    end

    # loop through edges
    w = copy(weights)
    for e in edges(G)
        u, v = src(e), dst(e)
        y = get(weights, (u, v), get(weights, (v, u), 1.0)) # check both directions; default 1.0
        w[(u, v)] = y   
        w[(v, u)] = y
    end

    E = Set{Tuple{Int, Int}}()
    for u ∈ Lset
        for v ∈ neighbors(G, u)
            if v ∈ Rset
                push!(E, (u, v))
            end
        end
    end

    isempty(E) && return 0.0, Tuple{Int,Int}[]

    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, 0 ≤ x[E] ≤ 1)

    # matching constraints
    for v ∈ vertices(G)
        if v ∈ Lset
            @constraint(model, sum(x[(v, u)] for u in neighbors(G, v) if (v, u) in E) ≤ 1)
        else
            @constraint(model, sum(x[(u, v)] for u in neighbors(G, v) if (u, v) in E) ≤ 1)
        end
    end

    @objective(model, Max, sum(w[(u, v)] * x[(u, v)] for (u, v) in E))
    optimize!(model)

    return objective_value(model), [(u, v) for (u, v) in E if value(x[(u, v)]) ≥ 1- 1e-6]
end

precompile(weighted_bipartite_matching, (SimpleGraph, Vector{Int}, Vector{Int}))