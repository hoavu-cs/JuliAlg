using JuMP
using HiGHS
using Graphs

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
            @constraint(model, sum(x[(v, u)] for u in neighbors(G, v) if (v, u) in E) <= 1)
        else
            @constraint(model, sum(x[(u, v)] for u in neighbors(G, v) if (u, v) in E) <= 1)
        end
    end

    @objective(model, Max, sum(w[(u, v)] * x[(u, v)] for (u, v) in E))
    optimize!(model)

    return objective_value(model), [(u, v) for (u, v) in E if value(x[(u, v)]) ≥ 1- 1e-6]
end

precompile(weighted_bipartite_matching, (SimpleGraph, Vector{Int}, Vector{Int}))