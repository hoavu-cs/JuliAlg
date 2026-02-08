using Graphs, GraphsFlows

"""
    Create auxiliary graph H for the Goldberg algorithm.
"""
function create_aux_graph(G::AbstractGraph, λ::Float64)
    n = nv(G)
    N = n + 2
    s = n + 1
    t = n + 2

    H = DiGraph(N)
    cap = zeros(Float64, N, N)

    for v ∈ 1:n
        add_edge!(H, s, v)
        cap[s, v] = degree(G, v)
        add_edge!(H, v, t)           
        cap[v, t] = 2.0 * λ
    end

    for e ∈ edges(G)
        u, v = src(e), dst(e)
        add_edge!(H, u, v); cap[u, v] = 1.0
        add_edge!(H, v, u); cap[v, u] = 1.0
    end

    return H, cap, s, t
end

"""
    Compute the density of a subgraph S in G, defined as |E(S)|/|S|.
"""
function density(G::AbstractGraph, S::Vector{Int})
    if isempty(S)
        return 0.0
    end
    Sset = Set(S)

    eS = 0  
    for v ∈ S
        for u ∈ neighbors(G, v)
            if u ∈ Sset
                eS += 1
            end
        end
    end

    eS = eS / 2.0  
    return eS / length(S)
end

"""
Goldberg algorithm for densest subgraph. 

Returns (bestS, bestλ, best_density_est)
- bestS: vertex ids (1..n) of the best subgraph found
- bestλ: last feasible λ from the search
- best_density_est: density(bestS) = |E(S)|/|S|
"""
function densest_subgraph(G::AbstractGraph, num_iterations::Int = 40, algorithm=:goldberg)
    n = nv(G)
    m = ne(G)

    low = 0.0
    high = maximum(degree(G)) / 2.0

    bestS = collect(1:n)
    bestλ = 0.0

    for _ ∈ 1:num_iterations
        mid = (low + high) / 2.0

        H, cap, s, t = create_aux_graph(G, mid)

        part_S, part_T, cut_value = GraphsFlows.mincut(H, s, t, cap, PushRelabelAlgorithm())
        S = [v for v in part_S if 1 ≤ v ≤ n]

        if cut_value ≤ 2.0 * m + 1e-9
            low = mid
            bestλ = mid
            if !isempty(S)
                bestS = S
            end
        else
            high = mid
        end
    end

    return bestS, bestλ, density(G, bestS)
end

precompile(create_aux_graph, (SimpleGraph{Int}, Float64))
precompile(density, (SimpleGraph{Int}, Vector{Int}))
precompile(densest_subgraph, (SimpleGraph{Int}, Int, Symbol))
