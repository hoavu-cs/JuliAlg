using Graphs, GraphsFlows
using Combinatorics
using SparseArrays

"""
    Create auxiliary graph H for the Goldberg algorithm.
"""
function create_aux_graph(G::AbstractGraph, λ::Float64)
    n = nv(G)
    N = n + 2
    s = n + 1
    t = n + 2

    H = DiGraph(N)
    cap = spzeros(Float64, N, N)

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

    Returns (best_S, best_density)
    - best_S: vertex ids (1..n) of the best subgraph found
    - best_density: density(best_S) = |E(S)|/|S|
"""
function densest_subgraph(G::AbstractGraph, num_iterations::Int = 40, algorithm=:goldberg)
    if is_directed(G)
        throw(ArgumentError("densest_subgraph only supports undirected graphs"))
    end

    n = nv(G)
    m = ne(G)

    low = 0.0
    high = maximum(degree(G)) / 2.0

    best_S = collect(1:n)
    best_λ = 0.0

    while high - low ≥ 1/(n * (n - 1))
        mid = (low + high) / 2.0

        H, cap, s, t = create_aux_graph(G, mid)

        part_s, part_t, cut_value = GraphsFlows.mincut(H, s, t, cap, PushRelabelAlgorithm())
        S = [v for v in part_s if 1 ≤ v ≤ n]

        if cut_value ≤ 2.0 * m + 1e-9
            low = mid
            best_λ = mid
            if !isempty(S)
                best_S = S
            end
        else
            high = mid
        end
    end

    return best_S, density(G, best_S)
end

"""
    Charikar's peeling algorithm for densest subgraph. 
    1/2-approximation.
"""
function densest_subgraph_peeling(G::AbstractGraph)
    if is_directed(G)
        throw(ArgumentError("densest_subgraph_peeling only supports undirected graphs"))
    end

    H = copy(G)
    n = nv(H)

    active = Set(1:n)
    best_S = copy(active)
    best_density = density(H, collect(best_S))
    remaining = n

    Δ = maximum(degree(H))
    B = Dict(d => Set(v for v in vertices(H) if degree(H, v) == d) for d in 0:Δ)
    d = 0

    while remaining > 0
        current_density = ne(H) / remaining
        if current_density > best_density
            best_density = current_density
            best_S = copy(active)
        end
        
        # Find the next minimum degree d
        d > Δ && break
        while isempty(B[d])
            d += 1
            d > Δ && break
        end
        d > Δ && break

        v = pop!(B[d]) # Remove a minimum degree vertex v
        for u in collect(neighbors(H, v))
            du = degree(H, u)
            rem_edge!(H, u, v)
            delete!(B[du], u) 
            push!(B[du - 1], u)
        end
        delete!(active, v)
        d = max(d-1, 0)
        remaining -= 1
    end

    return best_S, best_density
end

"""
    Densest at-most-k-subgraph problem: find a subset S of vertices with |S| ≤ k that maximizes density(S) = |E(S)|/|S|.
    We try all subsets of size ≤ k and pick the one with the highest density. 
    However, one can prune the original graph and speed up the search as follows.

    Iteratively remove vertices whose current degree is smallest until only k vertices remain. Let the remaining graph be H.
    Starting with d' = density(H), we iteratively remove nodes with the smallest degree in H until no node remains.
    In the process, we update d' to keep track of the maximum density seen so far. 
    We iteratively remove vertices whose current degree is strictly smallerthan d′, until no such vertex remains.
    Claim: Any vertex removed by this procedure cannot belong to an optimal densest
    subgraph H⋆.
"""
function densest_at_most_k_subgraph(G::AbstractGraph, k::Int)
    if is_directed(G)
        throw(ArgumentError("densest_at_most_k_subgraph only supports undirected graphs"))
    end

    n = nv(G)
    if k ≥ n
        return densest_subgraph(G)
    end
    
    # Step 1: prune the graphs by iteratively remove the lowest degree vertices
    # Store the best snapshot density once the number of vertices is at most k in dprime
    H = copy(G)
    Δ = maximum(degree(H))
    B = Dict(d => Set(v for v in vertices(H) if degree(H, v) == d) for d in 0:Δ)
    d = 0
    dprime = -Inf
    remaining = n

    while remaining > 0
        if remaining ≤ k
            dprime = max(dprime, ne(H) / remaining)
        end

        # Find the next minimum degree d 
        d > Δ && break
        while isempty(B[d])
            d += 1
            d > Δ && break
        end
        d > Δ && break

        v = pop!(B[d]) # Remove a minimum degree vertex v
        for u in collect(neighbors(H, v))
            du = degree(H, u)
            rem_edge!(H, u, v)
            delete!(B[du], u) 
            push!(B[du - 1], u)
        end

        d = max(d-1, 0)
        remaining -= 1
    end

    # Step 2: Remove nodes with degree < dprime until no such node remains
    # Remove nodes will be isolated in the remaining graph. 
    H = copy(G)
    while nv(H) > 0
        removes = [v for v in vertices(H) if 0 < degree(H, v) < dprime]
        isempty(removes) && break
        for v in removes
            for u in collect(neighbors(H, v))
                rem_edge!(H, u, v)
            end
        end
    end

    # Step 3: Brute force search on the remaining graph H with at most k vertices
    best_S = Int[]
    best_density = 0.0
    H = [v for v in vertices(H) if degree(H, v) > 0]
    for size in 1:min(k, length(H))
        for S in combinations(H, size)
            dS = density(G, S)
            if dS > best_density
                best_density = dS
                best_S = S
            end
        end
    end

    return best_S, best_density
end

precompile(create_aux_graph, (SimpleGraph{Int}, Float64))
precompile(density, (SimpleGraph{Int}, Vector{Int}))
precompile(densest_subgraph, (SimpleGraph{Int}, Int, Symbol))
precompile(densest_at_most_k_subgraph, (SimpleGraph{Int}, Int))