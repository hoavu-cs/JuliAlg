using Graphs
using SparseArrays
using DataStructures

"""
    bw_centrality(G::AbstractGraph; kwargs...)

Compute betweenness centrality scores for vertices of graph `G`.

Betweenness centrality measures the fraction of shortest paths that pass through each vertex.
For a vertex v, its betweenness centrality is defined as:

    BC(v) = ∑_{s≠v≠t} σ_{st}(v) / σ_{st}

where σ_{st} is the total number of shortest paths from s to t, and σ_{st}(v) is the number
of those paths that pass through v.

# Arguments
- `G`: Graph (directed or undirected)
- `normalized`: If `true`, normalize scores by (n-1)(n-2) for directed graphs or 
                (n-1)(n-2)/2 for undirected graphs, where n is number of vertices.
                Default: `true`.
- `endpoints`: If `true`, include endpoints (s and t) in the paths. Default: `false`.
- `weights`: Optional dictionary mapping `(u, v)` edges to weights (default `nothing`).
             If provided, shortest paths are computed using edge weights (higher weight = 
             lower cost). If `nothing`, all edges have equal weight (unweighted).

# Returns
- Vector of betweenness centrality scores (one per vertex)

# Examples
```julia
using Graphs, JuliAlg

# Unweighted undirected graph
g = SimpleGraph(5)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 4, 5)
bc = bw_centrality(g)

# Weighted directed graph
g = SimpleDiGraph(4)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 1, 4)
weights = Dict((1, 2) => 2.0, (2, 3) => 1.0, (3, 4) => 3.0, (1, 4) => 5.0)
bc = bw_centrality(g, weights=weights)
```

# Algorithm
The betweenness centrality algorithm will be implemented using Brandes' algorithm,
which runs in O(nm) time for unweighted graphs and O(nm + n² log n) for weighted graphs,
where n is number of vertices and m is number of edges.

For unweighted graphs: Use BFS from each source vertex.
For weighted graphs: Use Dijkstra's algorithm from each source vertex.

The algorithm accumulates dependencies to avoid explicit path counting.
"""
function bw_centrality(
    G::AbstractGraph;
    normalized::Bool=true,
    endpoints::Bool=false,
    weights::Union{Dict{Tuple{Int, Int}, Float64}, Nothing}=nothing
)
    n = nv(G)
    BC = zeros(Float64, n)  

    if n < 2
        return BC
    end

    processed_weights = average_undirected_weights(G, weights)

    bc_lock = ReentrantLock()

    if weights === nothing
        # Unweighted graph: Use BFS
        Threads.@threads for s in vertices(G)
            P = [Int[] for _ in 1:n]  # Predecessors on shortest paths
            σ = zeros(Int, n)  # Number of shortest paths from s to each vertex
            d = fill(typemax(Int), n)  # Distance from s to each vertex
            S = Int[]
            Q = Int[]
            σ[s] = 1
            d[s] = 0
            push!(Q, s)

            while !isempty(Q)
                v = popfirst!(Q)
                push!(S, v)
                for w in neighbors(G, v)
                    if d[w] == typemax(Int)
                        push!(Q, w)
                        d[w] = d[v] + 1
                    end
                    if d[w] == d[v] + 1
                        σ[w] += σ[v]
                        push!(P[w], v)
                    end
                end
            end

            δ = zeros(Float64, n)  # Dependency scores
            # back-propagation of dependencies
            while !isempty(S)
                w = pop!(S)
                for v ∈ P[w]
                    if v != s
                        δ[v] += (σ[v] / σ[w]) * (1 + δ[w])
                    end
                end
            end

            lock(bc_lock) do
                for v in vertices(G)
                    BC[v] += δ[v]
                end
            end
        end
    else
        # Weighted graph: Use Dijkstra's algorithm
        Threads.@threads for s in vertices(G)
            P = [Int[] for _ in 1:n]  # Predecessors on shortest paths
            σ = zeros(Int, n)  # Number of shortest paths from s to each vertex
            d = fill(Inf, n)   # Distance from s to each vertex
            S = Int[]
            pq = DataStructures.PriorityQueue{Int, Float64}()

            σ[s] = 1
            d[s] = 0.0
            push!(pq, s => 0.0)

            while !isempty(pq)
                v = popfirst!(pq).first
                push!(S, v)
                for w in neighbors(G, v)
                    edge_key = (v, w)
                    if !haskey(processed_weights, edge_key)
                        continue
                    end
                    dist_vw = d[v] + processed_weights[edge_key]
                    if dist_vw < d[w]
                        d[w] = dist_vw
                        σ[w] = σ[v]
                        P[w] = [v]
                        pq[w] = dist_vw
                    elseif dist_vw ≈ d[w]
                        σ[w] += σ[v]
                        push!(P[w], v)
                    end
                end
            end

            δ = zeros(Float64, n)  # Dependency scores
            # back-propagation of dependencies
            while !isempty(S)
                w = pop!(S)
                for v ∈ P[w]
                    if v != s
                        δ[v] += (σ[v] / σ[w]) * (1 + δ[w])
                    end
                end
            end

            lock(bc_lock) do
                for v in vertices(G)
                    BC[v] += δ[v]
                end
            end
        end
    end

    if !is_directed(G)
        BC ./= 2
    end
    
    if normalized
        if n > 2
            denominator = is_directed(G) ? (n - 1) * (n - 2) : (n - 1) * (n - 2) / 2
            return BC ./ denominator
        else
            return BC
        end
    else
        return BC
    end
end

"""
    bw_centrality(G::AbstractGraph, weights::Dict{Tuple{Int, Int}, Float64}; kwargs...)

Convenience method for weighted betweenness centrality.
"""
function bw_centrality(
    G::AbstractGraph,
    weights::Dict{Tuple{Int, Int}, Float64};
    normalized::Bool=true,
    endpoints::Bool=false
)
    return bw_centrality(G; normalized=normalized, endpoints=endpoints, weights=weights)
end

# Precompile for common use cases
precompile(bw_centrality, (SimpleDiGraph{Int},))
precompile(bw_centrality, (SimpleGraph{Int},))
precompile(bw_centrality, (SimpleDiGraph{Int}, Dict{Tuple{Int, Int}, Float64}))