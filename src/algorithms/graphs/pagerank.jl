using Graphs
using SparseArrays

"""
    pagerank(G::AbstractGraph; kwargs...)

Compute PageRank scores for vertices of graph `G`.

# Arguments
- `G`: Graph (directed or undirected)
- `damping`: Damping factor (default 0.85)
- `maxiter`: Maximum number of iterations (default 100)
- `tol`: Convergence tolerance (default 1e-6)
- `weights`: Optional dictionary mapping `(u, v)` edges to weights (default `nothing`).
             If provided, transition probabilities are proportional to edge weights.
             If `nothing`, all edges have equal weight.

# Returns
- Vector of PageRank scores (sums to 1)

# Examples
```julia
using Graphs, JuliAlg

# Unweighted directed graph
g = SimpleDiGraph(3)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 1)
r = pagerank(g)

# Weighted graph
weights = Dict((1, 2) => 2.0, (2, 3) => 1.0, (3, 1) => 1.0)
r = pagerank(g, weights=weights)
```

# Algorithm
PageRank computes the stationary distribution of a random walk with damping factor α.
With probability α, the walk follows an outgoing edge (weighted proportionally to edge weights).
With probability 1-α, it jumps to a random vertex uniformly.

For weighted graphs, transition probability from u to v is:
    w(u, v) / sum_{k ∈ outneighbors(u)} w(u, k)

Dangling nodes (no outgoing edges) are handled by redistributing their probability uniformly.
"""
function pagerank(
    G::AbstractGraph; 
    α::Float64=0.85,
    maxiter::Int=100,
    tol::Float64=1e-6,
    weights::Union{Dict{Tuple{Int, Int}, Float64}, Nothing}=nothing
)
    n = nv(G)
    r = fill(1.0 / n, n)
    out_weights = Vector{Float64}(undef, n)
    has_outgoing = fill(false, n)
    
    # For undirected graphs with weights, ensure symmetry
    processed_weights = weights
    if weights !== nothing && !is_directed(G)
        # Check for asymmetric weights and issue warning if needed
        asymmetric_found = false
        for (u, v) in keys(weights)
            if has_edge(G, u, v)  # Only check edges that exist in the graph
                w_uv = get(weights, (u, v), 0.0)
                w_vu = get(weights, (v, u), 0.0)
                if w_uv != w_vu
                    asymmetric_found = true
                    break
                end
            end
        end
        
        if asymmetric_found
            @warn "Undirected graph has asymmetric weights. Taking average of (u,v) and (v,u) as common weight."
            # Create a new symmetric weight dictionary
            processed_weights = Dict{Tuple{Int, Int}, Float64}()
            for e in edges(G)
                u = src(e)
                v = dst(e)
                w_uv = get(weights, (u, v), 0.0)
                w_vu = get(weights, (v, u), 0.0)
                avg_weight = (w_uv + w_vu) / 2.0
                processed_weights[(u, v)] = avg_weight
                processed_weights[(v, u)] = avg_weight
            end
        end
    end
    
    # Compute out-weights and identify dangling nodes
    for v in 1:n
        neighbors_out = outneighbors(G, v)
        if isempty(neighbors_out)
            has_outgoing[v] = false
            out_weights[v] = 0.0
        else
            has_outgoing[v] = true
            if processed_weights !== nothing
                total = 0.0
                for u in neighbors_out
                    total += get(processed_weights, (v, u), 0.0)
                end
                out_weights[v] = total
            else
                out_weights[v] = length(neighbors_out)
            end
        end
    end
    
    # Compute dangling node contribution (uniform redistribution)
    dangling_nodes = findall(x -> (x == false), has_outgoing)
    
    # Pull-based power iteration (each node computes its own rank from in-neighbors)
    # This is naturally parallelizable since each thread writes to a separate r_new[u]
    r_new = Vector{Float64}(undef, n)

    for iter in 1:maxiter
        # Compute dangling node contribution (sequential reduction)
        dangling_sum = 0.0
        for v in dangling_nodes
            dangling_sum += r[v]
        end
        base = (1.0 - α) / n + α * dangling_sum / n

        # Each node u pulls rank from its in-neighbors
        @Threads.threads for u in 1:n
            rank = base
            for v in inneighbors(G, u)
                if has_outgoing[v]
                    if processed_weights !== nothing
                        w = get(processed_weights, (v, u), 0.0)
                        if w > 0.0
                            rank += α * r[v] * (w / out_weights[v])
                        end
                    else
                        rank += α * r[v] / out_weights[v]
                    end
                end
            end
            r_new[u] = rank
        end

        # Check convergence using average absolute difference
        diff = 0.0
        for i in 1:n
            abs_diff = abs(r_new[i] - r[i])
            diff += abs_diff
        end
        diff /= n
        r, r_new = r_new, r

        if diff < tol
            return r
        end
    end
    
    @warn "PageRank did not converge after $maxiter iterations"
    return r
end

"""
    pagerank(G::AbstractGraph, weights::Dict{Tuple{Int, Int}, Float64}; kwargs...)

Convenience method for weighted PageRank.
"""
function pagerank(
    G::AbstractGraph,
    weights::Dict{Tuple{Int, Int}, Float64};
    α::Float64=0.85,
    maxiter::Int=100,
    tol::Float64=1e-6
)
    return pagerank(G; α=α, maxiter=maxiter, tol=tol, weights=weights)
end

# Precompile for common use cases
precompile(pagerank, (SimpleDiGraph{Int},))
precompile(pagerank, (SimpleGraph{Int},))
precompile(pagerank, (SimpleDiGraph{Int}, Dict{Tuple{Int, Int}, Float64}))

