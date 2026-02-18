using Graphs
using SparseArrays

"""
    pagerank(G, weights=nothing, α=0.85, maxiter=100, tol=1e-6)

Compute PageRank scores for vertices of directed graph `G`.

# Arguments
- `G`: Directed graph
- `weights`: Optional dictionary mapping `(u, v)` edges to weights (default `nothing`).
             If provided, transition probabilities are proportional to edge weights.
             If `nothing`, all edges have equal weight.
- `α`: Damping factor (default 0.85)
- `maxiter`: Maximum number of iterations (default 100)
- `tol`: Convergence tolerance (default 1e-6)

# Returns
- Vector of PageRank scores (sums to 1)

# Examples
```julia
using Graphs, JuliAlg

g = SimpleDiGraph(3)
add_edge!(g, 1, 2); add_edge!(g, 2, 3); add_edge!(g, 3, 1)

r = pagerank(g)

weights = Dict((1, 2) => 2.0, (2, 3) => 1.0, (3, 1) => 1.0)
r = pagerank(g, weights)
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
    G::AbstractGraph,
    weights::Union{Dict{Tuple{Int, Int}, Float64}, Nothing} = nothing,
    α::Float64 = 0.85,
    maxiter::Int = 100,
    tol::Float64 = 1e-6
)
    n = nv(G)
    r = fill(1.0 / n, n)
    out_weights = Vector{Float64}(undef, n)
    has_outgoing = fill(false, n)

    # Compute out-weights and identify dangling nodes
    for v in 1:n
        neighbors_out = outneighbors(G, v)
        if isempty(neighbors_out)
            has_outgoing[v] = false
            out_weights[v] = 0.0
        else
            has_outgoing[v] = true
            if weights !== nothing
                total = 0.0
                for u in neighbors_out
                    total += get(weights, (v, u), 0.0)
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

    for _ in 1:maxiter
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
                    if weights !== nothing
                        w = get(weights, (v, u), 0.0)
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
    pagerank(G::SimpleGraph, weights=nothing, α=0.85, maxiter=100, tol=1e-6)

Undirected-graph overload: converts `G` to a bidirectional directed graph
and delegates to the directed `pagerank`.
"""
function pagerank(
    G::SimpleGraph,
    weights::Union{Dict{Tuple{Int, Int}, Float64}, Nothing} = nothing,
    α::Float64 = 0.85,
    maxiter::Int = 100,
    tol::Float64 = 1e-6
)
    if weights !== nothing
        dg, directed_weights = _to_bidirectional_digraph(G, weights)
    else
        dg, directed_weights = SimpleDiGraph(G), nothing
    end
    return pagerank(dg, directed_weights, α, maxiter, tol)
end

# Precompile for common use cases
precompile(pagerank, (SimpleDiGraph{Int},))
precompile(pagerank, (SimpleGraph{Int},))
precompile(pagerank, (SimpleDiGraph{Int}, Dict{Tuple{Int, Int}, Float64}))