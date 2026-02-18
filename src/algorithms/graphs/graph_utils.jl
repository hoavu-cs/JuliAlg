using Graphs

"""
    _to_bidirectional_digraph(g::SimpleGraph, weights::Dict{Tuple{Int,Int},Float64})

Convert an undirected graph to a directed graph by replacing each undirected edge
with two directed edges. Uses `weights[(u,v)]` for both directions.

The caller is responsible for ensuring `weights[(u,v)] == weights[(v,u)]` for every edge.
"""
function _to_bidirectional_digraph(
    g::SimpleGraph,
    weights::Dict{Tuple{Int, Int}, Float64}
)
    dg = SimpleDiGraph(nv(g))
    directed_weights = Dict{Tuple{Int, Int}, Float64}()

    for e in edges(g)
        u, v = src(e), dst(e)
        add_edge!(dg, u, v)
        add_edge!(dg, v, u)
        w = weights[(u, v)]
        directed_weights[(u, v)] = w
        directed_weights[(v, u)] = w
    end

    return dg, directed_weights
end

"""
    check_graph_integrity(G::AbstractGraph, weights=nothing)

Validate graph and weights consistency:
- If `weights` is provided, checks that every edge has a weight defined.
- If `G` is undirected and `weights` is provided, checks that weights are symmetric,
  i.e. `w(u,v) == w(v,u)` for every edge.

Throws an `ArgumentError` describing the first violation found.
"""
function check_graph_integrity(
    G::AbstractGraph,
    weights::Union{Dict{Tuple{Int, Int}, Float64}, Nothing} = nothing
)
    if weights !== nothing
        for e in edges(G)
            u, v = src(e), dst(e)
            if !haskey(weights, (u, v))
                throw(ArgumentError("Missing weight for edge ($u, $v)"))
            end
            if !is_directed(G) && !haskey(weights, (v, u))
                throw(ArgumentError("Missing weight for edge ($v, $u)"))
            end
        end

        if !is_directed(G)
            for e in edges(G)
                u, v = src(e), dst(e)
                w_uv = weights[(u, v)]
                w_vu = weights[(v, u)]
                if w_uv != w_vu
                    throw(ArgumentError(
                        "Asymmetric weights for undirected edge ($u, $v): " *
                        "w($u,$v)=$w_uv ≠ w($v,$u)=$w_vu"
                    ))
                end
            end
        end
    end
end
