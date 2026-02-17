using Graphs

"""
    average_undirected_weights(
        G::AbstractGraph,
        weights::Union{Dict{Tuple{Int, Int}, Float64}, Nothing}
    )

For weighted undirected graphs, checks whether reciprocal edge weights are symmetric.
If asymmetric weights are found, emits a warning and replaces each reciprocal pair
`(u, v)` / `(v, u)` by their average value.

For directed graphs, or when `weights === nothing`, returns `weights` unchanged.
"""
function average_undirected_weights(
    G::AbstractGraph,
    weights::Union{Dict{Tuple{Int, Int}, Float64}, Nothing}
)
    if weights === nothing || is_directed(G)
        return weights
    end

    asymmetric_found = false
    for e in edges(G)
        u = src(e)
        v = dst(e)
        w_uv = get(weights, (u, v), 0.0)
        w_vu = get(weights, (v, u), 0.0)
        if w_uv != w_vu
            asymmetric_found = true
            break
        end
    end

    if !asymmetric_found
        return weights
    end

    @warn "Undirected graph has asymmetric weights. Taking average of (u,v) and (v,u) as common weight."
    averaged_weights = Dict{Tuple{Int, Int}, Float64}()
    for e in edges(G)
        u = src(e)
        v = dst(e)
        w_uv = get(weights, (u, v), 0.0)
        w_vu = get(weights, (v, u), 0.0)
        avg_weight = (w_uv + w_vu) / 2.0
        averaged_weights[(u, v)] = avg_weight
        averaged_weights[(v, u)] = avg_weight
    end

    return averaged_weights
end
