using Graphs

"""
    k_core_decomposition(G::AbstractGraph)

    Returns a Dict mapping each vertex to its core number.
    The core number of a vertex v is the largest k such that v belongs to a k-core,
    where the k-core is the maximal subgraph in which every vertex has degree at least k.
"""
function k_core_decomposition(G::AbstractGraph)
    if nv(G) == 0
        return Dict{Int, Int}()
    end

    k = 0
    H = copy(G)
    core = Dict{Int, Int}(v => degree(H, v) for v in vertices(H))
    Δ = maximum(degree(H))
    B = Dict(d => Set(v for v in vertices(H) if degree(H, v) == d) for d in 0:Δ)
    remaining = nv(H)

    while remaining > 0
        k > Δ && break
        while isempty(B[k])
            k += 1
            k > Δ && break
        end
        k > Δ && break

        v = pop!(B[k])
        for u in collect(neighbors(H, v))
            if degree(H, u) > k
                du = degree(H, u)
                rem_edge!(H, u, v)
                delete!(B[du], u)
                push!(B[du - 1], u)
            end
        end
        core[v] = k
        remaining -= 1
    end

    return core
end

precompile(k_core_decomposition, (SimpleGraph{Int},))