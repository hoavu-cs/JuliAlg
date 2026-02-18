using Graphs

"""
    Simulate the Independent Cascade (IC) model of influence spread.
    Given a directed graph `g`, edge influence probabilities `weights`,
    and an initial seed set `seed_set`, runs `num_simulations` simulations
    and returns the average number of activated nodes.
"""
function simulate_ic(
    g::AbstractGraph,
    weights::Dict{Tuple{Int, Int}, Float64},
    seed_set::Vector{Int};
    n_simulations::Int = 10_000
)

    results = Vector{Int}(undef, n_simulations)

    Threads.@threads for i ∈ 1:n_simulations
        activated = Set(seed_set)
        newly_activated = Set(seed_set)

        while !isempty(newly_activated)
            next_activated = Set{Int}()
            for u ∈ newly_activated
                for v ∈ outneighbors(g, u)
                    if v ∉ activated && rand() ≤ get(weights, (u, v), 0.0)
                        push!(next_activated, v)
                    end
                end
            end

            newly_activated = next_activated
            union!(activated, newly_activated)
        end

        results[i] = length(activated)
    end

    return sum(results) / n_simulations
end


"""
    Influence Maximization in Directed Graphs in the Independent Cascade Model.
    Given a directed graph `g`, edge influence probabilities `weights`,
    and a budget `k`, selects `k` nodes to maximize the expected spread of influence.
    n_simulations_small is used for a quick estimate of marginal gain during the greedy selection,
    while n_simulations is used for a more accurate estimate of the final spread.
    Returns a ≈ (1 - 1/e) approximate solution using a greedy algorithm.
"""
function influence_maximization_ic(
    g::AbstractGraph,
    weights::Dict{Tuple{Int, Int}, Float64},
    k::Int;
    n_simulations_small::Int = 1_000,
    n_simulations::Int = 10_000
)

    n = nv(g)
    solution = Int[]
    current_spread = 0.0
    remaining = Set(1:n)

    for _ ∈ 1:min(k, n)
        if isempty(remaining)
            break
        end

        u = -1
        Δ = -Inf

        for v ∈ remaining
            push!(solution, v)

            gain_small_sim = simulate_ic(g, weights, solution; n_simulations=n_simulations_small) - current_spread
            if gain_small_sim > Δ
                gain_large_sim = simulate_ic(g, weights, solution; n_simulations=n_simulations) - current_spread
                if gain_large_sim > Δ
                    Δ = gain_large_sim
                    u = v
                end
            end

            pop!(solution)
        end

        u == -1 && break
        delete!(remaining, u)
        push!(solution, u)
        current_spread += Δ
    end

    final_spread = simulate_ic(g, weights, solution; n_simulations=n_simulations)
    return solution, final_spread
end


"""
    simulate_ic(g::SimpleGraph, weights, seed_set, n_simulations) -> Float64

Undirected-graph overload: converts `g` to a bidirectional directed graph
(each undirected edge becomes two directed edges) and delegates to the
directed `simulate_ic`.
"""
function simulate_ic(
    g::SimpleGraph,
    weights::Dict{Tuple{Int, Int}, Float64},
    seed_set::Vector{Int};
    n_simulations::Int = 10_000
)
    dg, directed_weights = _to_bidirectional_digraph(g, weights)
    return simulate_ic(dg, directed_weights, seed_set; n_simulations=n_simulations)
end


"""
    influence_maximization_ic(g::SimpleGraph, weights, k, ...) -> (solution, spread)

Undirected-graph overload: converts `g` to a bidirectional directed graph
(each undirected edge becomes two directed edges) and delegates to the
directed `influence_maximization_ic`.
"""
function influence_maximization_ic(
    g::SimpleGraph,
    weights::Dict{Tuple{Int, Int}, Float64},
    k::Int;
    n_simulations_small::Int = 1_000,
    n_simulations::Int = 10_000
)
    dg, directed_weights = _to_bidirectional_digraph(g, weights)
    return influence_maximization_ic(dg, directed_weights, k; n_simulations_small=n_simulations_small, n_simulations=n_simulations)
end


precompile(simulate_ic,(SimpleDiGraph, Dict{Tuple{Int, Int}, Float64}, Vector{Int}, Int))
precompile(simulate_ic,(SimpleGraph, Dict{Tuple{Int, Int}, Float64}, Vector{Int}, Int))
precompile(influence_maximization_ic,(SimpleDiGraph, Dict{Tuple{Int, Int}, Float64}, Int, Int, Int))
precompile(influence_maximization_ic,(SimpleGraph, Dict{Tuple{Int, Int}, Float64}, Int, Int, Int))
