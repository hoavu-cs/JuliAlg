using Graphs

"""
    simulate_ic(g::AbstractGraph, weights::Dict{Tuple{Int, Int}, Float64}, seed_set::Vector{Int}; n_simulations::Int = 10_000)

Simulate the Independent Cascade (IC) model of influence spread.

# Arguments
- `g::AbstractGraph`: A directed graph representing the influence network.
- `weights::Dict{Tuple{Int, Int}, Float64}`: Dictionary mapping edges to influence probabilities.
- `seed_set::Vector{Int}`: Initial set of activated nodes (seed nodes).
- `n_simulations::Int`: Number of Monte Carlo simulations to run (default: 10,000).

# Returns
- `Float64`: The average number of activated nodes across all simulations.

# Description
The Independent Cascade model works as follows: starting from the seed set, each newly
activated node has a single chance to activate each of its neighbors with the given
probability. This process continues until no new nodes are activated.
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
    influence_maximization_ic(g::AbstractGraph, weights::Dict{Tuple{Int, Int}, Float64}, k::Int; n_simulations_small::Int = 1_000, n_simulations::Int = 10_000)

Find the k most influential nodes in a directed graph using the Independent Cascade model.

# Arguments
- `g::AbstractGraph`: A directed graph representing the influence network.
- `weights::Dict{Tuple{Int, Int}, Float64}`: Dictionary mapping edges to influence probabilities.
- `k::Int`: Number of seed nodes to select.
- `n_simulations_small::Int`: Number of simulations for quick marginal gain estimation (default: 1,000).
- `n_simulations::Int`: Number of simulations for accurate final spread estimation (default: 10,000).

# Returns
- `Tuple{Vector{Int}, Float64}`: A tuple containing:
  - `solution`: Vector of k node indices forming the seed set.
  - `final_spread`: The expected number of activated nodes.

# Description
This function uses a greedy algorithm to select k seed nodes that maximize the expected
influence spread under the Independent Cascade model. It provides a (1 - 1/e) approximation
to the optimal solution. The algorithm uses fewer simulations during the greedy selection
phase for efficiency, and more simulations for the final spread estimate.
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
    simulate_ic(g::SimpleGraph, weights::Dict{Tuple{Int, Int}, Float64}, seed_set::Vector{Int}; n_simulations::Int = 10_000)

Simulate the Independent Cascade model on an undirected graph.

# Arguments
- `g::SimpleGraph`: An undirected graph (converted to bidirectional directed graph internally).
- `weights::Dict{Tuple{Int, Int}, Float64}`: Dictionary mapping edges to influence probabilities.
- `seed_set::Vector{Int}`: Initial set of activated nodes (seed nodes).
- `n_simulations::Int`: Number of Monte Carlo simulations to run (default: 10,000).

# Returns
- `Float64`: The average number of activated nodes across all simulations.

# Description
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
    influence_maximization_ic(g::SimpleGraph, weights::Dict{Tuple{Int, Int}, Float64}, k::Int; n_simulations_small::Int = 1_000, n_simulations::Int = 10_000)

Find the k most influential nodes in an undirected graph using the Independent Cascade model.

# Arguments
- `g::SimpleGraph`: An undirected graph (converted to bidirectional directed graph internally).
- `weights::Dict{Tuple{Int, Int}, Float64}`: Dictionary mapping edges to influence probabilities.
- `k::Int`: Number of seed nodes to select.
- `n_simulations_small::Int`: Number of simulations for quick marginal gain estimation (default: 1,000).
- `n_simulations::Int`: Number of simulations for accurate final spread estimation (default: 10,000).

# Returns
- `Tuple{Vector{Int}, Float64}`: A tuple containing:
  - `solution`: Vector of k node indices forming the seed set.
  - `final_spread`: The expected number of activated nodes.

# Description
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
