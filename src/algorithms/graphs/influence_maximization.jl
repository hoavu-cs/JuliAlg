using Graphs

"""
    Simulate the Independent Cascade (IC) model of influence spread.
    Given a directed graph `g`, edge influence probabilities `weights`,
    and an initial seed set `seed_set`, runs `num_simulations` simulations
    and returns the average number of activated nodes.
"""
function simulate_ic(
    g::SimpleDiGraph, 
    weights::Dict{Tuple{Int, Int}, Float64}, 
    seed_set::Vector{Int}, 
    n_simulations::Int = 10_000
)

    total_activated = 0
    
    for _ ∈ 1:n_simulations
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

        total_activated += length(activated)
    end

    return total_activated / n_simulations
end


"""
    Influence Maximization in Directed Graphs in the Independent Cascade Model.
    Given a directed graph `g`, edge influence probabilities `weights`,
    and a budget `k`, selects `k` nodes to maximize the expected spread of influence.
    Returns a ≈ (1 - 1/e) approximate solution using a greedy algorithm.
"""
function influence_maximization_ic(
    g::SimpleDiGraph,
    weights::Dict{Tuple{Int, Int}, Float64},
    k::Int,
    n_simulations_small::Int = 1_000,   
    n_simulations::Int = 10_000
)

    n = nv(g)
    solution = Int[]
    current_spread = 0.0
    remaining = Set(1:n)

    for _ ∈ 1:min(k, n)
        isempty(remaining) && break

        u = -1
        Δ = -Inf

        for v ∈ remaining
            push!(solution, v)

            gain_small_sim = simulate_ic(g, weights, solution, n_simulations_small) - current_spread
            if gain_small_sim > Δ
                gain_large_sim = simulate_ic(g, weights, solution, n_simulations) - current_spread
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

    final_spread = simulate_ic(g, weights, solution, n_simulations)
    return solution, final_spread
end


precompile(
    simulate_ic, 
    (SimpleDiGraph, Dict{Tuple{Int, Int}, Float64}, Vector{Int}, Int)
)

precompile(
    influence_maximization_ic, 
    (SimpleDiGraph, Dict{Tuple{Int, Int}, Float64}, Int, Int, Int)
)
