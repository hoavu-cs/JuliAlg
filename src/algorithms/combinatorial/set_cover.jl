"""
    set_cover(universe::Vector{Int}, subsets::Vector{Vector{Int}}, costs::Vector{Float64})

    Greedy approximation for the weighted set cover problem.
    Selects subsets to cover every element in `universe` with minimum total cost.
    Returns the total cost and the indices of the selected subsets.
    Approximation guarantee: `O(ln(n))` where `n = |universe|`.
"""
function set_cover(subsets::Vector{Vector{Int}}, costs::Vector{Float64})
    m = length(subsets)
    
    universe = Set{Int}()
    for s in subsets
        union!(universe, s)
    end
    
    uncovered = copy(universe)
    subset_sets = [Set(s) for s in subsets]
    selected = Int[]
    total_cost = 0.0
    used = falses(m)

    while !isempty(uncovered) || count(used) == m
        ratios = zeros(m)

        Threads.@threads for i in 1:m
            used[i] && continue
            num_newly_covered = length(uncovered ∩ subset_sets[i])
            ratios[i] = num_newly_covered / costs[i]
        end

        best_idx = argmax(ratios)
        ratios[best_idx] <= 0.0 && break

        push!(selected, best_idx)
        used[best_idx] = true
        total_cost += costs[best_idx]
        uncovered = setdiff(uncovered, subset_sets[best_idx])
    end

    return total_cost, sort!(selected)
end

precompile(set_cover, (Vector{Vector{Int}}, Vector{Float64}))
