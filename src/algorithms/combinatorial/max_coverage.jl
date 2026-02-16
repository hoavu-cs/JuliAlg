"""
    max_coverage(subsets::Vector{Vector{Int}}, k::Int)

    Greedy approximation for the maximum coverage problem.
    Selects at most `k` subsets to maximize the number of covered elements.
    Returns the number of covered elements and the indices of the selected subsets.
    Approximation guarantee: `(1 - 1/e)` ≈ 0.632.
"""
function max_coverage(subsets::Vector{Vector{Int}}, k::Int)
    m = length(subsets)
    k = min(k, m)

    subset_sets = [Set(s) for s in subsets]
    covered = Set{Int}()
    selected = Int[]
    used = falses(m)

    for _ in 1:k
        gains = zeros(m)

        Threads.@threads for i in 1:m
            used[i] && continue
            gain = length(setdiff(subset_sets[i], covered))
            gains[i] = gain
        end

        best_idx = argmax(gains)
        gains[best_idx] == 0 && break

        push!(selected, best_idx)
        used[best_idx] = true
        union!(covered, subset_sets[best_idx])
    end

    return length(covered), sort!(selected)
end

precompile(max_coverage, (Vector{Vector{Int}}, Int))