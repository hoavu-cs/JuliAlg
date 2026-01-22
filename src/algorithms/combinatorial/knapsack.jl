using OffsetArrays

"""
    exact_knapsack(W::Int64, weights::Vector{Int64}, values::Vector{Int64})

    Exact value-based 0/1 knapsack.
    Returns the maximum achievable value within capacity `W`
    and the indices of the selected items.
    Time/space complexity: `O(n * sum(values))`.
"""
function exact_knapsack(W::Int64, weights::Vector{Int64}, values::Vector{Int64})
    n = length(weights)
    V = sum(values)
    INF = div(typemax(Int64), 2)
    
    dp = OffsetArray(fill(INF, n + 1, V + 1), 0:n, 0:V)
    dp[0, 0] = 0
    
    for i in 1:n
        wi = weights[i]
        vi = values[i]
        for v in 0:V
            dp[i, v] = dp[i - 1, v]  # not take
            if vi <= v
                dp[i, v] = min(dp[i, v], dp[i - 1, v - vi] + wi)
            end
        end
    end
    
    best_v = findlast(v -> dp[n, v] <= W, 0:V)
    best_v = isnothing(best_v) ? 0 : best_v - 1  # adjust for 0-indexing
    
    # Backtrack
    items = Int64[]
    i, v = n, best_v
    while i > 0 && v > 0
        wi, vi = weights[i], values[i]
        if vi <= v && dp[i, v] == dp[i - 1, v - vi] + wi
            push!(items, i)
            v -= vi
        end
        i -= 1
    end
    
    reverse!(items)
    return best_v, items
end

"""
    ptas_knapsack(W::Int64, epsilon::Float64,
                  weights::AbstractVector{Int64},
                  values::AbstractVector{Int64})

    Value-scaling PTAS for 0/1 knapsack.
    Returns a `(1 - ε)`-approximate solution value
    and the indices of the selected items.
"""
function ptas_knapsack(W::Int64, epsilon::Float64, weights::AbstractVector{Int64}, values::AbstractVector{Int64})
    n = length(weights)
    if n == 0
        return 0, Int64[]
    end

    max_value = maximum(values)
    K = (epsilon * max_value) / n # Scaling factor

    if !(K > 0.0)  
        return exact_knapsack(W, weights, values)
    end

    scaled_values = [floor(Int64, v / K) for v in values]

    # if scaling collapsed everything to zero, fallback
    if sum(scaled_values) == 0
        return exact_knapsack(W, weights, values)
    end

    opt_scaled, items = exact_knapsack(W, weights, scaled_values)
    actual_value = sum(values[i] for i in items)
    return actual_value, items
end

# precompile(exact_knapsack, (Int64, Vector{Int64}, Vector{Int64}))
# precompile(ptas_knapsack, (Int64, Float64, Vector{Int64}, Vector{Int64}))