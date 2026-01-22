using OffsetArrays

function weighted_interval_scheduling(start_times::Vector{Int64}, end_times::Vector{Int64}, weights::Vector{Int64})
    jobs = [(i, start_times[i], end_times[i], weights[i]) for i in 1:length(start_times)]
    n = length(jobs)
    sort!(jobs, by = x -> x[3])  
    
    # p[i] = rightmost job compatible with job i (0 if none)
    p = OffsetArray(zeros(Int64, n), 1:n)
    
    for i in 2:n
        s = 1
        e = i - 1
        pi = 0
        # binary search for the rightmost job that ends before jobs[i] starts
        while s <= e 
            m = div(s + e, 2)
            if jobs[m][3] <= jobs[i][2]
                pi = m
                s = m + 1
            else
                e = m - 1 
            end
        end
        p[i] = pi
    end
    
    # dp[i] = max weight using first i jobs (0-indexed: 0 to n)
    dp = OffsetArray(zeros(Int64, n + 1), 0:n)
    dp[0] = 0
    
    for i in 1:n
        dp[i] = max(dp[i-1], dp[p[i]] + jobs[i][4])
    end
    
    # Reconstruct solution
    S = Int64[]
    i = n
    while i >= 1 
        take = dp[p[i]] + jobs[i][4]
        skip = dp[i-1]
        if take >= skip
            push!(S, jobs[i][1])
            i = p[i]
        else
            i -= 1
        end
    end
    
    reverse!(S)
    return dp[n], S
end

# precompile(weighted_interval_scheduling, (Vector{Int64}, Vector{Int64}, Vector{Int64}))
