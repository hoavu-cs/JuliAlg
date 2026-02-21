
using DataStructures

"""
    lpt_makespan(jobs::Vector{Float64}, m::Int)

LPT (Longest Processing Time) heuristic for the makespan problem.
Approximation ratio: `4/3 - 1/(3m)` (Graham, 1966).

Returns a tuple `(makespan, assignments)` where:
- `makespan`: maximum load across all machines
- `assignments`: vector of machine assignments for each job
"""
function lpt_makespan(jobs::Vector{Float64}, m::Int)
    sorted_jobs = sortperm(jobs, rev=true) # sort jobs in decreasing order
    loads = PriorityQueue{Int, Float64}()
    assignments = zeros(Int, length(jobs)) # to store machine assignments for each job

    for i in 1:m
        enqueue!(loads, i, 0.0)  
    end

    for idx in sorted_jobs
        machine, load = peek(loads)
        dequeue!(loads)
        new_load = load + jobs[idx]
        enqueue!(loads, machine, new_load)
        assignments[idx] = machine
    end
    
    return maximum(values(loads)), assignments 
end

precompile(lpt_makespan, (Vector{Float64}, Int))
