using Test
using JuliAlg

# Brute force makespan - tries all possible assignments of jobs to machines
function brute_force_makespan(jobs::Vector{Float64}, m::Int)
    n = length(jobs)
    best_makespan = Inf
    
    # Iterate through all m^n possible assignments
    num_assignments = m^n
    for assignment in 0:(num_assignments - 1)
        loads = zeros(Float64, m)
        temp = assignment
        for job_idx in 1:n
            machine = (temp % m) + 1
            temp = div(temp, m)
            loads[machine] += jobs[job_idx]
        end
        makespan = maximum(loads)
        if makespan < best_makespan
            best_makespan = makespan
        end
    end
    
    return best_makespan
end

# Generate random test cases within brute force limit
function generate_brute_force_test_cases(max_brute_force::Int=2000000)
    test_cases = []
    
    # Different (n, m) combinations where m^n <= max_brute_force
    # n = number of jobs, m = number of machines
    configs = [
        (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), 
        (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7),
        (4, 2), (4, 3), (4, 4), (4, 5), (4, 6), (4, 7),
        (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7),
        (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 7),
        (7, 2), (7, 3), (7, 4), (7, 5), (7, 6), (7, 7)
    ]
    
    for (n, m) in configs
        # Check if m^n <= max_brute_force
        if m^n <= max_brute_force
            # Generate multiple random test cases for each config
            for _ in 1:10
                # Random job sizes between 1 and 20
                jobs = Float64[rand(1:20) for _ in 1:n]
                push!(test_cases, (jobs, m))
            end
        end
    end
    
    return test_cases
end


@testset "lpt_makespan approximation ratio" begin
    # Generate test cases within brute force limit (m^n <= 2000000)
    test_cases = generate_brute_force_test_cases(2000000)
    
    @testset "brute force verification for $(length(test_cases)) test cases" begin
        for (jobs, m) in test_cases
            # Compute LPT makespan
            lpt_result, _ = lpt_makespan(jobs, m)
            
            # Compute optimal (ground truth) via brute force
            optimal = brute_force_makespan(jobs, m)
            
            # Verify approximation ratio: makespan <= 4/3 * optimal
            @test lpt_result <= (4/3) * optimal + 1e-6
        end
    end
end

