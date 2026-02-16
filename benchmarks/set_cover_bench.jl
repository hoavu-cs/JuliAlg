using JuliOpt

if !isempty(ARGS) && ARGS[1] == "--worker"
    println("Threads: ", Threads.nthreads())

    # Warmup
    set_cover([Int[1, 2], Int[2, 3]], Float64[1.0, 1.0])

    times = Float64[]
    for _ in 1:5
        t = @elapsed include(joinpath(@__DIR__, "..", "test", "set_cover_test.jl"))
        push!(times, t)
    end
    sort!(times)
    println("\nMedian: $(round(times[3] * 1000, digits=2)) ms  Min: $(round(times[1] * 1000, digits=2)) ms")
else
    println("=== set_cover: thread scaling benchmark ===\n")

    thread_counts = [1, 2, 4, 8]
    for nt in thread_counts
        println("--- Threads: $nt ---")
        run(`julia --project --threads=$nt $(@__FILE__) --worker`)
        println()
    end
end