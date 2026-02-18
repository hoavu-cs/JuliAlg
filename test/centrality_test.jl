using Test
using Graphs
using JuliAlg

# Ground truth values precomputed using naive BFS (unweighted) / Dijkstra (weighted)
# that enumerate all shortest paths and count intermediary vertices per (s,t) pair.

@testset "Betweenness Centrality Tests" begin

    @testset "Empty and trivial graphs" begin
        @testset "Empty graph" begin
            g = SimpleGraph(0)
            bc = bw_centrality(g)
            @test isempty(bc)
        end

        @testset "Single vertex" begin
            g = SimpleGraph(1)
            bc = bw_centrality(g)
            @test bc == [0.0]
        end

        @testset "Two vertices, one edge" begin
            g = SimpleGraph(2)
            add_edge!(g, 1, 2)
            bc = bw_centrality(g, nothing, false)
            @test bc ≈ [0.0, 0.0] atol=1e-10
        end
    end

    @testset "Unweighted undirected graphs" begin
        @testset "Path graph P5" begin
            g = SimpleGraph(5)
            for i in 1:4
                add_edge!(g, i, i + 1)
            end
            bc = bw_centrality(g, nothing, false)
            @test bc ≈ [0.0, 3.0, 4.0, 3.0, 0.0] atol=1e-10
            @test argmax(bc) == 3
        end

        @testset "Star graph S5" begin
            g = SimpleGraph(5)
            for i in 2:5
                add_edge!(g, 1, i)
            end
            bc = bw_centrality(g, nothing, false)
            @test bc ≈ [6.0, 0.0, 0.0, 0.0, 0.0] atol=1e-10
            @test argmax(bc) == 1
        end

        @testset "Cycle graph C5" begin
            g = SimpleGraph(5)
            for i in 1:5
                add_edge!(g, i, mod1(i + 1, 5))
            end
            bc = bw_centrality(g, nothing, false)
            @test bc ≈ [1.0, 1.0, 1.0, 1.0, 1.0] atol=1e-10
            @test all(isapprox.(bc, bc[1]; atol=1e-10))
        end

        @testset "Complete graph K4" begin
            g = SimpleGraph(4)
            for i in 1:4, j in i+1:4
                add_edge!(g, i, j)
            end
            bc = bw_centrality(g, nothing, false)
            @test bc ≈ [0.0, 0.0, 0.0, 0.0] atol=1e-10
        end

        @testset "Diamond graph" begin
            g = SimpleGraph(4)
            add_edge!(g, 1, 2)
            add_edge!(g, 1, 3)
            add_edge!(g, 2, 4)
            add_edge!(g, 3, 4)
            bc = bw_centrality(g, nothing, false)
            @test bc ≈ [0.5, 0.5, 0.5, 0.5] atol=1e-10
        end

        @testset "Disconnected graph" begin
            g = SimpleGraph(4)
            add_edge!(g, 1, 2)
            add_edge!(g, 3, 4)
            bc = bw_centrality(g, nothing, false)
            @test bc ≈ [0.0, 0.0, 0.0, 0.0] atol=1e-10
        end
    end

    @testset "Unweighted directed graphs" begin
        @testset "Directed path" begin
            g = SimpleDiGraph(4)
            add_edge!(g, 1, 2)
            add_edge!(g, 2, 3)
            add_edge!(g, 3, 4)
            bc = bw_centrality(g, nothing, false)
            @test bc ≈ [0.0, 2.0, 2.0, 0.0] atol=1e-10
            @test bc[2] > 0
            @test bc[3] > 0
            @test bc[1] == 0.0
            @test bc[4] == 0.0
        end

        @testset "Directed cycle" begin
            g = SimpleDiGraph(4)
            add_edge!(g, 1, 2)
            add_edge!(g, 2, 3)
            add_edge!(g, 3, 4)
            add_edge!(g, 4, 1)
            bc = bw_centrality(g, nothing, false)
            @test bc ≈ [3.0, 3.0, 3.0, 3.0] atol=1e-10
        end

        @testset "Directed star outward" begin
            g = SimpleDiGraph(5)
            for i in 2:5
                add_edge!(g, 1, i)
            end
            bc = bw_centrality(g, nothing, false)
            @test bc ≈ [0.0, 0.0, 0.0, 0.0, 0.0] atol=1e-10
        end

        @testset "Directed graph with shortcut" begin
            g = SimpleDiGraph(4)
            add_edge!(g, 1, 2)
            add_edge!(g, 2, 3)
            add_edge!(g, 3, 4)
            add_edge!(g, 1, 4)
            bc = bw_centrality(g, nothing, false)
            @test bc ≈ [0.0, 1.0, 1.0, 0.0] atol=1e-10
        end
    end

    @testset "Weighted directed graphs" begin
        @testset "Weighted path where shortcut is longer" begin
            g = SimpleDiGraph(4)
            add_edge!(g, 1, 2)
            add_edge!(g, 2, 3)
            add_edge!(g, 3, 4)
            add_edge!(g, 1, 4)
            weights = Dict((1, 2) => 2.0, (2, 3) => 1.0, (3, 4) => 3.0, (1, 4) => 10.0)
            bc = bw_centrality(g, weights, false)
            @test bc ≈ [0.0, 2.0, 2.0, 0.0] atol=1e-10
            @test bc[2] > 0
            @test bc[3] > 0
        end

        @testset "Weighted path where shortcut is shorter" begin
            g = SimpleDiGraph(4)
            add_edge!(g, 1, 2)
            add_edge!(g, 2, 3)
            add_edge!(g, 3, 4)
            add_edge!(g, 1, 4)
            weights = Dict((1, 2) => 2.0, (2, 3) => 1.0, (3, 4) => 3.0, (1, 4) => 5.0)
            bc = bw_centrality(g, weights, false)
            @test bc ≈ [0.0, 1.0, 1.0, 0.0] atol=1e-10
        end

        @testset "Weighted triangle" begin
            g = SimpleDiGraph(3)
            add_edge!(g, 1, 2)
            add_edge!(g, 2, 3)
            add_edge!(g, 1, 3)
            weights = Dict((1, 2) => 1.0, (2, 3) => 1.0, (1, 3) => 5.0)
            bc = bw_centrality(g, weights, false)
            @test bc ≈ [0.0, 1.0, 0.0] atol=1e-10
            @test bc[2] > 0
        end

        @testset "Equal weights match unweighted" begin
            g = SimpleDiGraph(5)
            add_edge!(g, 1, 2)
            add_edge!(g, 2, 3)
            add_edge!(g, 3, 4)
            add_edge!(g, 4, 5)
            add_edge!(g, 1, 3)
            weights = Dict((1, 2) => 1.0, (2, 3) => 1.0, (3, 4) => 1.0, (4, 5) => 1.0, (1, 3) => 1.0)
            bc_w = bw_centrality(g, weights, false)
            bc_u = bw_centrality(g, nothing, false)
            expected = [0.0, 0.0, 4.0, 3.0, 0.0]
            @test bc_w ≈ expected atol=1e-10
            @test bc_u ≈ expected atol=1e-10
        end
    end

    @testset "Weighted undirected graphs" begin
        @testset "Weighted path graph" begin
            g = SimpleGraph(4)
            add_edge!(g, 1, 2)
            add_edge!(g, 2, 3)
            add_edge!(g, 3, 4)
            add_edge!(g, 1, 4)
            weights = Dict(
                (1, 2) => 1.0, (2, 1) => 1.0,
                (2, 3) => 2.0, (3, 2) => 2.0,
                (3, 4) => 5.0, (4, 3) => 5.0,
                (1, 4) => 7.0, (4, 1) => 7.0
            )
            bc = bw_centrality(g, weights, false)
            @test bc ≈ [0.0, 1.0, 1.0, 0.0] atol=1e-10
        end


        @testset "Weighted star" begin
            g = SimpleGraph(5)
            for i in 2:5
                add_edge!(g, 1, i)
            end
            weights = Dict{Tuple{Int,Int}, Float64}()
            for i in 2:5
                weights[(1, i)] = Float64(i)
                weights[(i, 1)] = Float64(i)
            end
            bc = bw_centrality(g, weights, false)
            @test bc ≈ [6.0, 0.0, 0.0, 0.0, 0.0] atol=1e-10
            @test argmax(bc) == 1
        end
    end

    @testset "Normalization" begin
        @testset "Normalized undirected" begin
            g = SimpleGraph(5)
            for i in 1:4
                add_edge!(g, i, i + 1)
            end
            bc_norm = bw_centrality(g)
            @test bc_norm ≈ [0.0, 0.5, 2/3, 0.5, 0.0] atol=1e-10
        end

        @testset "Normalized directed" begin
            g = SimpleDiGraph(4)
            add_edge!(g, 1, 2)
            add_edge!(g, 2, 3)
            add_edge!(g, 3, 4)
            bc_norm = bw_centrality(g)
            @test bc_norm ≈ [0.0, 1/3, 1/3, 0.0] atol=1e-10
        end
    end

    @testset "Random graphs" begin
        import Random

        @testset "Random undirected graph (20 vertices)" begin
            Random.seed!(42)
            g = SimpleGraph(20)
            for i in 1:20, j in i+1:20
                if rand() < 0.2
                    add_edge!(g, i, j)
                end
            end
            bc = bw_centrality(g, nothing, false)
            expected = [7.666666666666667, 8.952380952380953, 5.845238095238093,
                        4.309523809523809, 1.8928571428571428, 32.42857142857145,
                        9.833333333333332, 0.0, 0.0, 9.833333333333332,
                        3.6428571428571432, 36.95238095238096, 6.9523809523809526,
                        0.0, 17.642857142857135, 33.47619047619048,
                        20.749999999999996, 25.40476190476192, 5.416666666666666, 0.0]
            @test bc ≈ expected atol=1e-8
        end

        @testset "Random directed graph (20 vertices)" begin
            Random.seed!(123)
            g = SimpleDiGraph(20)
            for i in 1:20, j in 1:20
                if i != j && rand() < 0.15
                    add_edge!(g, i, j)
                end
            end
            bc = bw_centrality(g, nothing, false)
            expected = [11.199999999999998, 72.16666666666677, 16.583333333333336,
                        4.366666666666667, 10.733333333333333, 0.3333333333333333,
                        72.93333333333342, 59.38333333333338, 28.333333333333332,
                        8.899999999999999, 17.733333333333334, 28.866666666666667,
                        36.86666666666668, 24.41666666666665, 36.9, 1.5,
                        39.23333333333333, 0.0, 16.349999999999998, 33.199999999999996]
            @test bc ≈ expected atol=1e-8
        end

        @testset "Random weighted directed graph (15 vertices)" begin
            Random.seed!(77)
            g = SimpleDiGraph(15)
            weights = Dict{Tuple{Int, Int}, Float64}()
            for i in 1:15, j in 1:15
                if i != j && rand() < 0.2
                    add_edge!(g, i, j)
                    weights[(i, j)] = rand() * 10.0
                end
            end
            bc = bw_centrality(g, weights, false)
            expected = [76.0, 48.0, 12.0, 23.0, 60.0, 0.0, 25.0, 16.0, 4.0,
                        39.0, 18.0, 13.0, 0.0, 21.0, 0.0]
            @test bc ≈ expected atol=1e-8
        end

        @testset "Random weighted undirected graph (15 vertices)" begin
            Random.seed!(99)
            g = SimpleGraph(15)
            weights = Dict{Tuple{Int, Int}, Float64}()
            for i in 1:15, j in i+1:15
                if rand() < 0.25
                    add_edge!(g, i, j)
                    w = rand() * 10.0
                    weights[(i, j)] = w
                    weights[(j, i)] = w
                end
            end
            bc = bw_centrality(g, weights, false)
            expected = [0.0, 1.0, 27.0, 0.0, 1.0, 34.0, 13.0, 15.0, 40.0,
                        9.0, 7.0, 21.0, 0.0, 34.0, 1.0]
            @test bc ≈ expected atol=1e-8
        end

        @testset "Large graph for multithreading benchmark" begin
            Random.seed!(2024)
            g = SimpleGraph(10000)

            for _ in 1:20000
                u, v = rand(1:5000, 2)
                if u != v
                    add_edge!(g, u, v)
                end
            end

            for _ in 1:15000
                u, v = rand(5001:10000, 2)
                if u != v
                    add_edge!(g, u, v)
                end
            end

            bc = bw_centrality(g, nothing, false)
            @test length(bc) == 10000
        end
    end
end