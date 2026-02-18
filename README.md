# JuliAlg

A Julia package for combinatorial optimization and graph algorithms. It contains implementation of exact, approximate, and heuristic solutions for classic problems in optimization and network analysis. I aim to implement every  algorithms that I can get a good grasp of that are 
- theoretically sound
- practical with applications in mind
- those that require more understanding or optimization. 

I also utilize multi-threading for algorithms that can benefit from parallelism. See the benchmarks section for some results. There are some NP-Hard problems with no known polynomial-time approximation. For these, I try to come up with heuristics that help reduce the search space (such as in densest at-most-k-subgraph). Some are dealt with using parameterization.

Claude Code / DeepSeek + Cline is often used to generate test cases and some documentation. The core algorithms and implementations are mostly written and optimized by me. 

The goal is to use this package internally for my other projects, but I would also be happy if it can be useful to others. Pull requests and contributions are welcome. 


## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/<owner>/JuliAlg.jl")
```

Or for local development:

```bash
git clone <repo-url>
cd JuliAlg
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Algorithms

### Combinatorial Optimization

| Function | Problem | Method | Guarantee |
|---|---|---|---|
| `exact_knapsack(W, weights, values)` | 0/1 Knapsack | Value-based DP | Exact, O(n * sum(values)) |
| `ptas_knapsack(W, epsilon, weights, values)` | 0/1 Knapsack | Value scaling + DP | (1 - epsilon)-approx |
| `bin_packing(items, bin_capacity)` | Bin Packing | Best-Fit Decreasing | <= (11/9)OPT + 6/9 |
| `weighted_interval_scheduling(starts, ends, weights)` | Weighted Interval Scheduling | DP + binary search | Exact, O(n log n) |
| `set_cover(subsets, costs)` | Weighted Set Cover | Greedy | O(ln n)-approx |
| `max_coverage(subsets, k)` | Maximum Coverage | Greedy | (1 - 1/e)-approx |

### Graph Algorithms

| Function | Problem | Method | Guarantee |
|---|---|---|---|
| `pagerank(G, weights, α)` | PageRank | Power iteration | Exact (up to convergence) |
| `influence_maximization_ic(g, weights, k)` | Influence Maximization | Greedy + Monte Carlo IC | (1 - 1/e)-approx |
| `simulate_ic(g, weights, seed_set)` | IC Spread Estimation | Monte Carlo simulation | - |
| `densest_subgraph(G)` | Densest Subgraph | Goldberg's algorithm (binary search + max-flow) | Exact |
| `densest_subgraph_peeling(G)` | Densest Subgraph | Charikar's peeling algorithm | 1/2-approx |
| `densest_at_most_k_subgraph(G, k)` | Densest At-Most-k Subgraph | Degree-based pruning + brute force | Heuristic |
| `k_core_decomposition(G)` | K-Core Decomposition | Iterative peeling | Exact, O(m) |
| `bw_centrality(G, weights, normalized)` | Betweenness Centrality | Brandes' algorithm (BFS / Dijkstra) | Exact, O(nm) / O(nm + n² log n) |

`pagerank`, `bw_centrality`, `influence_maximization_ic`, and `simulate_ic` accept both directed (`SimpleDiGraph`) and undirected (`SimpleGraph`) graphs. For weighted undirected graphs, ensure `weights[(u,v)] == weights[(v,u)]` for all edges.

## Usage

All algorithms return a tuple of `(objective_value, selected_items)`.

### Knapsack

```julia
using JuliAlg

W = 10
weights = [2, 3, 4, 5]
values  = [3, 4, 5, 6]

# Exact solution
value, items = exact_knapsack(W, weights, values)

# PTAS with epsilon = 0.1
value, items = ptas_knapsack(W, 0.1, weights, values)
```

### Bin Packing

```julia
items = [7, 5, 3, 4, 2, 6]
capacity = 10

num_bins, bins = bin_packing(items, capacity)
```

### Weighted Interval Scheduling

```julia
starts  = [1, 3, 0, 5, 8]
ends    = [4, 5, 6, 7, 9]
weights = [3, 2, 4, 7, 2]

value, selected = weighted_interval_scheduling(starts, ends, weights)
```

### Set Cover

```julia
subsets = [[1, 2, 3], [2, 4], [3, 4, 5]]
costs   = [1.0, 2.0, 1.5]

cost, selected = set_cover(subsets, costs)
```

### Maximum Coverage

```julia
subsets = [[1, 2], [2, 3, 4], [4, 5]]
k = 2

covered, selected = max_coverage(subsets, Int64(k))
```

### Influence Maximization

```julia
using Graphs, JuliAlg

g = SimpleDiGraph(5)
add_edge!(g, 1, 2); add_edge!(g, 2, 3)
add_edge!(g, 3, 4); add_edge!(g, 4, 5)

weights = Dict((src(e), dst(e)) => 0.5 for e in edges(g))
k = 2

seeds, spread = influence_maximization_ic(g, weights, k)
```

### PageRank

```julia
using Graphs, JuliAlg

g = SimpleDiGraph(4)
add_edge!(g, 1, 2); add_edge!(g, 2, 3); add_edge!(g, 3, 1); add_edge!(g, 2, 4)

r = pagerank(g)          # unweighted, α=0.85
r = pagerank(g, nothing, 0.9)  # custom damping factor

# Weighted
weights = Dict((1,2) => 2.0, (2,3) => 1.0, (3,1) => 1.0, (2,4) => 3.0)
r = pagerank(g, weights)
```

### Betweenness Centrality

```julia
using Graphs, JuliAlg

g = SimpleGraph(5)
for i in 1:4; add_edge!(g, i, i+1); end

bc = bw_centrality(g)           # normalized
bc = bw_centrality(g, nothing, false)  # unnormalized

# Weighted (symmetric weights required for undirected)
weights = Dict((1,2)=>1.0, (2,1)=>1.0, (2,3)=>2.0, (3,2)=>2.0,
               (3,4)=>1.0, (4,3)=>1.0, (4,5)=>3.0, (5,4)=>3.0)
bc = bw_centrality(g, weights)
```

### Densest Subgraph

```julia
using Graphs, JuliAlg

g = complete_graph(5)
add_vertex!(g)
add_edge!(g, 1, 6)  # pendant vertex

S, density = densest_subgraph(g)
# S = [1, 2, 3, 4, 5], density = 2.0

# Densest subgraph with at most k vertices
S, d = JuliAlg.densest_at_most_k_subgraph(g, 3)
# Finds the 3-vertex subset with highest density
```

## Testing

```bash
# Run all tests
julia --project -e 'using Pkg; Pkg.test()'

# Run a single test file
julia --project -e 'using JuliAlg; using Test; include("test/knapsack_test.jl")'

# With threads (needed for influence maximization)
julia --threads=auto --project -e 'using Pkg; Pkg.test()'
```

## Benchmarks for Thread Scaling

### Influence Maximization

| Threads | Median (ms) | Min (ms) | Speedup |
|---------|-------------|----------|---------|
| 1       | 7,111       | 7,095    | 1.0x    |
| 2       | 3,788       | 3,725    | 1.9x    |
| 4       | 2,006       | 1,994    | 3.5x    |
| 8       | 1,695       | 1,578    | 4.2x    |

### PageRank

| Threads | Median (ms) | Min (ms) | Speedup |
|---------|-------------|----------|---------|
| 1       | 677.84      | 517.57   | 1.0x    |
| 2       | 586.70      | 483.27   | 1.2x    |
| 4       | 562.86      | 491.88   | 1.2x    |
| 8       | 535.50      | 475.12   | 1.3x    |

### Set Cover

| Threads | Median (ms) | Min (ms) | Speedup |
|---------|-------------|----------|---------|
| 1       | 5,854       | 5,763    | 1.0x    |
| 2       | 4,141       | 4,040    | 1.4x    |
| 4       | 3,191       | 3,132    | 1.8x    |
| 8       | 2,892       | 2,834    | 2.0x    |

### Max Coverage

| Threads | Median (ms) | Min (ms) | Speedup |
|---------|-------------|----------|---------|
| 1       | 4,240       | 3,963    | 1.0x    |
| 2       | 2,945       | 2,889    | 1.4x    |
| 4       | 2,284       | 2,230    | 1.9x    |
| 8       | 2,040       | 2,021    | 2.1x    |

### Betweenness Centrality

| Threads | Median (ms) | Min (ms) | Speedup |
|---------|-------------|----------|---------|
| 1       | 15,739      | 14,430   | 1.0x    |
| 2       | 7,571       | 7,288    | 2.0x    |
| 4       | 3,675       | 3,454    | 4.2x    |
| 8       | 3,745       | 3,300    | 4.4x    |

```bash
# Run benchmarks
julia --project benchmarks/influence_maximization_bench.jl
julia --project benchmarks/pagerank_bench.jl
julia --project benchmarks/set_cover_bench.jl
julia --project benchmarks/max_coverage_bench.jl
julia --project benchmarks/centrality_bench.jl
```

## Dependencies

Graphs, GraphsFlows, DataStructures, OffsetArrays, Combinatorics. See `Project.toml` for version constraints. Requires Julia >= 1.10.
