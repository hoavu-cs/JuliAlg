# JuliOpt

A Julia package for combinatorial optimization and graph algorithms. Implements exact, approximate, and heuristic solutions for classic problems in optimization and network analysis.

Claude is often used to generate test cases and some documentation, but the core algorithms and implementations are mostly written and optimized by the author. 

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/<owner>/JuliOpt.jl")
```

Or for local development:

```bash
git clone <repo-url>
cd JuliOpt
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
| `influence_maximization_ic(g, weights, k)` | Influence Maximization | Greedy + Monte Carlo IC | (1 - 1/e)-approx |
| `simulate_ic(g, weights, seed_set)` | IC Spread Estimation | Monte Carlo simulation | - |
| `densest_subgraph(G)` | Densest Subgraph | Goldberg's algorithm (binary search + max-flow) | Exact |

## Usage

All algorithms return a tuple of `(objective_value, selected_items)`.

### Knapsack

```julia
using JuliOpt

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
using Graphs, JuliOpt

g = SimpleDiGraph(5)
add_edge!(g, 1, 2); add_edge!(g, 2, 3)
add_edge!(g, 3, 4); add_edge!(g, 4, 5)

weights = Dict((src(e), dst(e)) => 0.5 for e in edges(g))
k = 2

seeds, spread = influence_maximization_ic(g, weights, k)
```

### Densest Subgraph

```julia
using Graphs, JuliOpt

g = complete_graph(5)
add_vertex!(g)
add_edge!(g, 1, 6)  # pendant vertex

S, lambda, density = densest_subgraph(g)
# S = [1, 2, 3, 4, 5], density = 2.0
```

## Testing

```bash
# Run all tests
julia --project -e 'using Pkg; Pkg.test()'

# Run a single test file
julia --project -e 'using JuliOpt; using Test; include("test/knapsack_test.jl")'

# With threads (needed for influence maximization)
julia --threads=auto --project -e 'using Pkg; Pkg.test()'
```

## Dependencies

Graphs, GraphsFlows, DataStructures, OffsetArrays, Combinatorics. See `Project.toml` for version constraints. Requires Julia >= 1.10.
