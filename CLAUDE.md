# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

JuliAlg is a Julia package for combinatorial optimization and graph algorithms. It implements exact and approximate solutions for knapsack, bin packing, weighted interval scheduling, and influence maximization problems.

## Common Commands

```bash
# Run all tests
julia --project -e 'using Pkg; Pkg.test()'

# Run a single test file (e.g., knapsack only)
julia --project -e 'using JuliAlg; using Test; include("test/knapsack_test.jl")'

# Run tests with multiple threads (needed for influence maximization)
julia --threads=auto --project -e 'using Pkg; Pkg.test()'

# Install/update dependencies
julia --project -e 'using Pkg; Pkg.instantiate()'

# Interactive development with Revise
julia --project -e 'using Revise; using JuliAlg'
```

## Architecture

All source code lives in `src/`, organized by problem domain:

- **`src/JuliAlg.jl`** — Module entry point. Exports 6 public functions and includes all algorithm files.
- **`src/algorithms/combinatorial/`** — Knapsack (exact DP + PTAS), bin packing (Best-Fit Decreasing), weighted interval scheduling (DP + binary search).
- **`src/algorithms/graphs/`** — Influence maximization using the Independent Cascade model with greedy seed selection and Monte Carlo simulation. Densest subgraph (Goldberg's exact algorithm) and densest at-most-k-subgraph (degree-based pruning + brute force).

Key design decisions:
- All algorithms return a tuple of `(objective_value, selected_items)`.
- `simulate_ic` uses `Threads.@threads` for parallel Monte Carlo simulations.
- `ptas_knapsack` calls `exact_knapsack` internally as a subroutine with value-scaled inputs.
- OffsetArrays are used for 0-indexed DP tables in knapsack. SortedDict is used for efficient bin selection in bin packing.
- All exported functions have `precompile` directives at the end of `JuliAlg.jl`'s included files.

## Testing

Tests are in `test/` with one file per algorithm area. `test/runtests.jl` includes all test files inside a top-level `@testset`. Influence maximization tests use `Combinatorics.jl` for exhaustive brute-force validation on small graphs, including Erdos-Renyi random graphs.

## Dependencies

Core: `Graphs`, `DataStructures`, `OffsetArrays`, `Combinatorics`, `Revise`. See `Project.toml` for version constraints.