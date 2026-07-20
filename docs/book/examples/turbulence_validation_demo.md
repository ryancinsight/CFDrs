# Example: turbulence_validation_demo

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example turbulence_validation_demo`  
**Source**: [`examples/turbulence_validation_demo.rs`](../../examples/turbulence_validation_demo.rs)

## What This Example Demonstrates

Runs the comprehensive turbulence validation suite for k-ε, k-ω SST, and
Spalart-Allmaras via `run_turbulence_validation::<f64>()`.

| Validation | Method | Target |
|---|---|---|
| Homogeneous turbulence decay | Analytical k-ε decay | < 1 % error |
| Wall boundary conditions | Near-wall k, ω, νt | Match literature |
| Numerical stability | NaN/inf check over iteration | None |
| Convergence | Monotonic physical evolution | Pass |

## Key Code Snippet

```rust
use cfd_2d::physics::turbulence::run_turbulence_validation;

// Runs all four validation passes in sequence
run_turbulence_validation::<f64>();
```

## Physics Background

**Homogeneous turbulence decay** follows `k(t) = k₀ (1 + t/t₀)^(−n)` where n
depends on the turbulence model constants. The k-ε model with standard
constants should reproduce n ≈ 1.2 – 1.4 over the initial decay period.

**Near-wall boundary conditions**: at the wall, k = 0 (no-slip),
`ω_wall = 6ν / (β₁ Δy²)` (k-ω SST), and `νt_wall = 0`.

## Book Chapter

[← Turbulence Models and Cavitation](../turbulence_multiphase.md)

