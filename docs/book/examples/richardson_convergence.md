# Example: richardson_convergence

**Crate**: `cfd-validation`
**Run**: `cargo run -p cfd-validation --example richardson_convergence`
**Source**: [`crates/cfd-validation/examples/richardson_convergence.rs`](../../../crates/cfd-validation/examples/richardson_convergence.rs)

## What This Example Demonstrates

Richardson extrapolation on the lid-driven cavity benchmark at `Re = 100` across
three grid levels (16, 32, 64). Estimates the observed order of accuracy `p`,
the extrapolated exact solution, and the Grid Convergence Index (GCI).

| Quantity | Definition |
|---|---|
| Observed order `p` | `ln(|(f3 − f2) / (f2 − f1)|) / ln(r)` |
| Extrapolated exact | `f_rich = f1 + (f1 − f2) / (r^p − 1)` |
| GCI | `1.25 · |f1 − f2| / (f1 · (r^p − 1))` |
| Refinement ratio `r` | `2.0` (uniform halving) |

## Key Code Snippet

```rust
use cfd_core::error::Result;
use cfd_validation::benchmarks::{Benchmark, BenchmarkConfig, LidDrivenCavity};

// Re = 100 lid-driven cavity on three uniformly-refined grids.
let grids = [64, 32, 16];           // h1 < h2 < h3 ordering
let f = |n: usize| -> Result<f64> {
    let cavity = LidDrivenCavity::new(1.0, 1.0, 100.0);
    let cfg = BenchmarkConfig::for_grid(n);
    let result = cavity.run(&cfg)?;
    Ok(result.values[0])
};

let (f1, f2, f3) = (f(64)?, f(32)?, f(16)?;
let r = 2.0;
let p = ((f3 - f2) / (f2 - f1)).abs().ln() / r.ln();
let f_rich = f1 + (f1 - f2) / (r.powf(p) - 1.0);
let gci = 1.25 * (f1 - f2).abs() / f1 / (r.powf(p) - 1.0);
```

## Physics Background

Richardson extrapolation (Richardson & Glauert, 1927; later formalized by Roache
as the Grid Convergence Index) combines solutions from three uniformly-refined
grids to estimate both the *observed* order of accuracy and the continuum value
the finite-difference solution is converging toward. The GCI with the safety
factor 1.25 (Roache, 1994) provides a conservative asymptotic bound on the
discretization error; for a formally second-order scheme on a smooth benchmark
such as lid-driven cavity at `Re = 100`, the observed `p` should land near `2`,
and a small GCI certifies that the resolved grids lie in the asymptotic regime.
This example serves as CFDrs's code-verification oracle: a dip in `p` flags
implementation drift.

## Book Chapter

[← Validation Suite](../crate_validation.md)
