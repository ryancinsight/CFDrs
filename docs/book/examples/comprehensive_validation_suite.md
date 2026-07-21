# Example: comprehensive_validation_suite

**Crate**: `cfd-validation`
**Run**: `cargo run -p cfd-validation --example comprehensive_validation_suite`
**Source**: [`crates/cfd-validation/examples/comprehensive_validation_suite.rs`](../../crates/cfd-validation/examples/comprehensive_validation_suite.rs)

## What This Example Demonstrates

A single end-to-end runner that executes every implemented 2D and 3D benchmark
in `cfd-validation` and emits a unified validation report with convergence
status and per-benchmark metrics. Mathematical invariants are verified at the
end of the run.

| Benchmark | Dimension | Notes |
|---|---|---|
| LidDrivenCavity | 2D | Re = 100 canonical separator |
| BifurcationFlow | 2D | Symmetric + asymmetric Murray-law split |
| VenturiFlow | 2D | Bernoulli pressure-drop check |
| SerpentineFlow | 2D | Dean secondary flow |
| BifurcationFlow3D | 3D | Casson blood in Y-junction |
| VenturiFlow3D | 3D | Cavitation number check |
| SerpentineFlow3D | 3D | Dean number vs threshold |

## Key Code Snippet

```rust
use cfd_core::error::Result;
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_validation::benchmarks::{
    Benchmark, BenchmarkConfig, BenchmarkRunner,
    BifurcationFlow, BifurcationFlow3D, LidDrivenCavity,
    SerpentineFlow, SerpentineFlow3D, VenturiFlow, VenturiFlow3D,
};

// Heterogeneous benchmark set collected behind Box<dyn Benchmark<f64>>
// at the non-hot suite-entry boundary (per standards dyn-policy exception).
let benchmarks: Vec<Box<dyn Benchmark<f64>>> = vec![
    Box::new(LidDrivenCavity::new(1.0, 1.0, 100.0)),
    Box::new(BifurcationFlow::symmetric()),
    Box::new(BifurcationFlow3D::new(CassonBlood::normal_blood())),
    // ... VenturiFlow, SerpentineFlow, VenturiFlow3D, SerpentineFlow3D ...
];

let config = BenchmarkConfig { tolerance: 1e-4, ..Default::default() };
let results: Vec<_> = benchmarks.iter()
    .map(|b| b.run(&config))
    .collect::<Result<Vec<_>>>()?;
let report = BenchmarkRunner::generate_report(&results);

for r in &report.benchmarks {
    let converged = r.convergence.last().is_some_and(|&c| c < config.tolerance);
    println!("- {:<30} [{}]", r.name, if converged { "CONVERGED" } else { "NON-CONVERGED" });
}
```

## Physics Background

The suite spans canonical fluid-dynamics regimes — lubrication (lid-driven
cavity), branching-network haemodynamics (Murray- and Hagen-Poiseuille-driven
bifurcation split), inviscid-recoverable pressure drop (Bernoulli over Venturi
geometries), and curved-channel secondary flow (Dean vortices). Each benchmark's
convergence is judged against a tolerance on its last-iteration residual. The
closing invariant check at line 88 (`All mathematical invariants verified`)
asserts that the suite's physics laws — mass conservation, Murray's-law flow
split, and Bernoulli pressure-drop relations — hold across the combined runs.

## Book Chapter

[← Validation Suite](../crate_validation.md)
