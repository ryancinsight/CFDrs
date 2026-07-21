# Example: Validated Venturi Flow Benchmark

**Crate**: `cfd-suite` (workspace root)
**Run**: `cargo run --example venturi_validated`
**Source**: [`examples/venturi_validated.rs`](../../examples/venturi_validated.rs)

## What This Example Demonstrates

Full venturi flow validation across multiple operating points against published experimental data.

| API | Purpose |
|---|---|
| `cfd_1d::{evaluate_venturi_screening}` | Comparison of 1-D screening, 2-D CFD, and experimental pressure drop data |
| `cfd_2d::solvers::venturi_flow` | Comparison of 1-D screening, 2-D CFD, and experimental pressure drop data |
| `cfd_validation` | Comparison of 1-D screening, 2-D CFD, and experimental pressure drop data |

## Physics Background

Comparison of 1-D screening, 2-D CFD, and experimental pressure drop data. Discharge coefficient Cd validation.

## Book Chapter

[← Part VII — Performance and Atlas Integration](../performance_and_atlas.md)