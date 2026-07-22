# Example: Serpentine-Channel Mixing Analysis

**Crate**: `cfd-suite` (workspace root)
**Run**: `cargo run --example serpentine_mixing_comprehensive`
**Source**: [`examples/serpentine_mixing_comprehensive.rs`](../../../examples/serpentine_mixing_comprehensive.rs)

## What This Example Demonstrates

Analytical transverse-diffusion model + discretized flow/scalar-transport for a serpentine micromixer.

| API | Purpose |
|---|---|
| `cfd_2d::solvers::serpentine_flow::{SerpentineGeometry, SerpentineSolver2D, AdvectionDiffusionMixing}` | Mixing efficiency by Peclet number |

## Physics Background

Mixing efficiency by Peclet number. Discretized convection-diffusion vs analytical laminar mixing model.

## Book Chapter

[← Part IV — Biomedical and Specialized Flows](../biomedical_flows.md)