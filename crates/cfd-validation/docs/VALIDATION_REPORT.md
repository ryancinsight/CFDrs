# cfd-validation Validation Report

## Purpose
- Aggregates manufactured solutions, conservation checks, convergence studies, and benchmarking for the project.

## Components
- Manufactured solutions (diffusion/advection/advection–diffusion/Navier–Stokes, turbulent extensions).
- Convergence (Richardson), error metrics (norms/statistics), conservation (mass/momentum/energy), benchmarking suites.

## Tests
- Present modules include analytical benchmarks, error metric tests, conservation traits; extend harness to collect per-crate summaries.

## Invariants
- Accuracy thresholds per physics; conservation tolerances; convergence rate expectations.

## Assumptions
- Consistent units and boundary semantics across dependent crates; reproducible seeds for stochastic components.

## Reporting
- Markdown/JSON/HTML reporters to summarize coverage and results; integrate per-crate reports.
