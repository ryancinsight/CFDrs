# Example: Venturi Blood Flow Validation

**Crate**: `cfd-suite` (workspace root)
**Run**: `cargo run --example venturi_blood_flow_validation`
**Source**: [`examples/venturi_blood_flow_validation.rs`](../../examples/venturi_blood_flow_validation.rs)

## What This Example Demonstrates

Bernoulli + viscous-loss Venturi model validated against Casson and Carreau-Yasuda blood rheology.

| API | Purpose |
|---|---|
| `cfd_2d::solvers::venturi_flow::{BernoulliVenturi, ViscousVenturi, VenturiGeometry}` | Inlet and throat shear rates drive rheology |
| `CassonBlood` | Inlet and throat shear rates drive rheology |
| `CarreauYasudaBlood` | Inlet and throat shear rates drive rheology |

## Physics Background

Inlet and throat shear rates drive rheology. Continuity validated to 8-epsilon tolerance.

## Book Chapter

[← Part IV — Biomedical and Specialized Flows](../biomedical_flows.md)