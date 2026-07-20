# Example: turbulent_channel_flow

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example turbulent_channel_flow`  
**Source**: [`examples/turbulent_channel_flow.rs`](../../examples/turbulent_channel_flow.rs)

## What This Example Demonstrates

Fully-developed turbulent channel flow using the k-ω SST turbulence closure,
validated against DNS data of Moser, Kim & Mansour (1999) at Re_τ = 395.

| Concept | API |
|---|---|
| k-ω SST model creation | `KOmegaSSTModel::new(nx, ny)` |
| Near-wall k, ω initialization | Manual log-law initialisation |
| Model name introspection | `sst_model.name()` |

## Key Code Snippet

```rust
use cfd_2d::physics::turbulence::{KOmegaSSTModel, TurbulenceModel};

// Re_τ = u_τ · h / ν = 395  (DNS benchmark)
let re_tau = 395.0;
let h = 1.0; // half-height [m]
let nu = 1.5e-5; // kinematic viscosity [m²/s]
let u_tau = re_tau * nu / h;

let sst_model = KOmegaSSTModel::new(10, 100);
println!("Model: {}", sst_model.name());
```

## Physics Background

The **k-ω SST** (Shear Stress Transport) model blends the k-ω formulation in
the near-wall region with the k-ε model in the free stream:

```
Dk/Dt  = P_k − β*kω + ∇·[(ν + σ_k νt) ∇k]
Dω/Dt  = γ P_k/νt − β' ω² + ∇·[(ν + σ_ω νt) ∇ω] + CDkω
```

The cross-diffusion term CDkω activates in the transition zone.  
Reference: Menter (1994), AIAA J. 32(8), 1598–1605.

## Book Chapter

[← Turbulence Models and Cavitation](../turbulence_multiphase.md)

