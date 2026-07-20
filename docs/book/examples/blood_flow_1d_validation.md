# Example: blood_flow_1d_validation

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example blood_flow_1d_validation --no-default-features`  
**Source**: [`examples/blood_flow_1d_validation.rs`](../../examples/blood_flow_1d_validation.rs)

## What This Example Demonstrates

Four literature validation cases for 1D non-Newtonian blood flow, each exercising
a different rheology or network model:

| Case | Reference | What is validated |
|---|---|---|
| 1. Poiseuille flow — Casson model | Merrill et al. (1969) | Pressure drop vs analytical |
| 2. Symmetric bifurcation | Murray (1926) | Flow split obeys Murray's Law |
| 3. Asymmetric bifurcation | Hagen-Poiseuille | Pressure-balanced flow split |
| 4. Fåhræus-Lindqvist effect | Pries et al. (1992) | Viscosity reduction in small tubes |

## Key Code Snippet

```rust
use cfd_1d::domain::channel::Channel;

// Casson model: √τ = √τ_y + √(μ_∞ · γ̇)
// τ_y ≈ 0.0056 Pa (yield stress, Ht=45%)
// μ_∞ ≈ 0.00345 Pa·s (infinite-shear viscosity)

// Carreau-Yasuda for full shear-rate range:
// μ(γ̇) = μ_∞ + (μ_0 − μ_∞) · [1 + (λγ̇)^a]^((n-1)/a)
```

## Physics Background

Blood is a **shear-thinning** non-Newtonian fluid. At high shear rates
(arteries, cardiac cycle peak) it approaches Newtonian behaviour (μ ≈ 3.5 mPa·s).
At low shear rates (venules, stenosed vessels) it exhibits:

- **Yield stress** from rouleaux formation
- **Fåhræus-Lindqvist effect**: apparent viscosity drops for tube diameter < 300 μm

## Book Chapter

[← Blood Flow and Rheology Workflows](../biomedical_flows.md)

