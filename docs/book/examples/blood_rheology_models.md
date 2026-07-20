# Example: blood_rheology_models

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example blood_rheology_models`  
**Source**: [`examples/blood_rheology_models.rs`](../../examples/blood_rheology_models.rs)

## What This Example Demonstrates

Five blood-viscosity models compared over the full physiological shear-rate range
(0.01 – 1 000 s⁻¹), with predicted pressure drops in a 1 mm millifluidic channel.

| Model | Constitutive law | Reference |
|---|---|---|
| Newtonian | μ = const | — |
| Casson | √τ = √τ_y + √(μ_∞ γ̇) | Merrill et al. (1969) |
| Carreau-Yasuda | μ = μ_∞ + (μ₀−μ_∞)[1+(λγ̇)ᵃ]^((n-1)/a) | Cho & Kensey (1991) |
| Cross | μ = μ_∞ + (μ₀−μ_∞)/(1+(Kγ̇)ⁿ) | Johnston et al. (2004) |
| Power-law (W-S) | μ = K γ̇^(n-1) | Walburn & Schneck (1976) |

## Key Code Snippet

```rust
// Millifluidic channel: D = 1 mm, L = 30 mm, Q = 1 mL/min
// For each model, apparent viscosity at γ̇_w = 8V/D feeds Hagen-Poiseuille:
//   ΔP = 128 μ_app L Q / (π D⁴)
```

## Validation Data (Ht ≈ 45 %, T = 37 °C)

| γ̇ [s⁻¹] | μ_exp [mPa·s] |
|---|---|
| 1 | 18 – 25 |
| 10 | 7 – 9 |
| 100 | 4 – 5 |
| 1 000 | ≈ 3.5 |

The example generates a viscosity-vs-shear-rate plot and a pressure-drop comparison table.

## Book Chapter

[← Blood Flow and Rheology Workflows](../biomedical_flows.md)

