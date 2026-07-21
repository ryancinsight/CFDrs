# Example: blood_poiseuille_2d

**Crate**: `cfd-validation`
**Run**: `cargo run -p cfd-validation --example blood_poiseuille_2d`
**Source**: [`crates/cfd-validation/examples/blood_poiseuille_2d.rs`](../../crates/cfd-validation/examples/blood_poiseuille_2d.rs)

## What This Example Demonstrates

Steady non-Newtonian Poiseuille flow in a 2D channel (200 µm × 200 µm × 1 mm)
solved with the Casson blood model and compared against the Newtonian
analytical Poiseuille solution for the same pressure gradient.

| Aspect | Reference |
|---|---|
| Velocity profile | Newtonian analytical Poiseuille solution |
| Rheology | CassonBlood::normal_blood (yield stress + infinite-shear viscosity) |
| Geometry | 2D channel, µm scale |
| Comparison | `PoiseuilleFlow2D::analytical_solution(mu_newtonian)` vs numerical |

## Key Code Snippet

```rust
use cfd_2d::solvers::{PoiseuilleConfig, PoiseuilleFlow2D, BloodModel};
use cfd_core::physics::fluid::blood::CassonBlood;

// Newtonian reference viscosity: 3.5 mPa·s (canonical whole-blood)
let mu_newtonian = 3.5e-3;

// Casson fluid: tau_y ≈ 0.0056 Pa, mu_inf ≈ 0.00345 Pa·s
let blood = CassonBlood::normal_blood();

let flow = PoiseuilleFlow2D::new(&config, BloodModel::Casson(blood));
flow.solve().expect("Poiseuille solve failed");
let u_numerical = flow.max_velocity();
let u_analytical = flow.analytical_solution(mu_newtonian).max_velocity();
let l2_error = flow.velocity_profile()
    .l2_against(&flow.analytical_solution(mu_newtonian).velocity_profile());
```

## Physics Background

Under Newtonian rheology the Hagen-Poiseuille solution in a 2D channel yields a
parabolic velocity profile. The Casson model introduces a yield stress via

`√τ = √τ_y + √(μ_∞ · γ̇)`

with `τ_y ≈ 0.0056 Pa` (yield stress at Ht ≈ 45%) and `μ_∞ ≈ 0.00345 Pa·s`
(infinite-shear viscosity). The plug-flow core from the yield stress blunts the
profile relative to the Newtonian parabola, so a non-zero L2 error between
numerical (Casson) and analytical (Newtonian) profiles is expected and confirms
that the solver is faithfully reproducing the non-Newtonian constitutive law.

## Book Chapter

[← Validation Suite](../crate_validation.md)
