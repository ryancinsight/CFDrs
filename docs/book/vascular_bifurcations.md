# Vascular Bifurcations and Stenosis

The vascular bifurcation is the keystone of CFDrs's biomedical flow
support.  It is generated as a CSG composition:

```rust
use cfd_schematics::bifurcation::BifurcationBuilder;
use leto::Point3;

let bif = BifurcationBuilder::new()
    .parent(Point3::new(0.0, 0.0, 0.0), 2.5)            // 5 mm diameter, mm units
    .child_a(Point3::new(0.0, 4.0, 0.0), 1.8, 0.35)    // 3.6 mm, 35° angle
    .child_b(Point3::new(0.0, 4.0, 0.0), 1.8, -0.35)   // mirror child
    .stenosed_a(0.55)                                   // 55% area stenosis
    .build();
```

Boundary conditions follow the Murray-style split, with Reynolds-scaled
flow division between children A and B.  CFDrs composes this CSG with the
1-D reduced solver to compute pressure and shear stress along each
parent / child branch.

### Rheology Choices

```rust
use cfd_1d::blood_rheology::Rheology;

let rheo = Rheology::carreau {          // μ_0, μ_∞, λ, n
    mu_inf: 0.00345,
    mu_0:   0.056,
    lambda: 3.313,
    n:      0.3568,
};
```

Carreau, Carreau-Yasuda, Casson, and Power-law models are all supported.
The rheology is a trait the solver consumes per timestep — switching models
is a single field change.

## Examples Referenced by This Chapter

- [`bifurcation_2d_blood_validation`](examples/bifurcation_2d_blood_validation.md)
  — 2-D Frauhofer-Keller bifurcation with non-Newtonian blood.
- [`venturi_blood_flow_validation`](examples/venturi_blood_flow_validation.md)
  — venturi stenosis under Newtonian and Carreau rheology.
- [`cross_fidelity_branching`](examples/cross_fidelity_branching.md) —
  concurrently run a 1-D, 2-D, and 3-D branch of a vascular bifurcation
  and compare outputs.

## Further Reading

- For microfluidic mixing screens, FDA shear limits, and TPMS scaffolds,
  see [Non-Newtonian Fluids and Blood Rheology](biomedical_flows.md).
- [`cfd-schematics` source](../../crates/cfd-schematics/src/)
