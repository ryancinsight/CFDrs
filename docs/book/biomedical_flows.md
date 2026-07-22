# Chapter 4 — Blood Flow and Rheology Workflows

CFDrs supports a focused subset of physiologically relevant geometries —
**vascular bifurcations**, **venturi / stenosis regions**, **microfluidic
screens**, and **TPMS (triply-periodic minimal-surface) scaffolds**.
Each of these geometries couples a `cfd-schematics` CSG primitive to a
low-dimensional CFD solve (1-D, 2-D, or 3-D), with **non-Newtonian
rheology** for blood.

## The Vascular Bifurcation

A bifurcation is generated as a CSG composition:

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

## Microfluidic Mixing Screens

Microfluidic screens are serpentine or bifurcation arrays designed for
mixing, separation, or hemolysis screening.  CFDrs handles them as a chain
of bifurcations:

```rust
use cfd_schematics::serpentine::SerpentineBuilder;

let chip = SerpentineBuilder::new()
    .section_count(20)
    .channel_radius(0.25)
    .turning_radius(2.5)
    .build();
```

Each section is solved with [`cfd-3d`] time-stepping, then stitched
together for the chip-scale perfusion metric.

## FDA Shear-Limit Screening

CFDrs exposes a **shear-limited screening** API that walks a serpentine /
bifurcation chain and reports branches with shear stress above an FDA-style
threshold.  The API is intentionally minimal so that custom geometries can
plug in:

```rust
let report = cfd_1d::screening::shear_violations(&state, &geometries)?;
for branch in report.violations {
    eprintln!("branch {}: peak shear {} Pa (limit 9.0)", branch.id, branch.peak);
}
```

## TPMS Scaffolds

Tissue scaffolds modeled as TPMS (Schwarz P, Gyroid, ...) are CSG-built
and tiled for `cfd-3d` FEM assembly.  Boundary conditions are velocity-driven
at the inlet and outflow at the outlet.

## Examples Referenced by This Chapter

Part IV covers six example chapters:

- [`blood_flow_1d_validation`](examples/blood_flow_1d_validation.md) —
  1-D reduced model, validated against a reference phantom.
- [`bifurcation_2d_blood_validation`](examples/bifurcation_2d_blood_validation.md)
  — 2-D Frauhofer-Keller bifurcation with non-Newtonian blood.
- [`blood_rheology_models`](examples/blood_rheology_models.md) — three
  rheology models compared on the same patient-specific geometry.
- [`venturi_blood_flow_validation`](examples/venturi_blood_flow_validation.md)
  — venturi stenosis under Newtonian and Carreau rheology.
- [`serpentine_mixing_comprehensive`](examples/serpentine_mixing_comprehensive.md)
  — serpentine mixer design space, 1-D reduced to chip-scale.
- [`cross_fidelity_branching`](examples/cross_fidelity_branching.md) —
  concurrently run a 1-D, 2-D, and 3-D branch of a vascular bifurcation
  and compare outputs.

## Further Reading

- [`cfd-schematics` source](../../crates/cfd-schematics/src/)
- [`cfd-1d` source](../../crates/cfd-1d/src/)
- [Turbulence and Multiphase Flow](turbulence_multiphase.md) for cavitation
  modelling relevant to stenosed regions.
