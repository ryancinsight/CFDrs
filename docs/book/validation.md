# Chapter 7 — Validation and Benchmarking

Computational accuracy is established by comparing against **analytical
solutions**, **experimental data**, and **published benchmark results**.
CFDrs ships a `cfd-validation` crate dedicated to this discipline.

## Benchmark Hierarchy

| Level | Description | CFDrs coverage |
|-------|-------------|---------------|
| Analytical | Exact solution exists (Poiseuille, Stokes) | ✓ Full |
| Semi-analytical | Similarity / series solutions (Blasius, Graetz) | ✓ Partial |
| Well-documented numerical | Ghia lid-driven cavity; Moser DNS profiles | ✓ Full |
| Experimental | FDA nozzle benchmark; blood bifurcation | ✓ Blood flows |

## Convergence Studies

Richardson extrapolation establishes the **observed order of accuracy**:

```text
φ_exact ≈ φ_fine + (φ_fine - φ_coarse) / (r^p - 1)
```

where *r* is the mesh refinement ratio and *p* the expected order.

```rust
use cfd_validation::convergence::RichardsonExtrapolation;

let rich = RichardsonExtrapolation::new(&coarse, &fine, refinement_ratio: 2.0);
let (order, exact_estimate) = rich.compute()?;
```

The `richardson_convergence` example confirms **2nd-order accuracy** for the
finite-volume stencils and **spectral convergence** for the FFT-based spectral
Poisson solver.

## Lid-Driven Cavity (Ghia 1982)

The canonical incompressible Navier-Stokes benchmark: unit-square cavity with
Re = 100, 400, 1000, 3200.

CFDrs reproduces the Ghia *u*-velocity profile along the vertical centreline
to within 0.5% for Re ≤ 1000:

```rust
cargo run -p cfd-suite --example cavity_validation
```

## Hagen-Poiseuille (Pipe Flow)

For a circular pipe of radius *R* and pressure gradient dP/dx:

```text
u(r) = (1/(4μ)) · (-dP/dx) · (R² - r²)
```

```rust
cargo run -p cfd-suite --example pipe_flow_validation
```

## Blood-Flow FDA Benchmarks

CFDrs includes the FDA nozzle haemolysis model and hemolysis
(Giersiepen-Wurzinger) as validated metrics:

```rust
cargo run -p cfd-suite --example blood_flow_1d_validation
cargo run -p cfd-1d --example fda_shear_limit_screening
```

## Examples

- [comprehensive_validation_suite](examples/comprehensive_validation_suite.md)
- [richardson_convergence](examples/richardson_convergence.md)
- [cavity_validation](examples/cavity_validation.md)
- [pipe_flow_validation](examples/pipe_flow_validation.md)
- [blood_poiseuille_2d](examples/blood_poiseuille_2d.md)
- [turbulent_channel_flow](examples/turbulent_channel_flow.md)
