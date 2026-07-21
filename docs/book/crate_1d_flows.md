# 1-D Biomedical Flows

Reduced-order 1-D haemodynamic examples live in `cfd-1d`.  They cover the
core biomedical-flow library: bifurcations, microfluidics, FDA device
screening, and TPMS scaffolds.

| Example | Physics |
|---------|---------|
| `blood_bifurcation` | Murray's law + Casson model at a symmetric bifurcation |
| `microfluidic_chip` | Pressure-driven flow through a multi-channel chip |
| `fda_shear_limit_screening` | FDA haemolysis threshold screening (shear history) |
| `venturi_parallel_analysis` | Parallel Venturi network for pressure drop mapping |

```bash
cargo run -p cfd-1d --example blood_bifurcation
cargo run -p cfd-1d --example microfluidic_chip
cargo run -p cfd-1d --example fda_shear_limit_screening
cargo run -p cfd-1d --example venturi_parallel_analysis
```
