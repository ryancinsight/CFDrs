# 1-D Biomedical Flows

Reduced-order 1-D haemodynamic examples live in `cfd-1d`.  They cover the
core biomedical-flow library: bifurcations, microfluidics, FDA device
screening, TPMS scaffolds, cavitation analysis, and hemolysis assessment.

| Example | Physics |
|---------|---------|
| `blood_bifurcation` | Murray's law + Casson model at a symmetric bifurcation |
| `microfluidic_chip` | Pressure-driven flow through a multi-channel chip |
| `fda_shear_limit_screening` | FDA haemolysis threshold screening (shear history) |
| `venturi_parallel_analysis` | Parallel Venturi network for pressure drop mapping |
| `tpms_blood_1d` | Gyroid-inspired TPMS 1D network with Carreau-Yasuda blood |
| `cavitation_venturi_analysis` | SDT device design: σ per channel, regime classification |
| `medical_millifluidic_screening` | Full medical screening: hemolysis + WSS + FDA + cavitation |
| `hemolysis_serpentine_analysis` | Giersiepen–Wurzinger hemolysis model on serpentine geometry |

```bash
cargo run -p cfd-1d --example blood_bifurcation
cargo run -p cfd-1d --example tpms_blood_1d
cargo run -p cfd-1d --example cavitation_venturi_analysis
cargo run -p cfd-1d --example medical_millifluidic_screening
cargo run -p cfd-1d --example hemolysis_serpentine_analysis
```
