# Validation Suite

CFDrs ships three crate-level validation examples that run inside `cfd-validation`.
Together they confirm solver accuracy, convergence order, and haemodynamic physics.

| Example | Validates |
|---------|-----------|
| `comprehensive_validation_suite` | End-to-end solver across all supported flow regimes |
| `richardson_convergence` | Second-order spatial convergence via Richardson extrapolation |
| `blood_poiseuille_2d` | 2D Poiseuille blood-flow against Hagen-Poiseuille analytical solution |

Run all three with:
```bash
cargo run -p cfd-validation --example comprehensive_validation_suite
cargo run -p cfd-validation --example richardson_convergence
cargo run -p cfd-validation --example blood_poiseuille_2d
```
