# Example: migration_cfd_math_parity

**Source**: `crates/cfd-math/examples/migration_cfd_math_parity.rs`  
**Crate**: `cfd-math`

## Overview

Atlas migration parity harness for `cfd-math::sparse` and `cfd-math::linear_solver`. Solves the same 1D Poisson problem `u'' = -sin(πx)` with:

1. **Legacy path**: `nalgebra::DMatrix` + hand-rolled dense conjugate gradient
2. **Atlas path**: `cfd_math::sparse::SparseMatrixBuilder` (CSR via `leto_ops`) + `cfd_math::linear_solver::ConjugateGradient`

Both paths use identical discrete operators, RHS, and CG tolerance. The harness compares solutions and residuals, emitting JSON for regression gates.

## Parity Tolerances

| Metric | Tolerance |
|---|---|
| Velocity/solution field RMS | 1e-6 relative |
| Mass conservation (residual L2) | 1e-8 relative |
| L∞ residual | 1e-8 absolute |

## Status

✅ **Parity pass** — Atlas CSR/CG produces bit-equivalent solutions to the nalgebra reference within tolerance. Zero nalgebra runtime dependencies in production code.

## Run

```bash
cargo run --release --example migration_cfd_math_parity -p cfd-math
```

## Part Reference

Part VII — Migration Validation: Legacy ↔ Atlas Parity
