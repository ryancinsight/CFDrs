# Example: nalgebra_parity (leto-ops canonical)

> **Canonical location**: `D:\atlas\repos\leto\crates\leto-ops\examples\nalgebra_parity.rs`  
> **Crate**: `leto-ops`  
> **Run**: `cargo run --release --example nalgebra_parity -p leto-ops`

The `nalgebra` ↔ `leto` parity harness lives in `leto-ops` — the direct
replacement for `nalgebra`'s linear algebra — not in downstream crates.

## Overview

Solves the 1-D Poisson problem `u'' = -sin(πx)` on `[0,1]` with:

1. **Legacy**: `nalgebra::DMatrix<f64>` + `nalgebra` partial-pivoting LU
2. **Atlas**: `leto_ops::CooMatrix → CsrMatrix` + `leto_ops::SparseLuSolver`

Both paths use identical discrete operators and RHS. Solutions and residuals
are compared to prove numerical equivalence of the Atlas stack.

## Parity Tolerances

| Metric | Tolerance |
|---|---|
| Solution L∞ agreement | ≤ 1e-6 |
| Residual L∞ (both paths) | ≤ 1e-8 |

## Status

✅ **Parity pass** — `SparseLuSolver` (Atlas) produces bit-equivalent
solutions to `nalgebra` LU. `cfd-math` has **zero** `nalgebra` dependencies.

## Part Reference

Part VII — Migration Validation: Legacy ↔ Atlas Parity
