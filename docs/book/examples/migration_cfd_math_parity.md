# Example: nalgebra_parity (Leto Ops)

**Canonical source**: [`leto-ops/examples/nalgebra_parity.rs`](https://github.com/ryancinsight/leto/blob/main/crates/leto-ops/examples/nalgebra_parity.rs)
**Owner**: `leto-ops`
**Run from the Leto checkout**: `cargo run --locked -p leto-ops --example nalgebra_parity`

The legacy `nalgebra` oracle belongs in `leto-ops`, the Atlas source layer
that owns the replacement sparse-matrix and solver APIs. CFDrs consumes those
APIs and does not duplicate the oracle downstream.

## Contract

The harness solves the manufactured one-dimensional Poisson system
`-u'' = sin(πx)` with homogeneous Dirichlet boundaries at 512 interior
points.

- The reference path uses a dense `nalgebra::DMatrix` and partial-pivoting LU.
- The Atlas path assembles Leto COO/CSR storage and solves it with
  `SparseLuSolver`.
- Both providers are checked against the exact discrete sine eigenmode.
- Backward-error limits use `γ(3n)`; forward-error limits use the exact
  infinity-norm condition number of the scaled tridiagonal operator.
- Provider agreement is bounded by the sum of the two forward-error budgets.

The executable prints each observed error and analytical bound, asserts all
five checks, and emits one JSON summary. It emits no timing claim; performance
requires a controlled benchmark.

## Part Reference

Part VII — Migration Validation: Legacy ↔ Atlas Parity
