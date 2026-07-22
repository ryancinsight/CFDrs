# Validation Report — CFDrs linear algebra provider parity

## Ownership decision

The `nalgebra` comparison harness is owned by
[`leto-ops/examples/nalgebra_parity.rs`](https://github.com/ryancinsight/leto/blob/main/crates/leto-ops/examples/nalgebra_parity.rs).
Leto Ops owns the sparse storage and solver seam that replaces the legacy
provider. CFDrs consumes that seam and does not maintain a second oracle.

A source and manifest scan of this checkout finds no direct `nalgebra`,
`ndarray`, `num-traits`, `num-complex`, `rustfft`, `rayon`, `tokio`,
or `packed_simd` dependency in CFDrs Rust targets. The active
`cfd-math` manifest consumes `eunomia`, `leto`, `leto-ops`, `moirai`,
and `hermes-simd`.

## Numerical contract

The canonical harness solves

\[
-u''(x) = \sin(\pi x),\qquad x\in(0,1),\qquad u(0)=u(1)=0
\]

on a centered finite-difference grid.

| Path | Storage | Solver |
|---|---|---|
| Reference | dense `nalgebra::DMatrix<f64>` | partial-pivoting LU |
| Atlas | Leto `CooMatrix<f64>` → `CsrMatrix<f64>` | `SparseLuSolver` |

The executable uses 512 interior points; its co-located regression test uses
64. Both paths assemble the same scaled tridiagonal operator and right-hand
side.

The five behavioral checks are:

1. normalized reference backward error ≤ `γ(3n)`;
2. normalized Atlas backward error ≤ `γ(3n)`;
3. reference error against the exact discrete sine eigenmode ≤ the
   condition-number-derived forward bound;
4. Atlas error against that eigenmode ≤ the same bound;
5. provider disagreement ≤ twice that forward bound.

The exact infinity-norm condition number is
`2 ceil(n/2) (n + 1 - ceil(n/2))`. The continuum discretization error is
reported separately and is not confused with provider error.

## Reproduction

From the Leto checkout:

```sh
cargo run --locked -p leto-ops --example nalgebra_parity
```

The runnable example prints the observed error and bound for every check, then
emits a JSON summary containing the problem order, check count, and pass state.
It intentionally emits no timing or memory-speedup claim. Such claims require a
controlled benchmark and allocation measurement.

## CFDrs evidence boundary

This report establishes the provider ownership and the mathematical parity
contract. The corresponding CFDrs migration contract remains documented in
[`docs/book/migration_validation.md`](../docs/book/migration_validation.md).
CFDrs verification separately checks its own callers and examples against the
pinned provider graph.
