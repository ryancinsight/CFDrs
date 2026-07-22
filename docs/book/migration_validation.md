# Chapter 17 — Migration Validation: Legacy ↔ Atlas Parity

The Atlas migration is governed by a **parity test harness** that runs
both the legacy and the Atlas ports of each solver on the same input and
demands agreement within a documented tolerance.  Until a solver passes
parity, it does not graduate from the legacy crate to the Atlas-only path.

## Parity Test Anatomy

A parity test has three components:

1. **Reference input**. A canonical scenario (Ghia cavity, Poiseuille pipe,
   serpentine hemolysis, Venturi cavitation).  Inputs are version-pinned in
   `cfd-validation/data/`.
2. **Legacy baseline**. The pre-migration sweep, run via the legacy
   `ndarray`-/`nalgebra`-based crate.  Outputs are checked into the
   `validation_reports/` directory with cryptographic hashes.
3. **Atlas candidate**. The migrated port running the **exact same** input
   through the new Atlas stack.  Outputs are byte-compared against the
   baseline hash under a tolerance budget.

## Tolerance Budgets

| Quantity | Default tolerance |
|---|---|
| Velocity field RMS | 1e-6 relative |
| Pressure field RMS | 1e-5 relative |
| Vorticity L∞ | 1e-4 absolute |
| Mass conservation | 1e-8 relative |
| Energy conservation | 1e-7 relative |
| Cavitation inception index | exact match |

The tolerance table is part of the migration contract — every PR that closes
a migration gap must update both the implementation and the validation
report.

## The CFDrs Parity Harness

```rust
// pseudo-Rust
fn parity_test<F: FloatElement>(scenario: &Scenario) -> ParityReport {
    let legacy_out = legacy::run(scenario.clone());
    let atlas_out  = atlas::run::<F, _>(scenario.clone());
    ParityReport::compare(&legacy_out, &atlas_out)
}
```

The harness lives in `cfd-validation` and is invoked by
[`comprehensive_validation_suite`](examples/comprehensive_validation_suite.md).

## Migration Status Updates

Each migration that **graduates** to legacy-prune status produces a
`validation_reports/migration_<crate>.md` document recording:

- The parity-test scenarios used.
- The legacy vs. Atlas timings.
- The error budget and observed error.
- The Atlas crates involved.
- The legacy crates that can now be removed.

## Validation Examples

- [`comprehensive_validation_suite`](examples/comprehensive_validation_suite.md) —
  the canonical parity sweep.
- [`richardson_convergence`](examples/richardson_convergence.md) —
  convergence order under Atlas linalg.
- [`cavity_validation`](examples/cavity_validation.md) — Re = 100
  lid-driven cavity, Ghia et al. oracle.
- [`pipe_flow_validation`](examples/pipe_flow_validation.md) — Hagen–
  Poiseuille.
- [`blood_poiseuille_2d`](examples/blood_poiseuille_2d.md) — non-Newtonian
  pipe flow.

## Further Reading

- [`cfd-validation` source](../../crates/cfd-validation/src/)
- [Atlas Dependency Map](appendix_dependencies.md)
- [Migration Overview](migration_overview.md)
