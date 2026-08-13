# cfd-optim

Design-space search for millifluidic layouts, part of the
[CFDrs](https://github.com/ryancinsight/CFDrs) simulation suite.

It takes a `cfd-schematics` blueprint, generates and scores candidate variants
against physical and manufacturing constraints, and reports the evidence behind the
selection.

## What it provides

- **`constraints`** — the public constraint set, including literature-anchored limits
  such as shear-induced platelet activation thresholds.
- Objectives, orchestration, and search strategies re-exported from the crate root:
  deterministic sweeps, a genetic search with explicit mutation operators, and an
  evaluated candidate pool.
- `robustness_sweep_blueprint` with `STANDARD_PERTURBATIONS` for sensitivity analysis.
- Reporting: evidence records, figure manifests, and narrative generation, so a
  design decision ships with its supporting artifacts rather than a bare score.

Stratified sampling is backed by `tyche` for reproducible design draws.

## Documentation

- [CFDrs book](https://ryancinsight.github.io/CFDrs/)
- [Workspace README](https://github.com/ryancinsight/CFDrs#readme)

## License

MIT OR Apache-2.0
