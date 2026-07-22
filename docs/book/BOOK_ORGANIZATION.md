# Appendix E — Book Organization Forward Roadmap

> **Authoritative TOC:** [`SUMMARY.md`](SUMMARY.md) is the single source of
> truth for the *current* book (the file mdBook renders).  This document is
> the *forward roadmap* — its Parts and Chapters describe intended future
> expansion and may run ahead of the chapters and examples that exist today.

## Overview

This document outlines the planned structure for **CFDrs: Computational
Fluid Dynamics, Vascular Flows, and Atlas-Stack Migration** — following
the kwavers and helios book models.

CFDrs is divided into three layers:

1. **CFDrs (this book)** — physics-anchored, example-driven CFD reference.
2. **Atlas stack** — the underlying trait/memory/array/SIMD/FFT layer;
   CFDrs is mid-flight in a bulk-then-cleanup migration to it (see
   [Atlas Dependency Map](appendix_dependencies.md)).
3. **Cross-domain atlases** — [`kwavers`](../../../kwavers/docs/book/README.md)
   (computational ultrasound) and [`helios`](../../../helios/docs/book/README.md)
   (radiotherapy simulation), both of which already consume the Atlas stack
   that CFDrs is migrating to.

## Book Title

**CFDrs: Computational Fluid Dynamics, Vascular Flows, and Atlas-Stack Migration.**

## Parts at a Glance

The authoritative TD summary is in [`SUMMARY.md`](SUMMARY.md).  The forward
parts are:

- **Part I — Foundations**: problem setup, core types, simulation lifecycle.
- **Part II — Core Flow Cases and Validation**: lid-driven cavity, Hagen–
  Poiseuille pipe, turbulent channel.
- **Part III — Turbulence, Multiphase, and Cavitation**: closure models,
  multiphase coupling, Rayleigh–Plesset cavitation.
- **Part IV — Biomedical and Specialized Flows**: bifurcations, venturi,
  serpentine mixing, FDA shear-limit screening, TPMS scaffolds.
- **Part V — Numerical Methods and Solvers**: Stencil / FEM / Spectral
  discretizations, matrix-free, Krylov solvers, time integrators.
- **Part VI — Geometry, Meshing, and CSG**: CSG primitives, mesh
  generation, fluid–solid mapping.
- **Part VII — Atlas Stack Integration (Migration Reference)** (this
  document's focus): 9 migration sub-chapters covering the bulk migration
  from `ndarray`/`nalgebra`/`tokio`/`rayon`/`rustfft` to
  `eunomia`/`leto`/`hermes`/`mnemosyne`/`themis`/`moirai`/`apollo`.
- **Part VIII — Crate-Level Examples** (cfd-1d / cfd-2d / cfd-3d /
  cfd-validation / cfd-optim): every example belongs to its producing
  crate.
- **Appendix**: Atlas dependency map, migration type reference, glossary,
  changelog.

## Implementation Status

### Currently Implemented (CFDrs)

- `cfd-core` — scalar seams, error types, boundary condition contracts.
- `cfd-math` —dense/sparse linalg, Krylov solvers, FFT bridge.
- `cfd-1d`, `cfd-2d`, `cfd-3d` — 1-D / 2-D / 3-D solvers.
- `cfd-validation` — orthogonal oracle harness (Ghia, Hagen–Poiseuille).
- `cfd-schematics` — CSG primitives (venturi, bifurcation, serpentine).
- `cfd-io` — HDF5 / NetCDF / CSV adapters.
- `cfd-optim` — optimization audit infrastructure.
- `cfd-python` — PyO3 bindings.

### Migration Status

- **eunomia**: COMPLETE — all numeric bounds ported.
- **leto arrays**: in-progress — `cfd-math`, `cfd-1d`, `cfd-2d`, `cfd-3d`.
- **leto geometry**: in-progress — `cfd-schematics` partial.
- **hermes-simd**: in-progress — `cfd-2d::stencil`, `cfd-3d::fem`.
- **mnemosyne**: in-progress — sub-arenas in `cfd-1d`, `cfd-3d`.
- **themis**: in-progress — `Placement::numa_aware()` wired in `cfd-3d`.
- **moirai**: in-progress — `cfd-2d::parallel`, `cfd-3d::parallel`.
- **apollo**: graduating — `cfd-3d::spectral`, `cfd-2d::spectral`.
- **coeus**: n/a — no autodiff in CFDrs.
- **hephaestus**: opt-in — GPU acceleration when wired.
- **ritk**: n/a — no image data in CFDrs.

### Examples Inventory (per crate)

Total CFDrs `.rs` examples: ~67

- `cfd-1d/examples/` — 10
- `cfd-2d/examples/` — 6
- `cfd-3d/examples/` — 4
- `cfd-validation/examples/` — 3
- `cfd-optim/examples/` — 6
- workspace `examples/` — ~38 (utilities, demos, integration)

The book documents all 67+ examples across Parts I–VIII.

## Phased Roadmap

### Phase 1: Foundation (Complete)

- [x] Workspace scaffolding (`cfd-core`, `cfd-math`, `cfd-{1,2,3}d`, etc.).
- [x] Legacy linalg surface (`ndarray` / `nalgebra`).
- [x] Reference validation harness (Ghia, Hagen–Poiseuille).
- [x] CSG primitive catalogue.
- [x] PyO3 boundary.

### Phase 2: Documentation Book (Current Focus)

- [x] Build [`SUMMARY.md`](SUMMARY.md) (8 parts + Appendix).
- [x] Each Part has a top-level chapter (`foundations.md`, `core_flows.md`,
      `turbulence_multiphase.md`, `biomedical_flows.md`,
      `numerics_and_solvers.md`, `geometry_and_meshing.md`,
      `performance_and_atlas.md`).
- [x] Each Part has its example chapter set.
- [x] Atlas-Stack Integration Part (this roadmap's anchor).
- [x] Atlas Dependency Map appendix.
- [x] Migration Type Reference appendix.
- [x] Glossary appendix.
- [f] Run `mdbook build` and resolve any dead-link warnings.

### Phase 3: Atlas Stack Migration (Mid-Flight)

- [f] Promote `cfd-math` `Linalg` crate to `leto::NdArray` everywhere.
- [f] Port `csr` sparse kernels to `leto::CsrMatrix` + `hermes-simd`
      lanes.
- [f] Wire all FFT calls through `apollo::FftPlan`.
- [f] Replace `Vec<T>` accumulators with `Arena<T>` sub-arenas.
- [f] Replace `rayon::par_iter` with `moirai::Executor::scope`.
- [f] Replace `Iterator<Item = &Tile<'a>>` with `LendingIterator`.
- [f] Bind `cfd-3d`, `cfd-2d` to `moirai` pools via `themis::Placement`.

### Phase 4: Validation and Pruning (Cleanup)

- [ ] Run parity harness across all canonical scenarios (Ghia, Hagen–
      Poiseuille, venturi, bifurcation, serpentine hemolysis).
- [ ] Document any residual gap in `validation_reports/` and close it.
- [ ] Benchmark Hermes-vectorized SpMV vs legacy baseline.
- [ ] Benchmark Apollo-FT vs `rustfft` baseline.
- [ ] Benchmark `moirai` parallel sweep vs `rayon` baseline.
- [ ] Prune unused `ndarray`/`nalgebra`/`rayon`/`tokio`/`rustfft` deps;
      record the cleanup commit.
- [ ] Re-publish the Atlas Dependency Map reflecting "no legacy crates
      anymore".

### Phase 5 (Forward): GPU and Structured Solvers (Optional)

- [ ] Wire `hephaestus` GPU kernels for `cfd-3d::stencil`.
- [ ] Wire `hephaestus` GPU kernels for `cfd-2d::stencil`.
- [ ] Document GPU vs CPU parity (matrix-free, SpMV, stencil).

## Build Instructions

```bash
# Build the book
mdbook build docs/book

# Serve with hot reload
mdbook serve docs/book
```

## Organization Principles

1. **Executable Examples**: each chapter cross-references `.rs` examples
   in `crates/*/examples/` or workspace `examples/`.
2. **SRP / SoC**: each solver crate owns a single fidelity layer; CSG
   geometry is split from solver kernels.
3. **SSOT / DRY**: Atlas crates own one trait frontier each — no
   parallelism between Atlas crates, no overlap.
4. **Zero-Cost Abstractions**: const generics, ZSTs, GATs, `Cow`
   encode capability and lifetime variance at the type level.
5. **DIP**: depend on [`eunomia`] traits, not concrete impls.

## Cross-References

- [`kwavers` book](../../../kwavers/docs/book/README.md) — Part VI
  Atlas migration is the model for `CFDrs`'s Part VII.
- [`helios` book](../../../helios/docs/book/README.md) — Part VIII
  Atlas migration is the model for `CFDrs`'s Part VII.
- [Atlas Dependency Map](appendix_dependencies.md)
- [Migration Type Reference](appendix_migration.md)
