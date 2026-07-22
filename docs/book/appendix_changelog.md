# CFDrs Changelog

## Unreleased

### Documentation
- Atlas-Stack Integration Part (Part VII) with nine migration sub-chapters
  (Eunomia, Arrays, Geometry, SIMD, Memory, Concurrency, FFT, GAT
  Tiling, Validation).
- Atlas glossary appendix (`appendix_glossary.md`) covering every Atlas
  crate.
- Migration type reference expanded with concrete call-site rewrites and
  validation checklists.
- Forward roadmap (`BOOK_ORGANIZATION.md`) — Phase 1 (Foundation)
  complete, Phase 2 (Documentation) current, Phase 3 (Atlas Migration)
  mid-flight, Phase 4 (Cleanup) gate, Phase 5 (optional GPU) forward.

### Atlas Migration Status
- Completed: eunomia traits across `cfd-core`, `cfd-1d`, `cfd-2d`,
  `cfd-3d`, `cfd-math`, `cfd-validation`, `cfd-schematics`.
- In-progress: `leto::NdArray` across `cfd-math`, `cfd-1d`,
  `cfd-2d`, `cfd-3d`.
- In-progress: `mnemosyne::Arena` sub-arenas in `cfd-1d`, `cfd-3d`.
- In-progress: `hermes-simd` lanes in `cfd-2d::stencil`,
  `cfd-3d::fem`.
- Graduating: `apollo::FftPlan` in `cfd-3d::spectral`.

### Migration Phase Plan
- Phase 3 § Migrations: continue bulk migration. See `BOOK_ORGANIZATION.md`.
- Phase 4 § Cleanup: prune legacy `ndarray`/`nalgebra`/`rayon`/`tokio`/
  `rustfft` after parity tests pass.

## v0.0.1

- Initial CFDrs scaffold.
- Atlas dependency graph drafted (see `appendix_dependencies.md`).
