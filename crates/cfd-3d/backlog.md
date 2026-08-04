# CFD-3D Update Backlog

## Closed historical compilation compatibility slice — 2026-07-31

The old `cfd-mesh` compatibility checklist is reconciled with the current
source. The listed trait, mesh-accessor, builder, matrix-bound, and iterator
diagnostics are no longer present. The source-level turbulence placeholder was
also removed; the public module now re-exports the canonical Eunomia-backed
implementations from `physics::turbulence`.

The current local locked build remains blocked before compilation by the shared
Atlas overlay requesting a Cargo.lock refresh for mutable provider patches. The
root cfd-3d audit records the pinned package-check evidence; this residual is
environmental lock coherence, not an unresolved cfd-3d compiler error.

The physical-unit typing of canonical turbulence outputs is closed by
CFDRS-AEQ-MET-44. `cfd-core::TurbulenceModel` now returns Aequitas
`KinematicViscosity<T>` and `SpecificEnergy<T>` values; every cfd-3d closure
wraps its real scalar result at the public boundary, while solver state remains
scalar and formulas extract explicitly. The provider supplies coherent `J/kg`
semantics, and no imaginary unit is introduced. The local locked compile
remains blocked before rustc by the shared Atlas overlay requesting a provider
lock refresh; hosted CI is the remaining integration gate.
