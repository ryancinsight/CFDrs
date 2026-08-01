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

The physical-unit typing of canonical turbulence state and viscosity outputs is
not silently marked complete by this cleanup. It is a separate Aequitas
boundary item because `cfd-core::TurbulenceModel` currently exposes scalar
fields and k-epsilon kinetic energy needs a provider dimension that is not yet
present in the current Aequitas quantity set.
