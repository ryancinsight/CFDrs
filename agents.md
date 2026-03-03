# CFDrs Workspace — Agent Reference

> **Owner**: Ryan Clanton (`@ryancinsight`)  
> **Current Sprint**: 1.90.0 — SDT Therapy + Leukapheresis Pipeline
> **Quality Gates**: Build 0 warnings ✅ · Tests 100% pass ✅ · Clippy 0 production ✅ · 0 technical debt ✅

---

## Purpose

`CFDrs` is a **modular Rust CFD simulation suite** targeting millifluidic / microfluidic device design — specifically Sonodynamic Therapy (SDT) devices constrained to a 96-well plate footprint (ANSI/SLAS 1-2004, 127.76 × 85.47 mm).

The hierarchy of concerns (in order of precedence):

```
Mathematical Proofs → Formal Verification → Empirical Validation → Production Deployment
```

---

## Agent Doc Audit Status

- Crate-level `agents.md` files (`crates/*/agents.md`) were audited against `Cargo.toml`, `src/lib.rs`, and top-level `src/` trees on **2026-02-26**.
- `crates/cfd-mesh/agents.md` was converted to a hybrid format: current agent-reference snapshot up top, historical rewrite proposal notes preserved below.
- `cfd-optim` documentation claims were refreshed to current code-truth counts (`DesignTopology` enum variants and `SdtMetrics` fields).
- Workspace warning cleanup completed on **2026-02-26**:
  - Removed non-root profile config from `crates/cfd-python/Cargo.toml`.
  - Corrected root profile package key to `cfd-python`.
  - `cargo check --workspace --no-default-features` now completes with **0 warnings**.

---

## Workspace Layout

```
Cargo.toml            workspace root; resolver = "2"
crates/
  cfd-core            Core traits, boundary conditions, GPU/CPU dispatch, fluid properties
  cfd-math            Linear solvers (GMRES/CG/BiCGSTAB/AMG), SIMD SpMV, spectral/DG/WENO
  cfd-mesh            Half-edge mesh, CSG Boolean ops, watertight validation, OpenFOAM I/O
  cfd-io              VTK/HDF5/CSV read-write, checkpointing, parallel I/O
  cfd-schematics      Design-time topology, geometry generation, 2D schematic rendering
  cfd-1d              1D lumped-network Hagen-Poiseuille solver for channel networks
  cfd-2d              2D incompressible N-S: SIMPLE/PISO, LBM, FDM/FVM
  cfd-3d              3D FEM/IBM/Level-Set, LES/DES turbulence, multiphase (VOF + LS)
  cfd-optim           SDT therapy + leukapheresis optimiser: 24 `DesignTopology` enum variants, wave-channel SVG, optional mesh-export pipeline
  cfd-validation      MMS, Richardson extrapolation, Ghia cavity, analytical benchmarks
  cfd-python          PyO3 Python bindings (cfd-python)
xtask/               Cargo xtask automation
docs/                ADR, SRS, PRD, sprint summaries, gap analyses, backlog, checklist
examples/            Root-package examples (pipe_flow_1d, blood_millifluidic, venturi, CSG, …)
benches/             Criterion benchmarks
tests/               Integration tests
.agent/rules/        persona.md — always-on agent persona
```

### Crate Dependency Graph (simplified, no cycles)

```
cfd-core
  └─► cfd-math
        └─► cfd-mesh
              └─► cfd-io
                    └─► cfd-schematics
                              ├─► cfd-1d
                              └─► cfd-2d
                                    └─► cfd-3d
                                          ├─► cfd-optim
                                          └─► cfd-validation
```

All crates depend on `cfd-core`. `cfd-python` wraps the workspace for PyO3.

---

## Active Phase — `cfd-mesh` Rewrite

`crates/cfd-mesh/agents.md` tracks the crate-level rewrite plan in detail.

| Phase | Status | Deliverable |
|-------|--------|-------------|
| 1–Foundation | ✅ DONE | `Mesh<'id>` + `with_mesh()` + GhostCell brand system |
| 2–Geometry | ✅ DONE | Exact Shewchuk predicates; unified 26-neighbour SnappingGrid |
| 3–Topology | ✅ DONE | Half-edge kernel (twin/next/prev invariants); SlotMap keys |
| 4–Validation | ✅ DONE | Manifold + Euler–Poincaré + winding + patch completeness |
| 5–NURBS | ✅ DONE | Cox–de Boor; adaptive tessellation |
| 6–Builders | ✅ DONE | VenturiMeshBuilder, SerpentineMeshBuilder, BranchingMeshBuilder → `IndexedMesh` |
| 7–I/O | ✅ DONE | STL ASCII/binary (read + write); VTK; cargo-fuzz harness |
| 8–Boolean | ✅ DONE | BVH broad-phase + exact T-T intersection + S-H clipping + reconstruction |
| 9–Docs | ✅ DONE | Rustdoc + aquamarine + 0 `missing_docs` warnings |
| 10–OpenFOAM | ✅ DONE | `write_openfoam_polymesh()` (IndexedMesh) + `_he()` (HalfEdgeMesh); all 5 polyMesh files; `PatchSpec`; `face_patch()` / `patch_info()` on `HalfEdgeMesh`; 113 tests |
| 11–Next | 🔲 TODO | See `backlog.md` for Phase 11 candidates |

---

## Recent Work — Phases 0–7

| Phase | Summary |
|-------|---------|
| 0 — Build Fixes | Fixed 15 examples across the workspace so they compile and run cleanly. |
| 1 — cfd-3d Examples | Added 4 new 3D examples (`venturi_3d_cavitation`, `bifurcation_3d_blood`, `serpentine_3d_dean`, `spectral_poisson_3d`) + `domain_solver_validation` integration test. |
| 2 — Cross-Fidelity | Added 3 cross-fidelity validation examples bridging 1D/2D/3D solvers. |
| 3 — External Validation | Added 2 external validation examples + integration test for literature comparisons. |
| 4 — cfd-optim Pipeline | Rewired `mesh-export` feature from `blue2mesh` to `cfd-mesh`; created `pipeline.rs` module (`DesignPipeline`, `DesignArtifacts`). |
| 5 — Optim Examples | Added `sdt_pipeline`, `sdt_top5_final`, and `sdt_2d_validation` examples in `cfd-optim`. |
| 7 — Integration Tests | Added 3 new integration test files for end-to-end validation. |

---

## Coding Conventions

### Module Size

Production modules: **< 500 lines** (hard limit). Test files: ≤ 565 lines (advisory).  
If a module approaches the limit, split by bounded context before adding new code.

### Quality Gate Checklist (run before every commit)

```sh
# Zero build warnings
cargo build --workspace --no-default-features 2>&1 | grep -c "^warning"   # expect 0

# Zero production clippy warnings
cargo clippy --workspace --no-default-features --lib --bins \
  -- -W clippy::all -W clippy::pedantic 2>&1 | grep -c "^warning"         # expect 0

# 100% tests pass
cargo test --workspace --no-default-features

# Feature-gated build
cargo build --workspace --features csg
cargo test  --workspace --features csg
```

### Technical Debt — Zero Tolerance

**Prohibited in all source files**:
- `TODO`, `FIXME`, `XXX` comments
- `unimplemented!()`, `todo!()` macros
- Placeholder / stub implementations
- Temporary workarounds documented as such

Every line of code must be mathematically justified and production-complete.

### Naming

- No `Enhanced*`, `Optimized*`, `New*`, `V2*` variants. Delete old code; replace it.
- Domain-relevant names that reveal architectural role.
- Ubiquitous language from the CFD / microfluidics domain.

---

## Architecture Principles

### Clean Architecture layers (unidirectional)

```
Domain      → pure business logic; entities, value objects, aggregates; no external deps
Application → use cases; command/query handlers
Infrastructure → repos, external adapters, file I/O, GPU kernels
Presentation → PyO3 bindings, CLI, examples
```

Crate boundaries **are** bounded contexts. No circular imports. No leaking internals.

### Type-System First

- Newtypes for physical quantities (pressure, velocity, viscosity)
- Typestates for builder patterns (prevent calling `.build()` before `.with_boundary()`)
- Trait-driven APIs; prefer generics over dynamic dispatch in hot paths
- `SlotMap` generational keys (not raw `u32`) for mesh entity handles

### Memory Safety Patterns

- Zero-copy: iterators, slices, `Cow`, `rkyv` serialisation
- Arena allocation for short-lived mesh temporaries
- `GhostCell<'id, T>` + `GhostToken<'id>` for branded compile-time access control
- Avoid `Arc`/`Mutex` in domain layer; prefer ownership transfer

### Parallelism

- `rayon` for data-parallel mesh operations and solver loops
- `tokio` async for I/O-bound tasks (checkpointing, external coupling)
- GPU via `wgpu` behind the `gpu` feature flag (default on)

---

## Mathematics & Verification

Every non-trivial algorithm must include:

1. **`# Theorem`** + **`**Proof sketch**`** in Rustdoc  
2. A **property-based test** (`proptest`) encoding the theorem  
3. Validation against an **analytical solution** (not just "no crash")

Key theorems already in the codebase:

| Location | Theorem |
|----------|---------|
| `cfd-mesh: topology/halfedge.rs` | Twin involution: `twin(twin(he)) == he` |
| `cfd-mesh: watertight/check.rs` | Euler–Poincaré: `V - E + F = 2(1-g)` |
| `cfd-mesh: topology/manifold.rs` | 2-manifold: 1–2 faces per edge; vertex stars are disks |
| `cfd-mesh: geometry/predicates.rs` | Shewchuk exact orient2d/3d |
| `cfd-mesh: geometry/nurbs/curve.rs` | Partition of unity: `Σ N_{i,p}(t) = 1` |
| `cfd-mesh: csg/boolean.rs` | Volume identity: `vol(A)+vol(B)=vol(A∪B)+vol(A∩B)` |
| `cfd-mesh: welding/snap.rs` | Welding idempotency: `snap(snap(p)) = snap(p)` |
| `cfd-core: physics/cavitation/` | Rayleigh-Plesset bubble dynamics |
| `cfd-core: physics/hemolysis.rs` | Giersiepen–Wurzinger: HI = 3.62×10⁻⁷·τ^2.416·t^0.785 |
| `cfd-math: linear_solver/gmres/` | GMRES minimises ‖r‖ over Kₘ(A, r₀) Krylov subspace |
| `cfd-math: linear_solver/multigrid/amg.rs` | AMG O(1) convergence for M-matrices |
| `cfd-math: high_order/weno.rs` | WENO5: 5th-order smooth / 3rd-order near discontinuities |
| `cfd-math: high_order/spectral/` | Chebyshev exponential convergence O(e^{-cN}) for analytic u |
| `cfd-1d: resistance/models/` | Hagen-Poiseuille: Q = πD⁴ΔP/(128μL) |
| `cfd-1d: vascular/murrays_law.rs` | Murray's cube law: D₀³ = Σ Dᵢ³ |
| `cfd-2d: solvers/lbm/bgk.rs` | Chapman-Enskog: BGK D2Q9 → N-S at O(Ma²) |
| `cfd-2d: discretization/tvd.rs` | TVD stability: no new extrema introduced |
| `cfd-3d: fem/stabilization.rs` | SUPG/PSPG restores LBB stability for equal-order elements |
| `cfd-3d: vof/advection.rs` | PLIC: VOF volume conserved to machine precision |
| `cfd-3d: level_set/sussman.rs` | Sussman redistancing converges to exact SDF |
| `cfd-validation` | Richardson extrapolation; MMS order-of-accuracy; GCI = 1.25|ε|/(r^p−1) |
| `cfd-validation: analytical/` | Poiseuille, Couette, Womersley, Taylor-Green exact solutions |

---

## Feature Flags

| Flag | Default | Enables |
|------|---------|---------|
| `gpu` | ✅ on | `wgpu` GPU compute; turbulence kernels |
| `csg` | off | CSG Boolean ops in `cfd-mesh`; CSG examples |
| `mpi` | off | MPI domain decomposition, parallel I/O |
| `simd` | off | SIMD SpMV paths (benchmarked slower; use rayon instead) |
| `millifluidic` | off | Millifluidic-specific builders in `cfd-mesh` |

Most tests and CI run with `--no-default-features` to avoid the `wgpu` GPU overhead.

---

## Testing Strategy

| Layer | Tool | Scope |
|-------|------|-------|
| Unit (inline) | `#[test]` | Every public function; every data-structure invariant |
| Integration | `tests/` + nextest | Cross-crate pipelines; mesh build → validate → I/O round-trip |
| Property-based | `proptest` | Theorems (idempotency, volume identity, twin involution, …) |
| Debug invariants | `#[cfg(debug_assertions)]` | Mesh topology checks after every mutation |
| Fuzzing | `cargo-fuzz` | STL / JSON parsers (targets in `crates/cfd-mesh/fuzz/`) |
| Benchmarks | `criterion` | Solver throughput; SIMD vs scalar; scaling analysis |

**No mocks, stubs, or fixture-driven approximations.** Validate against closed-form mathematics.

---

## I/O Formats

| Format | Crate | Direction |
|--------|-------|-----------|
| STL ASCII / binary | `cfd-mesh` | read + write |
| VTK legacy ASCII | `cfd-mesh`, `cfd-io` | write |
| OpenFOAM polyMesh | `cfd-mesh` | write |
| HDF5 | `cfd-io` | read + write |
| CSV | `cfd-io` | read + write |
| cfd-schematics JSON | `cfd-mesh` (feature `scheme-io`) | read |
| Python (PyO3) | `cfd-python` | bindings |

OpenFOAM `PatchType` → OF type string mapping (in `io/openfoam.rs`):

```
Inlet    → "patch"    (physicalType: inlet)
Outlet   → "patch"    (physicalType: outlet)
Wall     → "wall"
Symmetry → "symmetry"
Periodic → "cyclicAMI"
```

---

## Key Files & Docs

| Path | Purpose |
|------|---------|
| `README.md` | Project overview; sprint history; build commands |
| `backlog.md` | Prioritised work items |
| `docs/adr.md` | Architecture Decision Record |
| `docs/srs.md` | System Requirements Specification |
| `docs/prd.md` | Product Requirements Document |
| `docs/checklist.md` | Current sprint tasks |
| `docs/gap_audit.md` | Gap analysis findings |
| `.agent/rules/persona.md` | Always-on agent persona and guidelines |
| `CHANGELOG.md` | Per-sprint change log |
| `ACCOMPLISHMENTS.md` | Completed milestones |

### Crate-Level Agent References

Each crate has its own `agents.md` with module structure, key APIs, theorems, and prohibited patterns:

| Path | Crate scope |
|------|-------------|
| `crates/cfd-core/agents.md` | Foundation: physics models, GPU/MPI/SIMD compute, plugin system |
| `crates/cfd-math/agents.md` | Numerical methods: GMRES/CG/AMG, spectral/DG/WENO, time steppers |
| `crates/cfd-mesh/agents.md` | Half-edge mesh, CSG, OpenFOAM I/O — phase plan Phases 1–10 |
| `crates/cfd-io/agents.md` | VTK/CSV/checkpoint/JSON I/O |
| `crates/cfd-schematics/agents.md` | Topology authority: NodeSpec/ChannelSpec, presets, visualisation |
| `crates/cfd-1d/agents.md` | Lumped-network Hagen-Poiseuille solver, resistance models, vascular |
| `crates/cfd-2d/agents.md` | 2D N-S: SIMPLE/PISO/LBM, turbulence zoo, Rhie-Chow |
| `crates/cfd-3d/agents.md` | 3D FEM/IBM/Level-Set/VOF/spectral, domain solvers |
| `crates/cfd-optim/agents.md` | SDT therapy + leukapheresis optimiser: 24 `DesignTopology` enum variants, GA, wave-channel SVG, optional mesh-export pipeline |
| `crates/cfd-validation/agents.md` | MMS, Richardson/GCI, benchmarks, conservation checks |
| `crates/cfd-python/agents.md` | PyO3 cfd-python bindings, all Python classes, build instructions |

---

## Sprint Workflow

```
Phase 1 (0-10%)  : Audit — README + codebase analysis → update backlog.md / checklist.md / gap_audit.md
Phase 2 (10-50%) : Execute — 50% audit, 50% atomic implementation (Spec → Test → Implement → Verify)
Phase 3 (50%+)   : Close — optimise, verify, sync all docs
```

Every sprint deliverable must pass all quality gates before the commit message is written.

---

## Commit Message Format

```
Phase N: <Short imperative summary>

- <Bullet describing each logical change>
- <Include test count and warning count>
- N tests pass, 0 warnings
```

---

## Prohibited Patterns

| Pattern | Reason |
|---------|--------|
| `unimplemented!()` / `todo!()` | Voids mathematical correctness guarantee |
| `#[allow(dead_code)]` on production items | Indicates unused code that should be removed |
| `clone()` in hot paths | Use GAT lending iterators or slice references instead |
| Epsilon-based Boolean classification | Use exact Shewchuk predicates |
| Parallel mesh types for the same representation | SSOT — one canonical implementation |
| External `ghost-cell` crate | `permission/` provides an owned, identical, zero-cost implementation |

---

*For crate-specific architecture, see each crate's `agents.md`.*
