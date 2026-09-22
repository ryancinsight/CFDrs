<!-- Compacted 2026-09-21 under the 1,000-line board budget: a board is a queue, not a ledger, so closed sections and closed item bodies are gone -- their record is the PR that closed them and its `Item:` trailer. Open items, anchors and live-marked residuals are kept. Recover any removed narrative with `git log -p -- <this file>`. -->
## Finding 2026-08-20: CFDrs scope-vs-delivery audit

### Measured baseline

| Quantity | Value | How measured |
|---|---|---|
| Workspace packages | 12 (11 libs + xtask) | `cargo metadata --offline --no-deps` |
| Library source | 1154 `.rs`, 268 175 lines under `crates/*/src` | `find`/`cat`/`wc -l` |
| All Rust incl. tests/benches/examples | 1393 `.rs`, 331 108 lines | same |
| Test functions | 3129 `#[test]` | `grep -rn '#\[test\]'` |
| `proptest!` blocks | 60 | `grep -rn 'proptest!'` |
| Cargo targets | 10 lib, 1 cdylib, 2 bin, 112 test, 49 example, 18 bench | `cargo metadata` |
| `todo!` / `unimplemented!` / TODO / FIXME | 0 / 0 / 0 | `grep` over `crates xtask tests benches examples` |
| Files > 500 lines | 135 | `find … -exec wc -l` (Atlas baseline records 134; +1 drift) |
| Crate-level `#![allow(...)]` in `lib.rs` | 292 across 10 of 11 crates | `grep -h '^#!\[allow' crates/*/src/lib.rs \| wc -l` |
| `#[allow(` sites | 97 | `grep -rn '#\[allow('` (Atlas baseline 89; **+8 regression**) |
| Crates with `#![deny(missing_docs)]` | **0 of 11** (all use `warn`) | `grep -l 'deny(missing_docs)' crates/*/src/lib.rs` |
| Production `.unwrap()` | 0 | `grep -rn '\.unwrap()' crates/*/src \| grep -v test` |
| `.expect(` sites in `crates/*/src` | 1701, of which **739 read `expect("expected value")`** | `grep -rho 'expect("[^"]*"' \| sort \| uniq -c` |
| Ignored tests | 32, of which 3 admit unimplemented behaviour | `grep -rn '#\[ignore'` |
| Book chapters | 19 numbered + 2 appendices, 1610 total lines incl. a 249-line glossary | `wc -l docs/book/*.md` |
| ADRs | 1 file (`docs/adr.md`, 1601 lines, 53 table rows), 0 numbered ADR files, no index | `ls docs/adr` → absent |

### Substrate consumption — clean at the manifest layer

### F-1 [arch][major] Ungated `#[global_allocator]` inside a library crate

### F-2 [verification][patch] 54 files / 10 543 lines sit outside every cargo target

### F-3 [arch][major] 20 563 lines of committed Python domain logic in `validation/`

### F-4 [arch][minor] Three in-repo SIMD implementations beside a 2-call-site hermes binding

- `crates/cfd-core/src/compute/simd/x86.rs` (286 lines) and `aarch64.rs` (197): hand-written `std::arch` intrinsics behind six `pub unsafe fn` (`advection_avx2`, `advection_sse41`, `diffusion_avx2`, `advection_neon`, `diffusion_neon`, `dot_product_neon`), each `#[target_feature(enable = …)]`.
- `crates/cfd-math/src/simd/` (972 lines: `cfd.rs`, `vector.rs`, `vectorization.rs`, `ops/`, `tests.rs`).
- `crates/cfd-2d/src/solvers/simd_kernels.rs` (562 lines).

### F-5 [correctness][patch] 739 content-free `expect` sites, including production solver paths

### F-6 [verification][patch] The flagship grid-convergence study asserts nothing

"pending migration of domain-specific multigrid code to leto-ops API".

### F-7 [verification][patch] The figure SSOT gate checks names, and its provenance is dangling

### F-8 [pm-hygiene][patch] Two parallel PM systems, a non-conformant ADR store, stale status

Its status marker reads "(in progress 2026-07-31; owner=Claude atlas session

### Conformance floor: the pedantic baseline is suppressed at crate level

### Completeness

## ATLAS-CFDRS-BACKWARD-STEP-108 — provider-owned geometry and shear (hosted closure pending 2026-08-19)

## Finding 2026-08-19: exact default-head Rust and figure gates pass

separate gates. The local standalone locked package path remains blocked by the

Atlas-locked compilation remains separately blocked by the shared overlay's

## ATLAS-CFDRS-BACKWARD-STEP-108 — field-derived reattachment (in progress 2026-08-17)

## CFDrs book and example gate — Apollo public-bound dependency

## Bounded Newton-Krylov recovery integration — CFDrs (2026-08-19)

- `crates/cfd-1d/src/solver/core/newton_fallback.rs` is now declared in the solver module graph and is called once when the existing bounded-amplitude stagnation detector classifies the Picard trajectory as stalled.
- The recovery budget derives its warm-up and Newton/Krylov limits from `SolverConfig.max_iterations`, so the fallback does not create a second unbounded solver attempt. `cfd-math::JfnkSolver::solve_checked` propagates residual callback failures as typed errors instead of panicking.
- Static source reachability, package formatting, and diff checks are clean. Locked local Cargo verification is blocked before compilation by the shared Atlas overlay/lock mismatch; provider hosted exact-head Rust and book-figure gates remain open until the source branch is checked.

## Lint-floor gate — partial closure (2026-08-06)

- The root workspace now owns the Atlas lint floor; every CFDrs workspace package and `xtask` inherit it.
- Passing evidence: `cargo clippy -p xtask --all-targets`, `cargo clippy -p cfd-core --lib`, `cargo nextest run -p cfd-core --lib` (246/246), `cargo test --doc -p cfd-core` (3/3), and `cargo run --manifest-path xtask/Cargo.toml -- legacy-migration-audit` (zero legacy dependencies, zero legacy source tokens, clean allowlist).
- Residual: the workspace all-target Clippy gate is not green because existing cfd-math unwrap/output sites, cfd-schematics missing docs, cfd-core test/bench lint debt, and unrelated format debt still need ratchet increments. This increment does not claim full workspace closure.

## cfd-math coarsening ordering slice — committed follow-up (2026-08-06)

- Multigrid coarsening no longer unwraps floating-point partial comparisons; finite values sort before unordered values, preserving deterministic diagnostics for NaN inputs. The regression test covers finite/NaN ordering.
- `cargo nextest run -p cfd-math --lib` passes 198/198. Focused cfd-math library Clippy remains red at 48 existing diagnostics after this slice.

## cfd-math hierarchy and storage invariants — follow-up (2026-08-06)

- Multigrid hierarchy/interpolation state now uses invariant-checked expectations instead of bare unwraps. JFNK and spectral operations use named C-contiguous storage helpers or explicit connectivity invariants.
- `cargo nextest run -p cfd-math --lib` passes 198/198. Focused cfd-math library Clippy remains red at 22 existing diagnostics after this slice.

## cfd-math diagnostics and DG output — closure slice (2026-08-06)

- Performance-monitor mutex accesses now carry invariant diagnostics and its calibration messages use structured tracing. DG progress, warnings, and completion metrics also use structured tracing instead of stdout writes.
- `cargo clippy -p cfd-math --lib` passes with the Atlas floor; cfd-math Nextest remains 198/198. Workspace all-target closure is still open on cfd-schematics docs, cfd-core test/bench lint debt, and format debt.

## cfd-schematics topology model documentation — partial closure (2026-08-07)

- The exported topology specification contract now documents its model types, fields, enum variants, aliases, and route lookup methods in `src/topology/model.rs`.
- `cargo clippy -p cfd-schematics --lib` still fails only on the existing documentation floor at this stage; diagnostics decrease from 712 to 611.
- Residual: 611 missing-documentation diagnostics remain across the package's other modules. No crate-wide lint suppression was added.

## cfd-schematics constants and config manifests — partial closure (2026-08-07)

- The `ConstantsRegistry` getters and fields, adaptive primitive constants, and public config module manifests now carry API documentation.
- `cargo clippy -p cfd-schematics --lib` still fails only on the existing documentation floor at this stage; diagnostics decrease from 611 to 534.
- Residual: 534 missing-documentation diagnostics remain across geometry, domain, interface, infrastructure, and topology modules outside this scope.

## cfd-schematics geometry builders — partial closure (2026-08-07)

- The public node and channel builder setters now document their domain effects, including names, geometry, visual roles, therapy zones, and Venturi metadata.
- The package documentation-floor residual decreases from 534 to 524 diagnostics. No runtime behavior or builder defaults changed.

## cfd-schematics geometry-generator entrypoints — partial closure (2026-08-07)

- The metadata configuration, generation entry points, and fluent builder now document their public contracts, including metadata, topology, lineage, and rendering inputs.
- The package documentation-floor residual decreases from 524 to 508 diagnostics. No runtime behavior or generation defaults changed.

## cfd-schematics linear geometry generators — partial closure (2026-08-07)

- The series and parallel geometry builders now document their specification- driven blueprint construction contracts.
- The package documentation-floor residual decreases from 508 to 506 diagnostics. No runtime behavior or generated geometry changed. The current working tree reports 492 because peer-owned documentation changes are also present in `domain/model/blueprint/analysis_impl.rs` and remain uncommitted.

## cfd-schematics selective-tree generator — partial closure (2026-08-07)

- The public selective-tree path specification, topology variants, request fields, and generator entrypoint now document their physical and topology contracts.
- The package documentation-floor residual decreases from 506 to 468 diagnostics. The current working tree reports 454 because peer-owned documentation changes remain uncommitted in `domain/model/blueprint/analysis_impl.rs`.

## Venturi geometry metadata metrics (CFDRS-AEQ-MET-47, 2026-08-05)

`2147749372`). A targeted retry was blocked by the shared target queue and

## Cross-consumer refresh and cfd-3d facade closure (2026-07-31)

compiles with `cargo check --tests`; its Nextest link step remains blocked by

## Schematic volume metric audit (2026-07-31)

## Legacy analytical benchmark consolidation (2026-07-31)

unresolved `leto-ops::{OwnedNumericLu, SymbolicLu, factor_symbolic}` imports

## Non-Newtonian analytical metric audit (2026-07-31)

peer-dirty `cfd-math::linear_solver::block_preconditioner` unresolved

## Blasius analytical metric audit (2026-07-31)

`cfd-math::linear_solver::block_preconditioner` unresolved

## Taylor-Green analytical metric audit (2026-07-31)

`cfd-math::linear_solver::block_preconditioner` unresolved

## Stokes analytical metric audit (2026-07-31)

`cfd-math::linear_solver::block_preconditioner` unresolved

## Couette and Poiseuille analytical metric audit (2026-07-31)

`cfd-math::linear_solver::block_preconditioner` unresolved

## Temperature-model metric refresh (2026-07-30, CFDRS-AEQ-MET-31)

## Microvascular blood metric refresh (2026-07-30, CFDRS-AEQ-MET-33)

## Ideal-gas metric refresh (2026-07-30, CFDRS-AEQ-MET-34)

## Non-Newtonian metric refresh (2026-07-30, CFDRS-AEQ-MET-32)

remains blocked by the unrelated peer `hephaestus-wgpu` bound error at

## Solid material metric refresh (2026-07-29, CFDRS-AEQ-MET-28)

## Eunomia complex compatibility refresh (2026-07-28)

## Blueprint cross-fidelity trace refresh (2026-07-28)

dependency path. The all-targets Clippy command remains blocked by 47

## 2D cell-tracking metric refresh (2026-07-29, CFDRS-AEQ-MET-27)

## Solver runtime refresh (2026-07-28)

## Aequitas public metric gap audit (2026-07-24)

### Verification refresh (2026-07-27, CFDRS-AEQ-MET-23)

### Verification refresh (2026-07-28, CFDRS-AEQ-MET-24)

### Verification refresh (2026-07-28, CFDRS-AEQ-MET-25)

### Verification refresh (2026-07-26)

### Verification refresh (2026-07-27)

also pass. The warning-denied library gate is currently blocked by a concurrent
that WIP produces unresolved imports before cfd-1d lint can complete. The

### Verification refresh (2026-07-27, CFDRS-AEQ-MET-18)

### Verification refresh (2026-07-27, CFDRS-AEQ-MET-19)

### Verification refresh (2026-07-27, CFDRS-AEQ-MET-20)

### Verification refresh (2026-07-27, CFDRS-AEQ-MET-21)

### Verification refresh (2026-07-27, CFDRS-AEQ-MET-22)

| ID | Evidence | Closure |
|---|---|---|
| `CFDRS-AEQ-MET-27` | `cfd-2d::solvers::cell_tracking` exposed physical positions, velocities, time, cell/material properties, and bifurcation geometry as raw scalars. | **IMPLEMENTED in this increment.** Public contracts use Aequitas `Length`, `Velocity`, `Time`, `MassDensity`, and `DynamicViscosity`; dimensionless routing uses `Dimensionless`; scalar extraction remains at interpolation, numerical, and Pries formula boundaries. Focused cfd-2d checks and cell-tracking Nextest provide the value-semantic gate. See [`cell-tracking-physical-metrics.md`](docs/atlas-migration/cell-tracking-physical-metrics.md). |
| `CFDRS-AEQ-MET-21` | Public transient composition and droplet simulation APIs accepted requested, calculated, and returned timepoints as `Vec<T>` after event/configuration time had been typed. | **IMPLEMENTED and package-verified.** Public timepoint vectors and timing accessors now use `Time<T>`; scalar conversion is private to sorting, tolerance, and solver boundaries. Test-target check passes; cfd-1d Nextest passes 736/736 with 3 skips; warning-denied all-target Clippy passes; no compatibility facade exists. Mixture fraction storage remains separately tracked as dimensionless representation. See [`transient-composition-metrics.md`](docs/atlas-migration/transient-composition-metrics.md). |
| `CFDRS-AEQ-MET-22` | Public transient mixture fraction storage, blood-hematocrit construction/accessors, weighted blends, tolerances, and node/edge concentration queries still exposed raw dimensionless scalars. | **IMPLEMENTED and package-verified.** Public maps, constructors, accessors, blends, and tolerances now use `Dimensionless<T>`; scalar extraction is limited to normalization, arithmetic, transport, and assertions. Test-target check passes; cfd-1d Nextest passes 736/736 with 3 skips; warning-denied all-target Clippy passes; doctests pass 8/8 with 3 ignored. No compatibility facade exists. Solver residuals remain equation-dependent. See [`transient-composition-metrics.md`](docs/atlas-migration/transient-composition-metrics.md). |
| `CFDRS-AEQ-MET-20` | Transient composition events, timing/control configuration, and snapshots exposed activation time, hematocrit, flow, pressure, CFL, and snapshot flow/time as raw scalars. | **IMPLEMENTED and focused-verified.** Public contracts now carry Aequitas `Time`, `Dimensionless`, `VolumetricFlowRate`, and `Pressure`; simulator conversion is confined to numerical boundaries. The cfd-1d test-target check passes; composition parity passes 21/21; droplet parity passes 9/9; literature validation passes 5/5; and the typed control regression proves value preservation. Solver residuals remain a separate equation-dependent classification; no compatibility facade exists. See [`transient-composition-metrics.md`](docs/atlas-migration/transient-composition-metrics.md). |
| `CFDRS-AEQ-MET-18` | `cfd-1d` `NodeProperties` and `NetworkMetadata` exposed pressure, temperature, total volume, and pressure/temperature ranges as raw scalar values after MET-17. | **IMPLEMENTED in this increment.** Public metadata contracts now carry Aequitas `Pressure`, `ThermodynamicTemperature`, and `Volume`; typed builder setters, default values, and range preservation are covered by value-semantic tests. `HashMap<String, T>` remains dimension-unknown metadata by contract. Locked `cfd-1d` check passes; Nextest passes 735/735 with 3 skips; doctests pass 8/8 with 3 ignored. Warning-denied Clippy reaches the peer-owned cfd-math missing-docs residual and is not a clean gate for this slice. |
| `CFDRS-AEQ-MET-19` | The transient droplet public boundary exposed droplet volume, injection/snapshot time, normalized positions, occupancy spans, and split thresholds as raw scalar values. | **IMPLEMENTED and focused-verified.** Public and internal droplet state now carries Aequitas `Volume`, `Time`, and `Dimensionless` quantities, with scalar extraction confined to transport formulas. The same revision passes `cargo check -p cfd-1d --tests --offline`; focused Nextest passes 9/9 droplet-parity tests and 5/5 literature-validation tests. Composition event/time contracts remain a separately tracked residual. See [`transient-droplet-metrics.md`](docs/atlas-migration/transient-droplet-metrics.md). |
| `CFDRS-AEQ-MET-17` | `cfd-1d` `Network` and `NetworkState` exposed pressure, volumetric-flow, and simulation-time state as raw scalar vectors; network analysis returned additional physical hydraulic metrics as raw scalars. Residual norms have model-dependent units. | **IMPLEMENTED in the current typed network-state increment; provider dimensions in `f19ba15`.** Public state and analysis contracts now use Aequitas quantities, and all in-tree callers are migrated without adapters. Residuals remain scalar under the documented equation-dependent classification. Locked library, test-target, and example checks pass; cfd-1d Nextest passes 731/731 with 3 skips in 23.458 s; doctests pass 8/8 with 3 ignored. Warning-denied library Clippy is pending the concurrent cfd-math module deletion being reconciled. See [`network-state-metrics.md`](docs/atlas-migration/network-state-metrics.md). |
| `CFDRS-AEQ-MET-16` | `cfd-1d` network `Edge` and `EdgeProperties` exposed flow rate, linear hydraulic resistance, quadratic hydraulic-loss coefficient, and parallel conductance as raw `T` values after geometry and vascular-result cutovers. `Network` pressure/flow/residual vectors were a separate solver-state boundary with additional unit semantics. | **VERIFIED in `a50a9e91`; provider dimensions in `f19ba15`.** Aequitas `VolumetricFlowRate`, `HydraulicResistance`, `QuadraticHydraulicResistance`, and `HydraulicConductance` now remain at public edge contracts; base scalars are extracted only inside resistance, validation, junction-loss, analysis, and matrix-assembly kernels. Locked `cargo check -p cfd-1d` passes; Nextest passes 731/731 with 3 skipped; cfd-1d doctests pass 8/8 with 3 ignored. The remaining solver-state boundary was audited and closed by MET-17. |
| `CFDRS-AEQ-MET-15` | The remaining vascular metric gap covered Murray optimal bifurcation geometry and Olufsen structured-tree terminal radius/impedance after the Womersley and network boundary migration. | **IMPLEMENTED in `e3b664e5`.** `OptimalBifurcation` now uses Aequitas `Length`, `Angle`, `VolumetricFlowRate`, `Dimensionless`, `DynamicViscosity`, and `Pressure`; `OlufsenParameters` uses `Length` and returns `HydraulicResistance`. Validation and PyO3 consumers convert at explicit boundaries. Provider aliases are from `446eb9f`. Locked `cfd-1d` check passes; Nextest passes 731/731 with 3 skipped; focused Womersley adversarial tests pass 2/2; locked `cfd-validation` check passes. See [`vascular-metrics.md`](docs/atlas-migration/vascular-metrics.md). |
| `CFDRS-AEQ-MET-14` | `cfd-1d` vascular `WomersleyNumber`, `WomersleyFlow`, `WomersleyProfile`, `VesselSegment`, `Bifurcation`, and `BifurcationNetwork` expose length, radius, pressure, density, viscosity, frequency, flow, and derived vascular results as raw scalar values. | **IMPLEMENTED.** Aequitas owns the named pressure-gradient, hydraulic-resistance, hydraulic-inertance, and compliance dimensions in `446eb9f`. The Womersley and bifurcation contracts now carry typed physical inputs/results and extract base scalars only at analytical kernels. Murray/Olufsen completion is recorded in `CFDRS-AEQ-MET-15`. |
| `CFDRS-AEQ-MET-13` | `cfd-1d` `ChannelType::Curved` and `Micromixer` exposed curvature radius, hydraulic diameter, and path length as raw SI scalars at public construction and storage boundaries. | **IMPLEMENTED.** Curved-channel radius and micromixer hydraulic diameter/path length now use Aequitas `Length`; scalar extraction remains at resistance and dynamic-parameter boundaries. In-tree constructors and value-semantic tests are migrated. `cargo check -p cfd-1d --tests` passes; Nextest passes 731/731 with 3 skipped; doctests pass 8/8 with 3 ignored; Rustdoc exits 0 with 11 pre-existing link warnings. Warning-denied Clippy remains blocked only by the pre-existing `cfd-math::matrix_zeros` dead-code warning, and no MET-13 test exceeds the committed runtime budget. See [`curved-micromixer-metrics.md`](docs/atlas-migration/curved-micromixer-metrics.md). Vascular and Womersley metrics remain outside this slice. |
| `CFDRS-AEQ-MET-12` | `cfd-1d` `CrossSection`, `ChannelGeometry`, `Edge`, and `EdgeProperties` exposed channel dimensions, area, hydraulic diameter, and length as raw SI scalars across public network construction and analysis. | **IMPLEMENTED.** Cross-section dimensions/custom area, channel length, edge area, and edge-property geometry now use Aequitas `Length` and `Area`. Scalar extraction remains at resistance, junction-loss, transient-transport, and analysis kernels; blueprint conversion, examples, and test fixtures use explicit adapters. The locked check passes and Nextest passes 729/729 with 3 skipped. Doctest, Rustdoc, Clippy, and runtime-budget limits remain recorded below. See [`channel-geometry-metrics.md`](docs/atlas-migration/channel-geometry-metrics.md). |
| CFDRS-AEQ-MET-11 | cfd-1d RectangularChannel, CircularChannel, PorousMembrane, OrganCompartment, and ChannelProperties exposed linear geometry, roughness, area, or volume as raw SI scalars; Component::volume returned Option<T>. | **IMPLEMENTED.** Component geometry now uses Aequitas Length, area methods return Area, Component::volume returns Volume, and ChannelProperties stores Length. Factory and setter scalar inputs convert immediately at their dynamic-parameter boundary; resistance models extract base scalars only for numerical formulas. cfd-1d check passes; Nextest passes 729/729 with 3 skipped; doctests pass 8/8 with 3 ignored; Rustdoc completes with 11 existing link warnings. Clippy remains blocked by the pre-existing cfd-math matrix_zeros dead-code warning. See docs/atlas-migration/component-geometry-metrics.md. |
| CFDRS-AEQ-MET-10 | cfd-1d SurfaceProperties stored roughness, contact angle, and surface energy as raw SI scalars. cfd-core WettingProperties, FluidSolidInterface, and InterfaceProperties likewise exposed surface tension and static/advancing/receding angles without Aequitas dimensions. | **IMPLEMENTED.** Public channel surface contracts now use Aequitas Length, Angle, and EnergyPerArea; material interfaces use SurfaceTension and Angle, and adhesion returns EnergyPerArea. Scalar extraction is confined to the Darcy resistance and cosine-law boundaries. cfd-core Nextest passes 259/259; cfd-1d Nextest passes 729/729 with 3 skipped; doctests pass 11 with 3 ignored; focused shape-factor validation passes 1/1. The non-Newtonian validation binary exceeds the committed 30-second budget, and Clippy remains blocked by the pre-existing cfd-math `matrix_zeros` warning. See docs/atlas-migration/surface-wetting-metrics.md. |
| `CFDRS-AEQ-MET-09` | `cfd-1d` cell-separation and kappa-aware cascade APIs exposed cell diameter, density, treatment/recovery hydraulic diameter, parent inflow velocity, Zweifach–Fung channel diameter, and optimization stage widths as raw SI scalars. The values cross public constructors, routing validation, and `cfd-optim` blueprint summaries. | **IMPLEMENTED in this increment.** `CellProperties`, `PeripheralRecovery`, `CascadeStage`, the public Zweifach–Fung functions, and `StageBlueprintSeparationSummary` now carry Aequitas `Length`, `MassDensity`, and `Velocity` values. Scalar extraction is confined to validation and numerical formula boundaries; all in-tree tests and blueprint construction are migrated. See [`cell-separation-physical-metrics.md`](docs/atlas-migration/cell-separation-physical-metrics.md). Touched-file rustfmt, scoped diff checks, residue scans, `cfd-1d` Nextest (728/728, 3 skipped), `cfd-optim` Nextest (137/137), and focused cell-separation validation Nextest (16/16) pass. The full `cfd-validation` package gate remains blocked by the committed 30-second budget: `test_venturi_flow_3d` and `microventuri_35um_case_produces_converged_informative_2d_result` timed out at approximately 31 seconds; those tests are outside this slice. |
| `CFDRS-AEQ-MET-23` | The remaining public cell-separation family exposed equilibrium position, residual force, Dean drag, direct fluid/geometry inputs, Fahraeus/CFL diameters and widths, viscosity/shear metrics, plasma-skimming diameters, and cross-junction geometry/flow inputs as raw SI scalars after MET-09. | **IMPLEMENTED; focused value-semantic gates pass.** Aequitas adds `Force`/`Newton`; the public cell-separation family now uses typed `Length`, `MassDensity`, `DynamicViscosity`, `Velocity`, `ReciprocalTime`, and `VolumetricFlowRate` contracts. Scalar extraction remains at validation and numerical formula boundaries; callers, tests, and cross-fidelity validation are migrated with no compatibility facade. Rustfmt, metadata, diff checks, residue scans, `cfd-1d` check, full cfd-1d Nextest (736/736, three skipped), focused cfd-validation Nextest (57/57), cfd-1d doctests (8/8, three ignored), and warning-denied cfd-1d Clippy pass. See [`cell-separation-force-metrics.md`](docs/atlas-migration/cell-separation-force-metrics.md). |
| `CFDRS-AEQ-MET-08` | `cfd-core` selective cavitation, `cfd-1d` Venturi screening, and `cfd-optim` Venturi placement/blueprint metrics discarded Aequitas types at public boundaries for pressure, density, velocity, length, viscosity, radius, and surface tension. | **FOCUSED VERIFIED.** The public physical contracts carry Aequitas quantities and all in-tree constructors are migrated. Scalar conversion remains only in numerical formula kernels and the documented serialized report DTO boundary. Provider support is in Aequitas commits `07e2252` and `6dc68c4`. The producer suite passes 1,127/1,127 with 3 skips; warning-denied all-targets Clippy and doctests pass for cfd-core/cfd-1d/cfd-optim. The broad cfd-3d/cfd-validation gate remains open on eight solver-heavy tests exceeding 30 seconds. See [`venturi-physical-metrics.md`](docs/atlas-migration/venturi-physical-metrics.md). |

| ID | Evidence | Closure |
|---|---|---|
| `CFDRS-AEQ-MET-07` | `cfd-1d/src/physics/hemolysis/mod.rs` exposed wall shear stress and exposure duration as raw `f64` arguments and fields, while the returned Giersiepen/Taskin indices were dimensionless. `cfd-1d` flow analysis, `cfd-optim` reporting, and `cfd-validation` passed those scalars directly. | **IMPLEMENTED and focused-verified.** Giersiepen and Taskin accept Aequitas `Pressure` and `Time`; `HemolysisExposure` stores the same typed inputs, all in-tree callers are migrated, and the formula owner remains cfd-core/local model code. The producer suite passes 1,127/1,127 with 3 skips. The broader cfd-validation gate retains the eight documented 30-second runtime timeouts, outside the typed hemolysis contract. See [`hemolysis-exposure-metrics.md`](docs/atlas-migration/hemolysis-exposure-metrics.md). |
| `CFDRS-AEQ-MET-06` | `cfd-3d::cascade` exposed channel geometry, flow rate, outlet pressure, wall shear, pressure drop, and maximum velocity as raw SI scalars. The inlet calculation already constructed Aequitas area, flow, and velocity internally, so the public boundary discarded the provider types. | **IMPLEMENTED and source-verified.** `CascadeChannelSpec`, `CascadeConfig3D`, `ChannelResult3D`, and `CascadeResult3D` carry Aequitas `Length`, `VolumetricFlowRate`, `Pressure`, and `Velocity`. Serde keeps the established SI scalar wire keys through explicit representation adapters; FEM and mesh code convert only at the scalar numerical boundary. Locked cfd-3d/cfd-validation check passes. The broader package suite remains runtime-blocked by eight solver-heavy tests at 30 seconds. See [`cascade-physical-metrics.md`](docs/atlas-migration/cascade-physical-metrics.md). |

### Verification refresh for CFDRS-AEQ-MET-07

- 2026-07-22 (resolved in CFD-BOOK-CLOSEOUT-1): the stale book commit expanded source-backed chapter indexes with public types and behavioral contracts that do not exist in CFDrs, and its SUMMARY linked directly to `../../../parity_artefacts/INDEX.md`. mdBook treated that path as a source page and overwrote tracked archive HTML outside the book on every build. The fix-forward retains the twelve non-duplicated pages backed by real examples, restores the expanded chapters to their prior source-grounded content, consolidates linear-algebra parity on Leto Ops' analytical oracle, and routes archive navigation through a local book page. Exact scans find no Rust definitions for the rejected contracts (`NonDimParams`, `BoundaryKind`, `GhiaOracle`, `ShearReport`, `ScreeningConfig`, `GiersiepenWurzinger`, `HemolysisPath`, or `ParetoFront`). The parity HTML blob is byte-identical before and after mdBook (`85af4889c39f6d03d78b0dfceeb217f5d260efb5`). Evidence tier: source-definition audit, cumulative-diff review, successful book rebuild, documented example compilation, warning-denied Clippy, configured Nextest 177/177, and doctests 16/16.

- 2026-07-22 (resolved in CFD-SCHEMATIC-PATH-1): the stale `codex/cfd-example-paths` branch contained a valid native-path boundary that never reached main. Current main still required ten lossy or fallible UTF-8 conversions around renderer calls. The recovery ports only that contract onto the current tree: renderer traits borrow `Path`, plotting facades accept `AsRef<Path>`, sidecar naming stays in `OsStr`, and every live caller passes its native path directly. Evidence tier: exact branch/content comparison, affected package/example compilation, warning-denied Clippy, and a clean renderer conversion scan. Configured Nextest passes all 177 `cfd-schematics` tests, including native non-UTF-8 path format detection.

- 2026-07-21 (resolved in CFD-IRIS-COLOR-1): `cfd-schematics` duplicated Iris's normalized color-law role with a consumer-owned enum and three local formulas. Each edge or node lookup also rebuilt a value vector and rescanned its full map. The duplicate laws and wrapper enum are deleted; callers use Iris `NamedColorMap` directly. `AnalysisOverlay` now lends or owns maps with `Cow`, rejects non-finite values at construction, and stores one finite range per map. The old render cost was `Theta(E^2 + V^2)` range work with `E + V` transient allocations and `Theta(max(E, V))` transient elements; the new cost is `Theta(E + V)` construction, zero transient range allocations, and expected `O(1)` map/color lookup. This is an asymptotic and allocation proof, not a measured speedup claim. Evidence tier: source-level ownership and residue audit; focused value-semantic Nextest 176/176; warning-denied all-target/all-feature Clippy; affected example compilation; 16 passing doctests; warning-denied Rustdoc; and an executed, visually inspected Venturi pressure overlay. Major SemVer classification was attempted but its isolated temporary graph cannot build cfd-core because existing CFDrs direct pins and Proteus/Hephaestus transitive pins select distinct Aequitas and Leto source identities. That provider-pin coherence gap is independent of Iris. Kwavers volume rendering remains a separately claimed consumer migration.

- 2026-07-20 (resolved in CFD-LAPLACIAN-PROVIDER-1): cfd-math directly implemented the two-dimensional CPU Laplacian and cfd-core carried another copy as a GPU test oracle, although Hephaestus already owned the WGPU stencil. The CPU solver evaluated `-∇²` while the GPU solver evaluated `∇²`. Leto now owns the validated spacing, boundary, polarity, and native-precision CPU operation; Hephaestus consumes that contract; both CFD solver operators select negative polarity. The local formulas are deleted. Evidence tier: provider type unification; exact full-grid CPU and real-WGPU regressions; configured Nextest 622/622; all-feature and CPU-only checks; warning-denied Clippy/Rustdoc; six runnable doctests; and the updated example check. Three-dimensional, variable-coefficient, and SIMD diffusion operators remain separate contracts outside this slice. `cargo-semver-checks` was attempted but blocked in nightly Rustdoc by a long-lived shared-target Leto IDE check; the public constructor break remains classified `[major]`.

- 2026-07-17 (resolved stale work; upstream gap open): removing rsparse by routing `DirectSparseSolver` through unpreconditioned GMRES is not a valid provider migration. The solver chain already attempts GMRES after its exact sparse-LU tier, while cfd-2d invokes the direct tier specifically after GMRES stagnation or breakdown; the substitution therefore destroys failure-mode independence and contradicts the public direct-solver contract. Leto 0.38 exposes sparse CG and GMRES but no sparse direct factorization. CFDrs retains rsparse until upstream item `LETO-SPARSE-DIRECT-1` provides a generic sparse direct API and differential conformance. Evidence tier: source-level dependency/call-graph inspection, exact tree equivalence to `main`, focused value-semantic Nextest (4/4 cfd-math and 1/1 cfd-2d), and warning-denied cfd-math Clippy.

- 2026-07-17 (resolved): `GpuContext` acquires and queries through Hephaestus's `ComputeDeviceAcquisition` and `ComputeDeviceCapabilities` seams after provider release 0.16.1 repaired typed downlevel acquisition. The derived seven-storage-binding request preserves the full downlevel descriptor; raw adapter and feature methods are deleted. Nextest serializes only provider-acquiring tests through `gpu-device`, eliminating process-level WGPU device races while CPU tests remain concurrent. Evidence tier: compile-time API removal and empty source scan; value-semantic typed-limit regression; cfd-core GPU 245/245, cfd-math GPU 362/362, cfd-2d GPU 570/570 (27 pre-existing skips), root integration 26/26; warning-denied touched targets; doctest/rustdoc; and SemVer's expected major-only classification. The root all-target example lint baseline is independently tracked by CFD-EXAMPLE-CLIPPY-1.

- 2026-07-17 (resolved): root `cfd-suite --all-targets` Clippy no longer reports the 29 diagnostics formerly distributed across seven validation examples. Four retained examples now execute provider-owned cfd-1d/cfd-2d calculations; three unreferenced static reports are deleted rather than presenting hardcoded validation outputs. Evidence tier: executable examples plus warning-denied all-target Clippy.

- 2026-07-17 (open): `BifurcationSolver3D` builds an unlabeled SDF volume mesh but integrates daughter flow only across `outlet_0` and `outlet_1` labels. The resulting zero daughter flows are deterministic, so the invalid root FEM example is deleted. CFD-3D-BIFURCATION-BOUNDARIES-1 owns the upstream mesh terminal-facet contract and cfd-3d flow regressions. Evidence tier: direct executable reproduction and source inspection of mesh construction and label-based integration.

- 2026-07-17: `cfd-core::compute::gpu::GpuContext::synchronize` now delegates completion to Hephaestus `ComputeDevice::synchronize`; `GpuContext` no longer exposes raw WGPU device, queue, or limit fields; and cfd-2d creates its Poisson solver through `GpuPoissonSolver::from_context`. Evidence tier: compile-time provider integration plus GPU-enabled value-semantic regression coverage (244/244 cfd-core and 2/2 accelerated cfd-2d nextest), warning-denied cfd-core all-target Clippy, and exact source audits with no `device.poll`, old Poisson constructor, context device/queue access, or public raw-buffer accessor. The remaining adapter/feature introspection risk is resolved by the 0.3.0 typed capability boundary above.

- **Verification closure (2026-07-17)**: `cargo doc -p cfd-core --no-deps --features gpu --locked` completes warning-clean after the final raw-buffer visibility change. The final source state is verified by cfd-core GPU nextest 244/244, cfd-2d accelerated nextest 2/2, warning-denied cfd-core/ cfd-2d Clippy, and the package documentation gate.

- **SemVer classification (2026-07-17)**: Git-baseline semver checks identify the intentional removals as three breaking API classes under a minor-change assumption. The explicit major-change classification passes. CFDrs remains pre-1.0, so the workspace advances from `0.1.0` to `0.2.0` and records the migration in the `0.2.0` changelog section without retaining a compatibility surface.

- 2026-07-17: Preserved the stale peer's valid Leto source revision `6aedde0c7835238867d6f3cd17b030f7e69cb6f2`, which is merged on Leto `main`, and advanced its Moirai companion pin to merged `main` `5ead788c70c728d971237d7afa0b915ea7cf87e3`. Locked metadata, all-feature `cfd-schematics` check, and warning-denied Clippy pass; evidence tier is compile-time integration.

- 2026-07-16: Updated the workspace Moirai source pin to merged `main` `5ead788c70c728d971237d7afa0b915ea7cf87e3`. Locked metadata resolves Moirai 0.4 and Themis 0.10; `cfd-schematics` all-feature check, warning-denied Clippy, focused nextest, doctests, and docs pass. The source- and test-level evidence is compile-time integration plus value-semantic test coverage.
- 2026-07-16: Removed the `cfd-schematics` strict-Clippy baseline in the touched test/example cone. Direct geometry/phase values use exact bit-pattern assertions, and the Venturi example exposes named physical fields instead of positional tuple entries. The workspace-wide formatter remains blocked by unrelated pre-existing formatting in `crates/cfd-schematics/src/error.rs`; touched files pass `rustfmt`.

- 2026-07-10: Removed stale allowlist entries for `cfd-1d/src/scalar.rs` and `cfd-3d/src/scalar.rs`. Both seams are provider-native and contain no legacy dependency tokens; the active `cfd-core` compute-dispatch diff is untouched.

## Sprint 2026-07-07: cfd-1d/cfd-3d Eunomia identity seam

- **Resolved direct scalar dependency**: `Cargo.toml`, `crates/cfd-1d/Cargo.toml`, and `crates/cfd-3d/Cargo.toml` no longer declare direct `num-traits` dependencies for the 1D/3D solver scalar seams.
- **Resolved identity ownership**: `Cfd1dScalar` and `Cfd3dScalar` now expose `zero()` and `one()` through the Eunomia `NumericElement` constants already required by the crate-local scalar contracts, removing `num_traits::{Zero,One}` as a supertrait requirement.
- **Evidence tier**: compile-time integration, empirical nextest coverage, touched-file formatting, and static source/manifest audit. In `D:/atlas/repos/CFDrs`, touched-file rustfmt passed; `rustup run nightly cargo check -p cfd-1d -p cfd-3d` passed; direct residue scan found no `num_traits` or direct `num-traits` hits in the touched manifests/scalar cones; and `rustup run nightly cargo nextest run -p cfd-1d -p cfd-3d --status-level fail` passed 1122/1122 with one existing slow 3D mesh-convergence validation.
- **Residual risk**: package-wide fmt and all-targets clippy remain blocked by pre-existing unrelated formatting/lint debt outside this slice. Lockfile `num-traits` entries, if present, are transitive provider dependencies owned by upstream crates rather than direct CFDrs scalar-seam dependencies.

---

## Sprint 2026-07-05: public sparse/linear-solver Leto boundary

- **Resolved public sparse storage boundary**: `crates/cfd-math/src/sparse/mod.rs` now exposes `leto_ops::CsrMatrix<T>` as the public `SparseMatrix<T>` alias. Sparse builders, assembly helpers, sparse operations, and sparse tests construct and operate on Leto CSR directly rather than converting through nalgebra-sparse.
- **Resolved public solver vector/matrix boundary**: `LinearOperator`, `Preconditioner`, `LinearSolver`, direct solver, solver-chain, and migrated solver fixtures use Leto CSR and Leto `Array1<T>` boundaries in the requested cone. `cfd-validation::numerical` now stores computed and analytical validation vectors as `leto::Array1<T>` and uses Leto CSR for linear-solver validation test cases.
- **Evidence tier**: compile-time provider integration, package empirical regression tests, clippy, rustdoc, and static source audit. In `D:/atlas/repos/CFDrs`, cfd-math check passed; cfd-math test-target check passed; cfd-math all-target clippy passed; cfd-math doc passed; cfd-math nextest passed 361/361; cfd-validation check passed; cfd-validation all-target clippy passed; cfd-validation doc passed after fixing a stale intra-doc link; the targeted residue scan found no `nalgebra_sparse::CsrMatrix`, public nalgebra-sparse re-export, `DVector`, `row_offsets()`, `try_from_csr_data`, `CooMatrix`, or `nalgebra::` matches under the migrated sparse/linear-solver/validation files.
- **Residual risk**: The requested sparse/linear-solver public boundary has no known nalgebra sparse/vector holdouts in the scanned cone. Full `cfd-validation` package nextest remains blocked by the existing venturi cross-fidelity convergence tests `option2_selected_45um_geometry_routes_to_fallback_and_converges` and `microventuri_35um_case_produces_converged_informative_2d_result`, which are outside this boundary. Broader CFDrs provider migration still has nalgebra residue in other crates and contexts outside this slice.

---

## Sprint 2026-07-04: Solver Chain and FEM Consumer Leto Vector Boundary

- **Resolved chain vector API**: `crates/cfd-math/src/linear_solver/chain.rs` now exposes `LinearSolverChain::solve` and `solve_with_guess` with `leto::Array1<T>` RHS/result vectors instead of nalgebra `DVector`.
- **Resolved 2D direct-fallback consumers**: `crates/cfd-2d/src/linear_solver_bridge.rs` is the single cfd-2d conversion boundary from nalgebra work vectors into the Leto-backed `DirectSparseSolver`; momentum and pressure fallback paths, plus the momentum regression test, route through it.
- **Resolved 3D FEM assembly consumers**: `crates/cfd-3d/src/fem/leto_bridge.rs` is the single FEM conversion boundary for nalgebra work vectors crossing `SparseMatrixBuilder::build_with_rhs` and `LinearSolverChain`. `FemSolver` and `ProjectionSolver` now use the crate-level `Cfd3dScalar` seam, which carries the Leto real-scalar provider bound.
- **Evidence tier**: compile-time provider integration, focused empirical nextest, clippy on touched library surfaces, static residue scan, and diff hygiene. In `D:/atlas/repos/CFDrs`, `rustup run nightly cargo fmt -p cfd-math -p cfd-1d -p cfd-2d -p cfd-3d --check`, `cargo check -p cfd-math --no-default-features --lib`, `cargo check -p cfd-1d --no-default-features --lib`, `cargo check -p cfd-2d --no-default-features --lib`, `cargo check -p cfd-3d --no-default-features --lib`, `cargo nextest run -p cfd-math --no-default-features chain direct_solver core_solver simple_gmres --status-level fail` (4/4), `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`, `cargo clippy -p cfd-2d --no-default-features --lib -- -D warnings`, and `cargo clippy -p cfd-3d --no-default-features --lib -- -D warnings` passed.
- **Residual risk**: `cargo clippy -p cfd-2d --no-default-features --all-targets -- -D warnings` is blocked in `cfd-validation`, not cfd-2d: validation benchmark modules still pass nalgebra `DVector` to public `cfd_math::sparse::spmv`, and generic validation/1D literature paths need the Leto scalar bound propagated after the cfd-1d network solver seam moved. The broader iterative solver/preconditioner traits still expose nalgebra `DVector` until that trait family moves to Leto arrays.

---

## Sprint 2026-07-04: cfd-1d Eunomia/Leto Scalar Boundary

- **Resolved**: `cfd-1d` no longer exposes scattered nalgebra `RealField` imports as its domain scalar contract. The crate now routes domain, network, component, resistance, vascular, solver, transient, and analysis generic bounds through `Cfd1dScalar`, which explicitly combines the remaining nalgebra linear-system backend requirement with the Eunomia scalar provider contract consumed by migrated `cfd-core` APIs.
- **Geometry contract corrected**: `NetworkDomain::contains_1d` now accepts Leto `Point1<T>`, matching the migrated `cfd-core::geometry::Domain` contract.
- **Evidence tier**: static source audit, compile-time integration, full empirical nextest, and scoped clippy. `rustup run nightly cargo fmt -p cfd-1d --check`, `cargo check -p cfd-1d --no-default-features --lib`, `cargo nextest run -p cfd-1d --no-default-features --status-level fail` (725/725, 3 skipped), and `cargo clippy -p cfd-1d --no-default-features --lib -- -D warnings` passed. A direct scan found nalgebra `RealField` only in `crates/cfd-1d/src/scalar.rs`, where it documents the remaining matrix backend boundary.
- **Residual risk**: `cargo check -p cfd-2d --no-default-features --features gpu --lib` now reaches `cfd-2d` and fails on `cfd-2d`'s own nalgebra `RealField` bounds around migrated `cfd-core` boundary/fluid APIs. `cargo clippy -p cfd-1d --no-default-features --all-targets -- -D warnings` remains blocked by pre-existing lint debt in tests, examples, and cell-separation modules outside this scalar-boundary slice.

---

## Sprint 2026-07-04: cfd-core GPU Poisson Hephaestus Kernels

- **Resolved**: `crates/cfd-core/src/compute/gpu/poisson_solver.rs` no longer owns raw WGPU compute pipelines, bind-group layouts, parameter buffers, staging buffers, manual `map_async` readback, or `futures`/`mpsc` mapping channels. Jacobi, red-black, and residual entry points now dispatch through Hephaestus `WgslMultiStorageKernel`.
- **Provider boundary tightened**: Poisson field/source/residual storage now uses `WgpuDevice`'s `ComputeDevice` upload, allocation, and download contracts. The public constructor still accepts the existing WGPU device and queue handles, then immediately wraps them in a Hephaestus provider.
- **Shape contract corrected**: The solver now stores `nx`, `ny`, `dx`, and `dy` from construction and validates `phi.len() == source.len() == nx * ny` before dispatch. The previous implementation inferred a square grid from `phi.len()`, ignoring the constructor geometry.
- **Evidence tier**: static source audit, compile-time integration, full empirical nextest, and scoped clippy. `rustup run nightly cargo fmt -p cfd-core --check`, `cargo check -p cfd-core --features gpu`, `cargo check -p cfd-core --no-default-features`, `cargo clippy -p cfd-core --features gpu --all-targets -- -D warnings`, and full `cargo nextest run -p cfd-core --features gpu --status-level fail` (231/231) passed. A direct scan of `poisson_solver.rs` found no `create_buffer_init`, `create_buffer(`, `params_buffer`, `ComputePipeline`, `BindGroupLayout`, `map_async`, `futures::channel`, `poll(wgpu::PollType`, or `std::sync::mpsc` residue.
- **Residual risk**: `cfd-2d` accelerated Poisson consumer verification now reaches `cfd-2d` and remains blocked by `cfd-2d` Eunomia/nalgebra trait-bound errors before the consumer reaches the GPU Poisson path. Broader CFDrs GPU cleanup still has raw WGPU orchestration in non-Poisson kernels and tests.

---

## Sprint 2026-07-04: cfd-validation Direct Vector2 Provider Closure
- **Resolved**: `crates/cfd-validation/src/manufactured/navier_stokes.rs`, `crates/cfd-validation/src/conservation/{momentum.rs,angular_momentum.rs,mod.rs}`, `crates/cfd-validation/src/benchmarks/vorticity_stream.rs`, and MMS/conservation validation tests no longer import or use nalgebra `Vector2`. These surfaces now use `leto::geometry::Vector2`; component access uses Leto `[0]`/`[1]` indexing.
- **Boundary**: This closes direct nalgebra `Vector2` ownership in cfd-validation source and tests. It does not remove nalgebra from cfd-validation because dense `DMatrix`, 3D `Vector3`, and geometry point/vector surfaces remain separate Leto/Gaia provider migrations.
- **Evidence tier**: static source audit, compile-time integration, focused empirical nextest, and scoped clippy. `cargo fmt -p cfd-validation --check` passed. `cargo check -p cfd-validation --no-default-features --tests` passed. `cargo clippy -p cfd-validation --no-default-features --lib -- -D warnings` passed. Focused `cargo nextest run -p cfd-validation --no-default-features manufactured mms conservation taylor momentum angular --status-level fail` passed 104/104 tests. A direct scan over `crates/cfd-validation/{src,tests}` found no `nalgebra::Vector2`, `use nalgebra::Vector2`, or `use nalgebra::{..., Vector2}` residue.
- **Residual risk**: Full `cargo clippy -p cfd-validation --no-default-features --lib --tests -- -D warnings` remains blocked by pre-existing unrelated lints in cross-fidelity, schematics, and other validation tests outside this provider slice. Remaining provider work should move cfd-validation dense matrices to Leto storage and geometry types to Gaia/Leto before nalgebra can be removed from its manifest.

---

## Sprint 2026-07-04: cfd-2d Direct Vector2 Provider Closure
- **Resolved**: `crates/cfd-2d/src/physics/vorticity_stream.rs`, `crates/cfd-2d/src/piso_algorithm/corrector.rs`, `crates/cfd-2d/src/solvers/lbm/solver.rs`, `crates/cfd-2d/src/physics/immersed_boundary.rs`, `crates/cfd-2d/examples/blood_venturi.rs`, and `crates/cfd-2d/benches/solver_benchmarks.rs` no longer import or use nalgebra `Vector2`. These surfaces now use `leto::geometry::Vector2`; component access uses Leto `[0]`/`[1]` indexing. `crates/cfd-validation/src/benchmarks/vorticity_stream.rs` was updated as the downstream public-API consumer of `VorticityStreamSolver::velocity_field`.
- **Boundary**: This closes direct nalgebra `Vector2` ownership for cfd-2d source/tests/examples/benches. It does not remove nalgebra from cfd-2d because scalar `RealField`, boundary `Vector3`, dense `DVector`/`DMatrix`, and nalgebra-sparse matrix surfaces remain separate provider migrations. `physics::immersed_boundary` still uses nalgebra `DMatrix` for force/velocity matrices pending a Leto dense-storage migration.
- **Evidence tier**: static source audit, compile-time integration, focused empirical nextest, and scoped clippy. `cargo fmt -p cfd-2d -p cfd-validation --check` passed. `cargo check -p cfd-2d --no-default-features --examples --benches` passed. `cargo check -p cfd-validation --no-default-features` passed. `cargo clippy -p cfd-2d --no-default-features --example blood_venturi --bench solver_benchmarks -- -D warnings` passed. `cargo clippy -p cfd-validation --no-default-features --lib -- -D warnings` passed. Focused `cargo nextest run -p cfd-2d --no-default-features vorticity corrector lbm immersed --status-level fail` passed 44/44 tests. A direct scan over cfd-2d source/tests/examples/benches plus the touched cfd-validation benchmark found no `nalgebra::Vector2`, `use nalgebra::Vector2`, or `use nalgebra::{..., Vector2}` residue.
- **Residual risk**: Full `cargo clippy -p cfd-2d --no-default-features --all-targets -- -D warnings` remains blocked by pre-existing unrelated lints in examples/tests/modules outside this provider slice. The provider migration still needs cfd-math Leto linear-solver storage, cfd-core boundary vector migration, and Eunomia scalar contract replacement before cfd-2d can remove direct nalgebra/nalgebra-sparse manifest ownership.

---

## Sprint 2026-07-04: cfd-2d Problem/Streamtube Atlas Provider Seam
- **Resolved**: `crates/cfd-2d/src/problem.rs` now stores incompressible problem and solution velocity fields with `leto::geometry::Vector2` and routes initial pressure, velocity-magnitude maxima, and pressure maxima through `crates/cfd-2d/src/scalar.rs`/Eunomia instead of local nalgebra vector storage or direct `T::zero()` folds. `crates/cfd-2d/src/physics/streamtube/partitioning.rs` no longer uses direct `num_traits::{Float,FromPrimitive}`, `T::from_f64(...).unwrap()`, `T::zero()`, `T::one()`, `Float::abs`, `Float::sqrt`, or scalar `.abs()` in the touched APIs and tests; constants, absolute values, and square roots now route through `eunomia::{FloatElement,NumericElement}` and the crate-local scalar adapter.
- **Boundary**: This slice is limited to the problem setup and streamtube partitioning scalar/vector-provider seam. `problem.rs` still carries `nalgebra::RealField` because `cfd_core::physics::{boundary,fluid}` types are still nalgebra-bound upstream; removing that bound requires an upstream cfd-core provider migration. It does not remove the cfd-2d manifest's direct `num-traits` dependency because direct residues remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed. `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo nextest run -p cfd-2d --no-default-features problem streamtube separating --status-level fail` passed 4/4 tests. `git diff --check` passed for the touched problem/streamtube and PM artifact files. A direct-provider scan over both touched files found no `num_traits`, `FromPrimitive`, `Float::`, `T::from_*`, `.to_f64()`, `T::zero()`, `T::one()`, scalar `.abs()`, or local nalgebra `Vector2` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d immersed-boundary tests, momentum setup/interpolation/boundary, turbulence validation, and f64-only/test scalar surfaces; full cfd-2d direct `num-traits` removal remains a larger crate-level Eunomia migration before the manifest dependency can be dropped. Full `problem.rs` nalgebra-bound removal is blocked by upstream `cfd-core` boundary/fluid contracts.

---

## Sprint 2026-07-04: cfd-3d Level-Set Eunomia Scalar Seam
- **Resolved**: `crates/cfd-3d/src/level_set/{weno,advection,solver}.rs` no longer import or bound direct `num_traits::{FromPrimitive,Float}`. WENO5-Z weights, SSPRK3 coefficients, transport input validation, narrow-band limits, and reinitialization/Godunov math now route through `cfd-3d::scalar`, backed by Eunomia `FloatElement` and `NumericElement`.
- **Boundary**: This is the level-set module scalar-provider cleanup. Direct cfd-3d `num-traits` ownership is now closed by the later root lib-test/manifest slice. This entry still preserves the current `nalgebra::Vector3` boundary pending the larger Leto/Gaia geometry migration.
- **Evidence tier**: static source audit, compile-time integration, focused empirical nextest, and lib clippy. `cargo fmt -p cfd-3d --check` passed. `cargo check -p cfd-3d --no-default-features` passed. Focused `cargo nextest run -p cfd-3d --no-default-features level_set --status-level fail` passed 13/13 tests. `cargo clippy -p cfd-3d --no-default-features --lib -- -D warnings` passed. A targeted scan over `crates/cfd-3d/src/level_set` and `crates/cfd-3d/src/scalar.rs` found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `num_traits::Float`, `T::zero()`, `T::one()`, or `Float::` residue.
- **Residual risk**: `cargo clippy -p cfd-3d --no-default-features --all-targets -- -D warnings` is still blocked by pre-existing lint debt in unrelated cfd-3d test/module code (`poiseuille_test`, `fem_tests`, `smagorinsky_test`, `blueprint_integration`, `vof_tests`, `robustness_tests`, `bifurcation`, `venturi`, `trifurcation`, and VOF modules). Full provider completion still requires the remaining cfd-3d direct `num-traits` cleanup plus Leto/Gaia replacement of nalgebra geometry and storage surfaces.

---

## Sprint 2026-07-04: cfd-1d Direct num-traits Removal and Resistance/Vascular Eunomia Cleanup
- **Resolved**: `cfd-1d` no longer declares or directly references `num-traits`. The resistance scalar contract now uses Eunomia `FloatElement`/`NumericElement` for scalar construction, diagnostics, and transcendental/math operations. The touched hydraulic resistance, serpentine, slug-flow, Bessel/Womersley, structured-tree, bifurcation, network blueprint/sink, solver-analysis, and package-test seams no longer import `num_traits`, use `FromPrimitive`/`ToPrimitive`, call `T::from_f64`/`T::from_usize`/`T::from_u32`, or bridge through `nalgebra::try_convert`.
- **Boundary**: This closes direct `num-traits` ownership for the `cfd-1d` crate. It does not remove transitive `num-traits` through `approx`, `nalgebra`, `nalgebra-sparse`, `half`/Eunomia/Leto/Hephaestus, `num-complex`, or other provider stacks. It also does not replace the remaining nalgebra/nalgebra-sparse storage and solve boundaries; those stay as Leto-backed dense/sparse migration work.
- **Evidence tier**: static source audit, compile-time integration, empirical nextest, and dependency-tree audit. `cargo fmt -p cfd-1d --check` passed. `cargo check -p cfd-1d` passed. `cargo nextest run -p cfd-1d` passed 725/725 tests with 3 skipped. A full scan over `crates/cfd-1d/Cargo.toml`, `crates/cfd-1d/src`, and `crates/cfd-1d/tests` found no direct `num_traits`, `num-traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`, `T::from_usize`, `T::from_u32`, `nalgebra::try_convert`, or `.to_f64().unwrap` residue.
- **Residual risk**: `cargo clippy -p cfd-1d --all-targets -- -D warnings` remains blocked by existing all-target lint debt outside this provider cleanup, including `blueprint_solve_trace.rs`, `adversarial_tests.rs`, `resistance_model_validation.rs`, `medical_millifluidic_screening.rs`, `geometry_integration_demo.rs`, cell-separation tests/modules, droplet-regime tests, entrance-model tests, matrix-assembly tests, and venturi coefficient tests. Broader Atlas migration work remains for Leto/nalgebra-sparse storage replacement and Hephaestus higher-level GPU kernel ownership.

---

## Sprint 2026-07-04: cfd-1d Domain Components Eunomia Boundary Cleanup
- **Resolved**: `cfd-1d` domain components no longer carry direct `num_traits` conversion/math bounds. The component trait pressure-drop calculation uses Eunomia `NumericElement::abs`; factory defaults and component constants use the existing Atlas provider conversion seam; and channel, membrane, mixer, pump, valve, and sensor implementations no longer import `FromPrimitive` or `Float`.
- **Boundary**: This slice covers `crates/cfd-1d/src/domain/components/{mod,channels,factory,membranes,mixers,pumps,sensors,valves}.rs`. It preserves the current nalgebra `RealField` and resistance-model boundaries for later Leto/Eunomia work.
- **Evidence tier**: compile-time integration, empirical nextest, static source audit, and lint regression audit for the touched file. `cargo fmt -p cfd-1d --check` passed. `cargo check -p cfd-1d` passed. `cargo nextest run -p cfd-1d` passed 725/725 tests with 3 skipped. A focused component scan found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`, `T::from_usize`, `Float::`, or generic `.abs()` residue. `cargo clippy -p cfd-1d --all-targets -- -D warnings` no longer reports the touched `channels.rs` item-order lint.
- **Residual risk**: Full `cfd-1d` all-target clippy remains blocked by existing lint debt outside this slice in tests/examples, domain-network, cell-separation, vascular, solver-core, and benches. Broader `cfd-1d` direct-provider residue remains in solver-core, domain-network, vascular Bessel/Womersley, tests, and remaining nalgebra/nalgebra-sparse storage boundaries.

---

## Sprint 2026-07-04: cfd-1d Channel/Branching/Analysis Eunomia Boundary Cleanup
- **Resolved**: The next coherent `cfd-1d` provider seam no longer depends on direct `num_traits` bounds. Channel flow-regime classification, channel flow-resistance constants and powers, channel geometry perimeter math, Poiseuille shape factors, branching network solver bounds, and network pressure/flow/resistance/performance analysis paths now route scalar construction, powers, square roots, absolute values, and scalar-to-f64 display/oracle conversion through `SafeFromF64`, Eunomia `FloatElement`, and Eunomia `NumericElement`.
- **Boundary**: This slice covers `crates/cfd-1d/src/domain/channel`, the `domain/junctions/branching` solver/physics/validation cone, and `solver/analysis` aggregates/analyzers. It does not remove the direct `cfd-1d` manifest `num-traits` dependency because other domains still use `FromPrimitive`, `ToPrimitive`, and `Float`.
- **Evidence tier**: compile-time integration, empirical nextest, and static source audit. `cargo check -p cfd-1d` passed. `cargo nextest run -p cfd-1d` passed 725/725 tests with 3 skipped. Focused source scans found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`, generic `.to_f64()`, or `Float::` residue in the touched channel/branching/analyzer cone or in `solver/analysis`.
- **Residual risk**: Full `cfd-1d` all-target clippy remains blocked by unrelated existing lint debt in examples, tests, and cell-separation/ resistance modules. The separate domain-components direct `num-traits` residue is closed by the later domain-components slice; broader direct provider residue remains in solver-core, domain-network, vascular Bessel/Womersley, and other package areas outside this slice.

---

## Sprint 2026-07-04: cfd-1d Murray's-Law Eunomia Boundary Cleanup
- **Resolved**: `cfd-1d` Murray's-law vascular geometry no longer depends on `num_traits::FromPrimitive`. `MurraysLaw` and `OptimalBifurcation` route scalar constants, power functions, absolute value, and inverse cosine through Eunomia `FloatElement`/`NumericElement`. This uses the current Eunomia `FloatElement::acos` surface, which already has value-semantic tests in the provider checkout.
- **Boundary**: This slice covers only `crates/cfd-1d/src/physics/vascular/murrays_law/{law,bifurcation}.rs` plus the required bound propagation in `physics/vascular/bifurcation.rs`. `cfd-1d` still declares `num-traits` because many other modules still use `FromPrimitive`, `ToPrimitive`, and `Float`.
- **Evidence tier**: static source audit plus formatting. `cargo fmt -p cfd-1d --check` passed, and a focused scan over the Murray's-law files found no `num_traits`, `FromPrimitive`, `ToPrimitive`, or unqualified generic `acos`/`powf`/`powi`/`abs` residue. The provider inverse-cosine contract was re-verified in the Eunomia checkout with `cargo nextest run -p eunomia acos` (2/2).
- **Residual risk**: `cargo check -p cfd-1d` remains blocked by unrelated dirty-tree errors outside this slice, including missing `SafeFromF64` bounds, generic `abs`/`powf` ambiguity, and stale `to_f64().unwrap_or(...)` call sites. Full `cfd-1d` nextest evidence is pending behind that cleanup.

---

## Sprint 2026-07-03: cfd-core Fåhræus-Lindqvist Eunomia Scalars
- **Resolved**: `FahraeuasLindqvist<T>` no longer depends on direct `num_traits::FromPrimitive`, direct generic `T::from_f64()`, direct `T::zero()`, direct `T::one()`, direct generic `powf()`, direct generic `exp()`, or direct generic `abs()` for local scalar construction and microvascular viscosity formulas. Pries/Secomb exponent formulas, `mu_45`, relative-viscosity clamping, and tube hematocrit now route scalar constants and math through Eunomia `FloatElement`/`NumericElement`.
- **Boundary**: This slice is limited to `crates/cfd-core/src/physics/fluid/blood/fahraeus_lindqvist.rs`. Together with the prior Cross/Casson/Carreau slices, local blood-model scalar `num_traits` construction is closed. The broader `Fluid<T>` trait still carries an inherited `nalgebra::RealField` boundary for later provider work.
- **Evidence tier**: compile-time provider integration, existing value-semantic blood tests, and static source audit. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core fahraeus_lindqvist` passed 3/3. `cargo nextest run -p cfd-core blood` passed 24/24. Touched-file rustfmt and touched-file `git diff --check` passed. Focused `fahraeus_lindqvist.rs` residue scan found no direct `num_traits`, `FromPrimitive`, generic `T::from_f64`, `T::zero`, `T::one`, generic `powf`, generic `exp`, or generic `abs` residue. Broader blood residue scan now matches only concrete `f64` helper/test expressions.
- **Residual risk**: `cargo clippy -p cfd-core --all-targets -- -D warnings` remains blocked by pre-existing unrelated lints in `crates/cfd-core/src/physics/boundary/applicator.rs` and `crates/cfd-core/src/physics/fluid_dynamics/rhie_chow.rs`.

## Sprint 2026-07-03: cfd-core Casson/Carreau Blood Eunomia Scalars
- **Resolved**: `CassonBlood<T>`, `CarreauYasudaBlood<T>`, and `BloodModel<T>` no longer require direct `num_traits::FromPrimitive` for local scalar construction or model dispatch. Casson constants, hematocrit scaling, temperature correction, square-root apparent-viscosity formula, and validation constants now route through Eunomia `FloatElement`/ `NumericElement`. Carreau-Yasuda constants, zero/one identities, and real powers remain on Eunomia after the dispatch-bound fix.
- **Boundary**: This slice is limited to `crates/cfd-core/src/physics/fluid/blood/{casson,carreau_yasuda,mod}.rs`. The concrete `temperature_viscosity_factor(f64)` helper remains intentionally concrete. Fåhræus-Lindqvist and the broader `Fluid<T>` trait still carry residual blood-fluid provider migration work.
- **Evidence tier**: compile-time provider integration, existing value-semantic blood tests, and static source audit. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core casson` passed 12/12. `cargo nextest run -p cfd-core carreau_yasuda` passed 4/4. `cargo nextest run -p cfd-core blood` passed 24/24. Touched-file rustfmt and touched-file `git diff --check` passed. Focused touched-file residue scan found no direct `num_traits`, `FromPrimitive`, generic `T::from_f64`, `T::zero`, `T::one`, generic `powf`, generic `sqrt`, or generic `exp` residue; the only remaining match is the intentionally concrete `f64` temperature helper.
- **Residual risk**: `cargo clippy -p cfd-core --all-targets -- -D warnings` remains blocked by pre-existing unrelated lints in `crates/cfd-core/src/physics/boundary/applicator.rs` and `crates/cfd-core/src/physics/fluid_dynamics/rhie_chow.rs`.

## Sprint 2026-07-03: cfd-core Cavitation Eunomia Scalar Cone
- **Resolved**: The touched cavitation Rayleigh-Plesset, biological damage, regime-analysis, cavitation-number, and material-damage surfaces no longer depend on `nalgebra::RealField`, direct `num_traits::FromPrimitive`, direct `T::zero()`, direct `T::one()`, or direct `T::from_f64()` construction. Scalar constants, powers, square roots, exponentials, finite checks, and scalar min/max now route through Eunomia `FloatElement`/`NumericElement`.
- **Boundary**: This slice is limited to `crates/cfd-core/src/physics/cavitation/{rayleigh_plesset,bio_damage,number,damage}.rs` and `crates/cfd-core/src/physics/cavitation/regimes/`. Remaining cavitation scalar holdouts at that point were `models.rs`, `venturi.rs`, `nuclei_transport.rs`, and `heterogeneous_nucleation.rs`; the later Venturi slice removed `venturi.rs` from the active residual list.
- **Evidence tier**: compile-time provider integration, value-semantic focused tests, and static source audit. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core cavitation` passed 35/35 tests. Focused residue scan across the migrated cavitation files found no `nalgebra`, `RealField`, `num_traits`, `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`, direct `T::from_u64`, direct `T::zero`, direct `T::one`, direct `powf`, or direct `powi` residue. Touched-file rustfmt and touched-file `git diff --check` passed. Full `cargo clippy -p cfd-core --all-targets -- -D warnings` remains blocked by unrelated existing lints in `physics/boundary/applicator.rs` and `physics/fluid_dynamics/rhie_chow.rs`.

## Sprint 2026-07-03: cfd-core Material Eunomia Traits
- **Resolved**: `SolidProperties`, `InterfaceProperties`, `ElasticSolid`, `WettingProperties`, and `FluidSolidInterface` no longer depend on `nalgebra::RealField`; scalar math and constant construction now use Eunomia `FloatElement`/`NumericElement`.
- **Boundary**: `MaterialDatabase` still requires `RealField` because it stores `Box<dyn Fluid<T>>`; fluid, hemolysis, fluid-dynamics, boundary, geometry, mesh, and solver nalgebra surfaces remain open.
- **Evidence tier**: compile-time integration plus empirical focused tests and static source audit. Focused material scan shows no `nalgebra::RealField`, `RealField`, `num_traits`, `FromPrimitive`, `ToPrimitive`, or `Float` residue in the migrated solid/interface files. `cargo check -p cfd-core` passed, and `cargo nextest run -p cfd-core material` passed 4/4 tests. No runtime performance claim is made.

## Sprint 2026-07-03: cfd-core Velocity Leto Vector
- **Resolved**: `cfd-core::physics::values::Velocity` now stores `leto::geometry::Vector3<T>` instead of `nalgebra::Vector3<T>`, and its generic contract is Eunomia `FloatElement`/`NumericElement` instead of `nalgebra::RealField`. `PhysicalParameters::gravity` now uses the same Leto vector type and no longer requires `RealField` for its own methods.
- **Provider extension**: Leto geometry now derives Serde for `Point2`, `Point3`, `Vector3`, `UnitVector3`, and `Isometry3`, preserving CFDrs' serialized value-object boundary while using provider-owned vector storage.
- **Boundary**: `ProblemAggregate` and `SimulationAggregate` still carry `RealField` because their `Domain<T>` and fluid-property contracts have not been migrated in this slice. Material, hemolysis, and fluid-dynamics nalgebra-bound surfaces remain open.
- **Evidence tier**: compile-time integration plus static source audit. Touched-file rustfmt passed. Focused residue scan shows no nalgebra `Vector3`/`RealField` in `Velocity` or `PhysicalParameters`. `cargo check -p cfd-core` passed after compiling the modified local Leto provider. `cargo nextest run -p cfd-core --lib` passed 183/183 tests. Downstream `cargo check -p cfd-2d`, `cargo check -p cfd-3d`, and `cargo check -p cfd-validation` passed. No runtime performance claim is made.

## Sprint 2026-07-03: cfd-core Physics Value Eunomia Scalars
- **Resolved**: Scalar-only physics value wrappers no longer use `nalgebra::RealField` or `nalgebra::ComplexField`. `Temperature`, `Pressure`, `ReynoldsNumber`, and `DimensionlessNumber` now depend on Eunomia `FloatElement`/`NumericElement` for scalar construction, zero, absolute value, and square root. Immediate aggregate owners now declare the same bound where they store these value wrappers.
- **Boundary**: `Velocity` remains on `nalgebra::Vector3`, so its `RealField` bound is intentionally preserved for the Leto/Gaia vector replacement slice. Additional `RealField` use in material, hemolysis, and broader aggregate contracts remains open.
- **Evidence tier**: compile-time integration plus static source audit. Touched-file rustfmt passed. A focused scan of the four scalar value-wrapper files returns no `nalgebra::RealField`, `RealField`, or `ComplexField` matches. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --lib` passed 183/183 tests. No runtime performance claim is made.
- **Residual risk**: Full nalgebra removal from `cfd-core` still requires replacing vector/matrix contracts and material/hemolysis scalar bounds with Atlas-owned providers.

## Sprint 2026-07-03: cfd-1d Non-Python ndarray Path Removal
- **Resolved**: `cfd-1d` no longer declares the unused `sprs` dependency that pulled `ndarray v0.17.2` into the active 1D/3D dependency graph. The root workspace no longer declares unused `ndarray` as a shared dependency. cfd-1d tests that construct `ConstantPropertyFluid::water_20c()` now declare the required `eunomia::FloatElement` bound explicitly.
- **Dependency audit**: `cargo tree -p cfd-1d -i ndarray` and `cargo tree -p cfd-3d -i ndarray` report no matching package. `cargo tree --workspace -i ndarray` shows the remaining workspace path as `ndarray v0.16.1 -> numpy v0.22.1 -> cfd-python`.
- **Evidence tier**: manifest/lock static audit plus compile-time integration and empirical focused tests. Touched-file rustfmt passed. `cargo update -p sprs` removed `sprs`, `ndarray v0.17.2`, and stale transitive packages. `cargo check -p cfd-1d` passed. `cargo nextest run -p cfd-1d` passed 725/725 tests with 3 skipped.
- **Residual risk**: cfd-1d still uses `nalgebra`/`nalgebra-sparse`; replacing those surfaces with Leto sparse/vector types remains open. Workspace `ndarray` remains through the Python `numpy` boundary, which needs a separate binding API decision.

## Sprint 2026-07-03: cfd-math Eunomia Multigrid Interpolation
- **Resolved**: `cfd-math::linear_solver::preconditioners::multigrid::interpolation` has been changed to remove direct `num_traits::{FromPrimitive, ToPrimitive}` imports and direct scalar conversion/fallback paths. Interpolation scalar constants, index-distance conversion, quality row-sum extraction, constant preservation error, sparsity ratio, and absolute-value dispatch now route through Eunomia helpers.
- **Boundary**: This slice preserves the existing nalgebra `DVector` and sparse/vector surfaces. GMG transfer, AMG residual bounds, and full Leto migration remain separate provider slices.
- **Evidence tier**: compile-time integration plus empirical focused interpolation tests and static source audit. Touched-file rustfmt passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math interpolation` passed 14/14 tests. Static scan found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`, direct `T::from_usize`, conversion fallback, `from_f64_or`, `SafeFromF64`, stale `rayon`, direct `as f64`, or `.to_f64()` fallback hits in `multigrid/interpolation.rs`.
- **Residual risk**: `multigrid/gmg` and `amg.rs` still retain direct provider residue; nalgebra sparse/vector surfaces remain open for Leto migration.

## Sprint 2026-07-03: cfd-math Eunomia Multigrid Smoothers
- **Resolved**: `cfd-math::linear_solver::preconditioners::multigrid::smoothers` no longer uses direct `T::from_f64(...).unwrap_or_else(...)` conversion fallbacks or nalgebra absolute-value dispatch for smoother diagonal thresholds, Chebyshev eigenvalue defaults, Chebyshev recurrence constants, or smoother update thresholds. The touched AMG owner paths now construct coarsening thresholds, smoother relaxation values, and complexity filters through Eunomia.
- **Boundary**: This slice preserves the existing nalgebra `DVector` and `SparseMatrix` surfaces and keeps AMG's `FromPrimitive` bound because deeper coarsening/interpolation contracts still require it.
- **Evidence tier**: compile-time integration plus empirical focused value-semantic smoother tests and static source audit. Touched-file rustfmt passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math test_gauss_seidel_smoother test_jacobi_smoother test_symmetric_gauss_seidel test_sor_smoother test_chebyshev_smoother` passed 5/5 tests. Static scan found no direct scalar-conversion fallback, stale `SafeFromF64`, `from_f64_or`, direct `T::from_usize`, or stale `rayon` hits in the touched smoother/AMG files, and no direct `num_traits` provider residue in `smoothers.rs`.
- **Residual risk**: `amg.rs` still imports and bounds direct `num_traits::FromPrimitive` for deeper coarsening/interpolation routines. Other multigrid modules, the raw GPU operator path, and nalgebra sparse/vector surfaces remain open provider-migration work.

## Sprint 2026-07-03: cfd-math Eunomia Stability Scalars
- **Resolved**: `cfd-math::time_stepping::stability` no longer imports `num_traits::ToPrimitive` or uses direct `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallback conversions for analyzer defaults, CFL thresholds, RK order checks, or von Neumann amplification outputs.
- **Boundary**: This slice preserves the existing nalgebra `DMatrix`/`DVector` Butcher-tableau API. Leto replacement for that matrix/vector surface remains a separate migration item.
- **Evidence tier**: compile-time integration plus empirical value-semantic stability tests and static source audit. Touched-file `rustfmt --check` passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math stability` passed 5/5 tests. `git diff --check` passed. Static scan found no `num_traits`, `ToPrimitive`, `FromPrimitive`, `T::from_f64`, or conversion fallback hits under `crates/cfd-math/src/time_stepping/stability`.
- **Residual risk**: Package-level `cargo fmt --package cfd-math --check` is still blocked by unrelated existing formatting drift outside the touched stability files. Broader cfd-math still contains direct `num-traits`, nalgebra matrix/vector surfaces, and GPU provider gaps.

## Sprint 2026-07-03: cfd-core Eunomia Physics Value Boundary
- **Resolved**: `cfd-core` physics value objects (`Velocity`, `Temperature`, `Pressure`, `ReynoldsNumber`, `DimensionlessNumber`) and their dependent management aggregates no longer use direct `num_traits::FromPrimitive` bounds or `T::from_f64(...).unwrap_or_else(T::zero/one)` scalar constant conversions.
- **Boundary**: `nalgebra::RealField` remains in this slice because the touched value objects still store `nalgebra::Vector3` values and use nalgebra vector/scalar methods. The cfd-core manifest still retains `num-traits` because `management::conversion` and other modules outside the touched cone still use direct `num-traits`.
- **Evidence tier**: compile-time integration plus empirical value-semantic regression tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo nextest run -p cfd-core --no-default-features` passed 144/144 tests. Touched-file `rustfmt --check` passed. Static scan over the touched cone found no `num_traits`, `FromPrimitive`, `T::from_f64`, or silent fallback conversion hits.
- **Residual risk**: Full cfd-core Eunomia migration still needs `management::conversion`, remaining numeric helper modules, and eventual Leto replacement for nalgebra-owned vector/matrix state. Package-level `cargo fmt --check --package cfd-core` is still blocked by pre-existing unrelated formatting drift in `error.rs`, `physics/boundary/ghost_cells.rs`, and `physics/fluid_dynamics/operations.rs`.

---

## Phase 1: Foundation (Current)
- Auditing `cfd-schematics`...
- Goal: Identify redundancy, shared component duplication, and structural integrity violations.
- Checking for placeholders, temporary workarounds, and approximations.

## Findings
### 1. Dual-Path Topology Generation (SSOT & DRY Violation)
- `interface::presets` (e.g., `bifurcation.rs`, `serpentine.rs`) directly instantiate `NetworkBlueprint`, manually inserting nodes and channels.
- `topology::presets` combined with `topology::factory::BlueprintTopologyFactory` does exactly the same thing but via declarative `BlueprintTopologySpec`.
- **Gap**: `interface::presets` must be completely removed or refactored into thin wrappers that call `topology::presets` and `BlueprintTopologyFactory`.

### 2. Leaking Physics Constants (SOC Violation)
- `BLOOD_MU` (3.5e-3) and `shah_london_resistance` / `hp_resistance` formulas are hardcoded into schematic generation (`cfd-schematics/src/interface/presets/*` and `factory.rs`).
- **Gap**: Schematic generation (topology/geometry mappings) should not calculate fluidic resistances. That is the domain of `cfd-1d`. The schematic layer should only provide geometry (`length_m`, `width_m`, `height_m`) and let downstream solvers compute resistances based on their own runtime fluid models.

---

# Finding 2026-07-10: cfd-1d double-trifurcation Picard non-convergence (test regression)

- **Symptom**: `cfd-suite` integration test `tests/cross_fidelity_blueprint.rs::cross_fidelity_blueprint_complex_branching` fails. `solve_reference_trace` → cfd-1d `Network2DSolver` returns `MaxIterationsExceeded: Maximum iterations (10000) exceeded` on the `double_trifurcation_cif_venturi_rect` network (blood ρ=1060, µ=3.5e-3, Q=1e-7 m³/s). The sibling `cross_fidelity_blueprint_bifurcation` passes.
- **Isolation**: Reproduces with `-p cfd-suite --no-default-features` (gpu OFF), so it is independent of the concurrent cfd-core `gpu/turbulence_compute` directory-module refactor that currently breaks the `gpu` feature build (`E0583 turbulence_compute` — peer WIP, not this finding).
- **Locus**: `crates/cfd-1d/src/solver/core/mod.rs` Picard fixed-point loop (line ~358) with Anderson acceleration + Newton fallback; `has_converged_dual` never satisfied within 10000 iters for this stiff nonlinear network. Solver correctly returns a typed error (no silent bad result).
- **Hypothesis (unverified)**: migration regression from the Leto-CSR assembly push (`d58d1fe3`/`1d768895`) subtly altering the assembled matrix, or a genuine conditioning/stiffness limit the current Picard+Anderson path cannot crack (cf. recent peer work `0d101352` "enhance Anderson QR collapse detection" — actively-evolving convergence path).
- **DoR to resolve**: (1) differential-test the assembled cfd-1d matrix + rhs for this network against the pre-migration commit (parent of `d58d1fe3`) to classify regression vs. genuine stiffness; (2) capture the Picard residual trajectory (diverge / plateau / slow-decay) via a bounded-iteration diagnostic; (3) if genuine stiffness, verify the Newton fallback engages on Picard stagnation. Blocked from immediate fix: convergence path is under active concurrent peer edit (`0d101352`); coordinate before touching `solver/core`.
- **Evidence tier**: empirical (reproduced deterministically, gpu-independent). Not test-gaming: the test asserts real mass-conservation physics; the fix must be in the solver/assembly, never a weakened tolerance or raised iteration cap.
