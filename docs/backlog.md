# CFD Suite Backlog

## Sprint 1.96.13: cfd-2d Dependency-Aware Physics Audit
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Audit `cfd-2d` and its direct internal dependencies (`cfd-core`, `cfd-math`, `cfd-mesh`, `cfd-io`, `cfd-1d`, `cfd-schematics`) for role alignment in 2D Navier-Stokes, LBM, and schematics-projected channel solves.
- Correct the highest-risk local physics defect found during the audit without changing the public solver API.

### Dependency Audit Findings
- `cfd-core`: appropriate for shared fluid physics, boundary condition types, errors, constants, and GPU feature plumbing; `cfd-2d` should not duplicate those foundation contracts.
- `cfd-math`: appropriate for sparse/iterative numerical kernels and high-order schemes; local discretization remains in `cfd-2d` where it depends on 2D stencil geometry.
- `cfd-mesh`: appropriate for mesh-backed geometry inputs and unstructured contexts; structured solver kernels correctly keep 2D grid storage local.
- `cfd-io`: appropriate for field output and visualization export; it is not part of the time-integration physics path.
- `cfd-1d`: appropriate as a reduced-order pressure/flow reference for channel projection and cross-fidelity seeding, not as a replacement for 2D PDE solves.
- `cfd-schematics`: appropriate as topology and layout authority; `cfd-2d` consumes blueprints through projection metadata instead of owning schematic generation.

### Sprint Backlog Items

#### 2D LBM Low-Mach Physics
- [x] **LBM-MACH-001 [patch]**: Enforce `Ma <= 0.1` during LBM initialization.
- [x] **LBM-MACH-002 [patch]**: Enforce `Ma <= 0.1` for Zou-He velocity and pressure boundary reconstruction.
- [x] **LBM-MACH-003 [patch]**: Propagate boundary-condition Mach violations through `LbmSolver::step`.
- [x] **LBM-MACH-004 [patch]**: Add value-semantic tests for high-Mach initialization and inlet rejection.
- [x] **LBM-MACH-005 [patch]**: Verify the touched `cfd-2d` LBM module with bounded Cargo check, unit test, clippy, and nextest runs.

## Sprint 1.96.12: cfd-1d Dependency-Aware Physics Audit
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Audit `cfd-1d` and its direct internal dependencies (`cfd-core`, `cfd-math`, `cfd-schematics`) for role alignment in the 1D millifluidic solver.
- Correct the highest-risk local physics defect found during the audit without changing public APIs.

### Dependency Audit Findings
- `cfd-core`: appropriate as the source for fluid properties, constants, boundary conditions, shared errors, and solver traits; `cfd-1d` uses it at physics and solver boundaries rather than duplicating fluid databases.
- `cfd-math`: appropriate for sparse and iterative linear algebra in network solves; `cfd-1d` also uses `nalgebra-sparse` directly for assembly, which is acceptable because it owns the graph-to-matrix mapping.
- `cfd-schematics`: appropriate as topology and cross-section authority; `cfd-1d` consumes `NodeSpec`, `ChannelSpec`, and render metadata but keeps hydraulic resistance and flow-state logic local.
- External crates: `petgraph` is appropriate for directed network topology, `nalgebra` for scalar/vector algebra, `serde` for persisted design state, `rayon` only where independent network paths can be parallelized.

### Sprint Backlog Items

#### 1D Serpentine Reverse-Flow Physics
- [x] **SERP-REV-001 [patch]**: Use flow-speed magnitude for serpentine shear rate, Reynolds number, Dean number, friction losses, and bend minor losses.
- [x] **SERP-REV-002 [patch]**: Preserve orientation outside scalar resistance coefficients by keeping resistance magnitudes invariant under `Q -> -Q`.
- [x] **SERP-REV-003 [patch]**: Add value-semantic tests for reverse-flow coefficient and analysis invariance.
- [x] **SERP-REV-004 [patch]**: Verify the touched `cfd-1d` serpentine module with bounded Cargo check, unit test, clippy, and nextest runs.

## Sprint 1.96.11: cfd-3d Sigma SGS Energy Physics
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Correct `cfd-3d` Sigma LES SGS kinetic-energy diagnostics so they use the same Yoshizawa eddy-viscosity inversion as the other algebraic LES closures.
- Consolidate SGS energy computation through the shared turbulence helper instead of retaining a Sigma-specific strain-rate expression with incorrect dimensions.

### Sprint Backlog Items

#### 3D LES Closure Physics
- [x] **SIGMA-SGS-001 [patch]**: Replace `nu_t |S| / Delta` Sigma SGS energy with `(nu_t / (C_k Delta))^2`.
- [x] **SIGMA-SGS-002 [patch]**: Route Sigma SGS energy through the shared `kinetic_energy_from_eddy_viscosity` SSOT.
- [x] **SIGMA-SGS-003 [patch]**: Add value-semantic regression coverage for the Yoshizawa relation and the former dimensional mismatch.
- [x] **SIGMA-SGS-004 [patch]**: Verify the touched `cfd-3d` Sigma module with bounded Cargo check, unit test, clippy, and nextest runs.

## Sprint 1.96.10: cfd-3d VOF Directional CFL Physics
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Correct `cfd-3d` VOF timestep selection so it targets the same summed directional CFL invariant enforced by geometric advection.
- Add value-semantic regression coverage for diagonal flow that distinguishes the corrected reciprocal summed-rate bound from the former norm/min-spacing estimate.

### Sprint Backlog Items

#### 3D VOF Stability Physics
- [x] **VOF-CFL-001 [patch]**: Replace speed/min-spacing timestep selection with `C / max(|u_x|/dx + |u_y|/dy + |u_z|/dz)`.
- [x] **VOF-CFL-002 [patch]**: Document the explicit VOF CFL theorem at the solver timestep API.
- [x] **VOF-CFL-003 [patch]**: Add regression coverage that accepts the corrected diagonal-flow timestep and rejects the former estimate.
- [x] **VOF-CFL-004 [patch]**: Verify the touched `cfd-3d` VOF module with bounded Cargo check, integration test, and nextest runs.

## Sprint 1.96.9: cfd-2d Upwind Coefficient Orientation
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Correct `cfd-2d` finite-volume first-order upwind coefficients so west-face convection uses the Patankar face-orientation sign convention.
- Add value-semantic tests that distinguish positive and negative west-face flux orientation.

### Sprint Backlog Items

#### 2D Convection Discretization Physics
- [x] **CONV-001 [patch]**: Replace `a_W = D_W + max(-F_W, 0)` with `a_W = D_W + max(F_W, 0)` in `FirstOrderUpwind`.
- [x] **CONV-002 [patch]**: Document the finite-volume upwind coefficient invariant for east and west faces.
- [x] **CONV-003 [patch]**: Add regression coverage for west-face positive and negative flux orientation.
- [x] **CONV-004 [patch]**: Verify the touched `cfd-2d` convection module with bounded Cargo check, unit test, and nextest runs.

## Sprint 1.96.8: cfd-1d Branch Reverse-Flow Physics
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Correct `cfd-1d` branch-junction solvers so reversed parent flow preserves daughter orientation instead of being rejected.
- Keep wall-shear and non-Newtonian viscosity inputs magnitude-based while pressure drops remain signed.

### Sprint Backlog Items

#### 1D Branch-Junction Orientation Physics
- [x] **BRC-001 [patch]**: Remove nonnegative parent-flow rejection from two-way and three-way branch solvers.
- [x] **BRC-002 [patch]**: Solve branch splits on `|Q_parent|` and reapply flow orientation to daughter flows.
- [x] **BRC-003 [patch]**: Compute wall-shear and apparent viscosity from flow magnitude while preserving signed pressure drops.
- [x] **BRC-004 [patch]**: Add value-semantic reverse-flow tests for pressure-balanced and prescribed two-way/three-way branches.
- [x] **BRC-005 [patch]**: Verify the touched `cfd-1d` branch module with bounded Cargo check, unit test, and nextest runs.

## Sprint 1.96.7: cfd-3d Venturi Pressure-Coefficient Physics
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Correct `cfd-3d` Venturi pressure coefficients so the solver follows its documented throat dynamic-pressure definition.
- Add value-semantic tests for coefficient scaling and undefined zero-flux rejection.

### Sprint Backlog Items

#### 3D Venturi Coefficient Physics
- [x] **V3D-001 [patch]**: Replace inlet dynamic-pressure scaling with throat dynamic pressure derived from face-integrated flow and throat area.
- [x] **V3D-002 [patch]**: Preserve typed rejection when coefficient scaling has no positive throat flow or area.
- [x] **V3D-003 [patch]**: Add direct regression tests for throat coefficient values and zero-flux rejection.
- [x] **V3D-004 [patch]**: Verify the touched `cfd-3d` Venturi module with bounded Cargo check, unit test, and nextest runs.

## Sprint 1.96.6: cfd-2d Explicit Stability Physics
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Correct `cfd-2d` explicit time-step physics so 2D advection and diffusion bounds follow summed-rate CFL and von Neumann stability conditions.
- Add value-semantic regression tests that verify `max_stable_dt` satisfies the same stability invariants reported by `advection_cfl` and `diffusion_number`.

### Sprint Backlog Items

#### 2D Stability Physics
- [x] **CFL-001 [patch]**: Replace componentwise minimum advection time-step bound with `1 / (|u|/dx + |v|/dy)`.
- [x] **CFL-002 [patch]**: Replace `0.5*min(dx²,dy²)/ν` diffusion bound with `0.5 / (ν(1/dx² + 1/dy²))`.
- [x] **CFL-003 [patch]**: Add regression tests for summed 2D advection and diffusion stability limits.
- [x] **CFL-004 [patch]**: Verify the touched `cfd-2d` stability module with bounded Cargo check, unit test, and nextest runs.

## Sprint 1.96.5: cfd-1d Venturi Reverse-Flow Physics
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Correct `cfd-1d` Venturi resistance physics so reverse-flow inputs use the same scalar loss magnitudes as forward flow.
- Preserve signed reported velocities while making Reynolds number, viscosity, friction factor, and pressure-loss coefficients depend on velocity magnitude.

### Sprint Backlog Items

#### 1D Venturi Resistance Physics
- [x] **VNT-001 [patch]**: Make `VenturiModel::calculate_coefficients` decompose losses with `|V|` and `|Q|` instead of treating negative inlet velocity as zero flow.
- [x] **VNT-002 [patch]**: Make detailed Venturi analysis use shear-rate and Reynolds magnitudes for reverse-flow cases.
- [x] **VNT-003 [patch]**: Add regression tests for symmetric-geometry coefficient and pressure-loss orientation invariance.
- [x] **VNT-004 [patch]**: Verify the touched `cfd-1d` venturi module with bounded Cargo check, unit test, and nextest runs.

## Sprint 1.96.4: Geometric Conservation Residual Verification
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Replace copy-through GCL evolution in `cfd-validation` with an executable conservative finite-volume residual.
- Verify Euler and Runge-Kutta constant-state preservation while proving non-constant fields are evaluated by the residual.

### Sprint Backlog Items

#### Conservation Validation Physics
- [x] **GCL-001 [patch]**: Implement a conservative second-order face-flux residual for the geometric conservation checker.
- [x] **GCL-002 [patch]**: Run Euler, midpoint, SSPRK3, and RK4 checks through real residual-based stage updates.
- [x] **GCL-003 [patch]**: Add value-semantic tests for unsupported RK stages and quadratic-field residual sensitivity.
- [x] **GCL-004 [patch]**: Verify the touched library target with bounded Cargo check, unit test, and nextest runs.

## Sprint 1.96.3: Womersley Analytical SSOT
**Status**: Completed
**Start Date**: May 3, 2026

### Sprint Objectives
- Replace the validation-local approximate Womersley profile with the canonical exact Bessel implementation from `cfd-1d`.
- Add value-semantic no-slip verification for the wall velocity invariant.

### Sprint Backlog Items

#### Analytical Validation SSOT
- [x] **WOM-001 [patch]**: Delegate `cfd-validation` Womersley velocity, wall shear stress, and flow rate to the `cfd-1d` exact `WomersleyProfile`.
- [x] **WOM-002 [patch]**: Add wall no-slip regression coverage for multiple phases.
- [x] **WOM-003 [patch]**: Verify the touched analytical module with bounded Cargo checks and tests.

## Sprint 1.96.2: Optimization Terminology Contract Cleanup
**Status**: Completed
**Start Date**: May 3, 2026

### Sprint Objectives
- Remove misleading unsupported/stub-like wording from executable boundary and optimization contracts without changing numerical behavior.
- Keep the `cfd-schematics` serpentine optimization path generator named by role rather than by reduced implementation detail.

### Sprint Backlog Items

#### Contract Terminology Cleanup
- [x] **TERM-001 [patch]**: Rename internal `generate_simplified_serpentine_path` to `generate_optimization_serpentine_path` and update all optimization call sites.
- [x] **TERM-002 [patch]**: Replace boundary stencil `"not implemented"` error text with explicit unsupported-order wording.
- [x] **TERM-003 [patch]**: Verify touched crates with bounded Cargo check and record nextest compile-bound status.

## Sprint 1.96.1: Workspace SSOT Cleanup
**Status**: Completed
**Start Date**: April 29, 2026

### Sprint Objectives
- Remove obsolete tracked root-source artifacts that duplicate canonical `gaia` CSG and mesh assembly code.
- Preserve the Cargo workspace authority: only files referenced by package manifests, module trees, examples, tests, docs, or report-generation assets remain as source artifacts.

### Sprint Backlog Items

#### Source Authority Cleanup
- [x] **SSOT-001 [patch]**: Delete unreferenced root historical Rust artifacts `old_assemble.rs`, `old_arrangement.rs`, `old_phase2.rs`, `old_operations.rs`, `old_indexed.rs`, `old_gwn_bvh.rs`, `old_seam.rs`, `old_phase4.rs`, and empty `csg_bi.rs`.
- [x] **SSOT-002 [patch]**: Verify no Cargo manifest, crate module, example, test, documentation, or report artifact references the deleted files.
- [x] **SSOT-003 [patch]**: Configure Windows GNU builds to use MSYS2 clang/lld and LLVM archive tools for Rust linking and C/C++ build scripts.
- [x] **SSOT-004 [patch]**: Remove misleading audit-trigger terminology from explicit unsupported-operation paths and already-bounded model comments.

## Sprint 1.96.0: Hydrodynamic Cavitation On-a-Chip (HCOC) Modeling
**Status**: Completed
**Start Date**: Next

### Sprint Objectives
- Enable mathematical modeling of distinct cellular injury zones (lysis, necrosis, permeabilization) from hydrodynamic cavitation in micro-orifices.
- Implement stiffness-coupled heterogeneous nucleation to predict altered cavitation inception times for Circulating Tumor Cells (CTCs).
- Provide unified mathematical specifications and zero-cost `GhostCell` topological integration for the new physics models.

### Sprint Backlog Items

#### Cellular Injury Modeling
- [x] **HCOC-001**: Implement `bio_damage.rs` with scalar fractional damage models representing cell membrane strain, mapping cavitation intensity to lysis/necrosis ratios.
- [x] **HCOC-002**: Derive and integrate property tests validating conservation of cellular mass across injury states.

#### CTC Nucleation Transport
- [x] **HCOC-003**: Implement `heterogeneous_nucleation.rs` replacing the linear $k_n$ scalar with a tensorial or struct-based multi-population interfacial tension approach.
- [x] **HCOC-004**: Validate modified effective vapor pressure inception offsets against analytical expected bounds.
- [x] **HCOC-005**: Bind Apollo-backed periodic DNS transforms to a validated reusable `FftPlan3D`.
- [x] **HCOC-006**: Correct `cfd-2d` MUSCL3/QUICK face reconstruction and document the polynomial reproduction theorem.
- [x] **HCOC-007**: Separate `cfd-1d` margination wall-induced and shear-gradient inertial lift scaling and document the reference-equilibrium proof.
- [x] **HCOC-008**: Replace `cfd-2d` turbulence benchmark placeholder dispatch with a closed supported-model registry.
- [x] **HCOC-009**: Restore coupled cfd-1d blood pressure-event solves with finite startup viscosity and row-equilibrated hydraulic linear systems.
- [x] **HCOC-010**: Replace compact `cfd-1d` plasma-skimming screening with the threshold-aware Pries phase-separation model and Murray-inferred sibling geometry.
- [x] **HCOC-011**: Replace `cfd-2d` serpentine mixing's exponential estimate with the Neumann transverse-diffusion eigenfunction solution.
- [x] **HCOC-012**: Replace `cfd-3d` LES turbulent-kinetic-energy viscosity aliases with a shared Yoshizawa SGS energy relation.
- [x] **HCOC-013**: Remove silent clamps from `cfd-2d` Pries plasma-skimming and add checked physical-envelope regression tests.
- [x] **HCOC-014**: Replace `cfd-2d` WALE boundary zero-gradient reductions with second-order one-sided gradient stencils.
- [x] **HCOC-015**: Replace `cfd-3d` Spalart-Allmaras all-zero TKE output with a wall-distance Yoshizawa diagnostic.
- [x] **HCOC-016**: Reject uninitialized `cfd-3d` k-epsilon state instead of synthesizing zero turbulence fields.
- [x] **HCOC-017**: Replace `cfd-2d` Smagorinsky LES zero SGS energy/dissipation placeholders with Yoshizawa diagnostics and second-order boundary strain.
- [x] **HCOC-018**: Remove the `cfd-1d` margination singular wall cutoff and silent lateral-position clamp while preserving inertial-lift scaling tests.
- [x] **HCOC-019**: Remove stored `cfd-1d` droplet occupied-channel state so finite-length occupancy spans are the only authoritative representation.
- [x] **HCOC-020**: Remove residual `cfd-2d` Smagorinsky validation SGS floors and keep LES validation on physical zero-floor defaults.
- [x] **HCOC-021**: Correct Milestone 12 cross-track therapy utility so ultrasound-only SDT receives acoustic-cavitation delivery credit from channel resonance, cancer-center enrichment, and active therapy-zone fraction.

## Sprint 1.95.1: CFD-MESH 3D Performance & Memory Optimization
**Status**: In Progress
**Start Date**: March 29, 2026

### Sprint Objectives
- Eliminate WASM memory exhaustion (OOM) during high-resolution 3D mesh generation.
- Refactor `BowyerWatson3D` kernel for strictly zero-allocation execution inside the point insertion loop.
- Flatten memory structures (O(N) contiguous buffers) for WASM/heap efficiency.
- Prove invariants mathematically via property tests, ensuring exact numerical determinism.

### Sprint Backlog Items

#### Core Memory Flattening
- [x] **MESH-PERF-001**: Lift `HashMap` and `Vec` allocations out of `BowyerWatson3D::insert_point`.
- [x] **MESH-PERF-002**: Pre-allocate `IndexedMesh` internal arrays in `build_volume` using AABB heuristic.
- [x] **MESH-PERF-003**: Optimize `SdfMesher` BCC lattice generation loops to prevent over-allocation.

## Sprint 1.91.0: Advanced Validation Framework Expansion
**Status**: Archived
**Start Date**: November 2, 2025

### Sprint Objectives
- Extend Method of Manufactured Solutions (MMS) test coverage for complex geometries
- Implement Richardson extrapolation automation for error analysis
- Add performance benchmarking framework for production scaling validation
- Enhance validation suite with comprehensive error metrics and convergence studies

### Sprint Backlog Items

#### MMS Test Coverage Expansion
- [x] **MMS-001**: Implement MMS for non-rectangular domains (circular, annular geometries) - Tasks 1.1, 1.2 & 1.3 COMPLETED
- [ ] **MMS-002**: Add MMS validation for complex boundary conditions (mixed BC types)
- [ ] **MMS-003**: Extend MMS to turbulent flow regimes with manufactured turbulence quantities
- [ ] **MMS-004**: Implement MMS for multi-physics coupling (heat transfer, species transport)

#### Richardson Extrapolation Automation
- [x] **RE-001**: Create Richardson extrapolation core library for grid convergence studies ✅ COMPLETED
- [x] **RE-002**: Implement automated grid refinement strategies (geometric progression) ✅ COMPLETED
- [x] **RE-003**: Add error estimation and convergence rate calculation ✅ COMPLETED
- [x] **RE-004**: Integrate Richardson extrapolation with existing MMS framework ✅ COMPLETED

#### Performance Benchmarking Framework
- [ ] **PB-001**: Design comprehensive benchmarking suite for CFD operations
- [ ] **PB-002**: Implement memory usage profiling and optimization tracking
- [ ] **PB-003**: Add parallel scaling analysis for multi-core/multi-GPU configurations
- [ ] **PB-004**: Create performance regression detection and alerting

#### Validation Infrastructure
- [ ] **VI-001**: Enhance validation test organization with clear taxonomy
- [ ] **VI-002**: Implement automated validation report generation (HTML/PDF)
- [ ] **VI-003**: Add validation metrics dashboard and visualization
- [ ] **VI-004**: Create validation configuration management system

### Completed Items
- [x] **PRE-001**: Resolve all build, test, and example errors (except known SIMPLEC convergence issues)
- [x] **PRE-002**: Clean up compiler warnings and code quality issues
- [x] **PRE-003**: Fix boundary condition handling in momentum solver

### Known Issues
- **SIMPLEC Convergence**: SIMPLEC/PIMPLE solver has matrix construction issues resolved, but pressure Poisson system remains singular due to incorrect reference pressure handling. Requires deeper investigation of boundary condition implementation and pressure-velocity coupling algorithm.
- **GPU Feature Gaps**: Some GPU-related code is unused due to feature flag configuration

### Sprint Metrics
- **Test Coverage**: Maintain >80% coverage across all crates
- **Performance**: No regression in benchmark performance
- **Code Quality**: Zero warnings, comprehensive documentation
- **Validation**: All new MMS tests pass with expected convergence rates
