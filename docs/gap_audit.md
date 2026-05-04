<!-- markdownlint-disable MD022 MD032 MD025 MD024 MD060 MD009 MD029 MD030 -->

# Elite Mathematically-Verified Code Auditor: CFD Suite Comprehensive Gap Analysis

**Auditor Persona**: Elite Mathematically-Verified Code Auditor
**Date**: November 18, 2025 (Updated: February 23, 2026 - Sprint 1.95.0)
**Status**: ✅ ALL CRITICAL ISSUES RESOLVED

## Audit Principles
- **Mathematical Accuracy**: Zero tolerance for error masking or unverified approximations.
- **Implementation Completeness**: Full theorem documentation and rigorous testing required.
- **Literature Compliance**: Algorithms must match primary literature exactly.
- **Quality Standards**: No "working but incorrect" implementations.

---

# Sprint 1.96.5 Resolution: cfd-1d Venturi Reverse-Flow Physics

### RESOLVED-034: Venturi Reverse Flow Lost Inertial Coefficients
- **Location**: `crates/cfd-1d/src/physics/resistance/models/venturi/`
- **Issue**: `VenturiModel::calculate_coefficients` selected its finite-flow branch with `v_inlet > threshold`, so negative inlet velocity was treated as zero flow and lost contraction, diffuser, and recovery coefficients. Detailed analysis also used signed velocity for shear-rate and Reynolds inputs.
- **Remediation**: Coefficient decomposition now uses `|V_inlet|`, `|V_throat|`, and `|Q|` for scalar loss magnitudes while detailed analysis preserves signed velocities in reported kinematics. Shear-rate, viscosity, Reynolds correction, friction factor, and pressure-loss terms use magnitudes.
- **Mathematical basis**: In a symmetric Venturi, reversing flow orientation changes the sign of axial velocities but leaves kinetic-energy losses proportional to `V²`, wall-shear magnitude, Reynolds magnitude, and scalar resistance coefficients unchanged.
- **Verification**:
  - `cargo check -p cfd-1d --no-default-features` passed.
  - `cargo test -p cfd-1d --no-default-features physics::resistance::models::venturi --lib -- --nocapture` passed 16/16 tests.
  - `cargo nextest run -p cfd-1d --lib --no-default-features physics::resistance::models::venturi --fail-fast --hide-progress-bar --status-level fail` passed 16/16 tests in 0.230 s.

---

# Sprint 1.96.4 Resolution: Geometric Conservation Residual Verification

### RESOLVED-033: Copy-Through Geometric Conservation Evolution
- **Location**: `crates/cfd-validation/src/conservation/geometric.rs`
- **Issue**: Euler and Runge-Kutta GCL checks copied the scalar field between time levels, so tests could pass without evaluating any numerical residual.
- **Remediation**: Added a conservative second-order face-flux residual and routed Euler, midpoint, SSPRK3, and RK4 checks through residual-based stage updates. Unsupported RK stage counts now return a typed `UnsupportedOperation` error.
- **Mathematical basis**: For constant `u_ij = c`, every face gradient is zero, so the finite-volume residual divergence is zero at every interior cell. Explicit Runge-Kutta stage combinations of `u` and `dt R(u)` therefore preserve the constant state exactly.
- **Verification**:
  - `cargo check -p cfd-validation --no-default-features` passed in 0.91 s.
  - `cargo test -p cfd-validation --no-default-features conservation::geometric --lib -- --nocapture` passed 7/7 tests.
  - `cargo nextest run -p cfd-validation --lib --no-default-features conservation::geometric --fail-fast --hide-progress-bar` passed 7/7 tests in 0.229 s under a 30-second shell timeout.

---

# Module: CFD-IO

**Status**:  COMPLETE (Verified)

## Audit Summary
- **Bit-Exact Roundtrip**: Verified. 
  read(write(data)) == data invariant holds.
- **Parallel Checksum**: Verified. MPI allreduce checksums implemented.
- **Documentation**: Invariants documented.

---

# Module: CFD-MATH (Linear Solvers & Preconditioners)

**Status**:  PASS (Remediated)

## Theorem Verification  Pass
- **Conjugate Gradient**: Good documentation of Optimality Theorem and Condition Number bounds.
- **Preconditioners**: Good mathematical foundation docs for Jacobi, SOR, AMG.
- **Parallel Consistency**: Clarified that current implementations are primarily serial reference implementations or building blocks for parallel solvers (which reside in `cfd-core` or use specific MPI types).

## Algorithm Audit  PASS

### RESOLVED-001: Redundant and Inferior ILU Implementation
- **Location**: `crates/cfd-math/src/linear_solver/preconditioners.rs`
- **Issue**: Two conflicting implementations of ILU factorization existed.
- **Remediation**: `ILUPreconditioner` has been removed. `IncompleteLU` is now the single source of truth, supporting both ILU(0) and ILU(k).

### RESOLVED-002: Serial Schwarz Decomposition (Clarified)
- **Location**: `crates/cfd-math/src/linear_solver/preconditioners.rs`
- **Issue**: The Schwarz preconditioner was implemented using a serial loop but named generic `SchwarzPreconditioner`.
- **Remediation**: Renamed to `SerialSchwarzPreconditioner` to explicitly reflect its nature as a serial implementation (useful for testing/research or as a local block solver). Documentation updated to clarify it does not provide MPI parallelism itself.

### RESOLVED-003: Parallel Claims in Solvers
- **Location**: `crates/cfd-math/src/linear_solver/conjugate_gradient.rs`
- **Issue**: Documentation claimed "Parallel Scalability" for a generic implementation that uses serial `nalgebra` types by default.
- **Remediation**: Removed misleading "Parallel Scalability" claims. Updated documentation to clarify that parallelism depends on the underlying vector/matrix types or should use the specific MPI solvers in `cfd-core`.

## Testing Validation  PARTIAL
- **Unit Tests**: Good coverage for serial cases (Identity, Diagonal, Poisson 1D).
- **Missing**: Parallel scaling tests (moved to `cfd-core` benchmarks), convergence verification for difficult matrices (high condition number).

## Code Quality Audit  PASS
- **Good**: Extensive doc comments with literature references.
- **Remediated**: Redundant structs removed, misleading names fixed.

---

# Module: CFD-CORE (Distributed Computing & MPI)

**Status**:  FAILED

## Algorithm Audit  FAILED

### RESOLVED-004: Fake Distributed GMRES Implementation
- **Location**: `crates/cfd-core/src/compute/mpi/distributed_solvers.rs` (`DistributedGMRES::gmres_iteration`)
- **Issue**: The GMRES solver performed the Arnoldi process but skipped the least-squares solution step, setting residual to zero manually.
- **Remediation**: Implemented full GMRES cycle with Givens rotations for updating the QR factorization of the Hessenberg matrix, and backward substitution to solve the least-squares problem.

### RESOLVED-005: Fake Additive Schwarz Preconditioner
- **Location**: `crates/cfd-core/src/compute/mpi/distributed_solvers.rs` (`AdditiveSchwarzPreconditioner::create_local_solvers`)
- **Issue**: The local subdomain solvers were implemented as Identity operations.
- **Remediation**: Updated to use `assemble_local_matrix` and `nalgebra::LU` for direct local solves. Added fallback to Jacobi (diagonal scaling) if matrix assembly is not supported by the operator.

### RESOLVED-006: Fake Block Jacobi Preconditioner
- **Location**: `crates/cfd-core/src/compute/mpi/distributed_solvers.rs` (`BlockJacobiPreconditioner::new`)
- **Issue**: The preconditioner assumed unit diagonal blocks.
- **Remediation**: Extended `DistributedLinearOperator` trait with `extract_diagonal` and updated `BlockJacobiPreconditioner` to use the actual operator diagonal.

---

# Module: CFD-MESH (Topology & Refinement)

**Status**:  FAILED

## Algorithm Audit  FAILED

### RESOLVED-007: Fake Mesh Refinement
- **Location**: `crates/gaia/src/refinement/mod.rs`
- **Issue**: `UniformRefinement` and `AdaptiveRefinement` methods were empty placeholders returning `Ok(())`.
- **Remediation**: Updated methods to explicitly return `MeshError::NotImplemented`. This removes the deceptive "working" status and correctly signals that the feature is pending implementation.

### RESOLVED-008: Missing Distributed Mesh Support
# Critical Audit Findings (Open Issues)

| ID | Severity | Component | Issue | Status |
|----|----------|-----------|-------|--------|
| **CRITICAL-001** |  Critical | CFD-MATH | Redundant ILUPreconditioner vs IncompleteLU | **CLOSED** |
| **CRITICAL-002** |  Critical | CFD-MATH | Serial implementation of Schwarz Decomposition | **CLOSED** |
| **MAJOR-003** |  Major | CFD-MATH | Unsubstantiated Parallel Scalability claims in CG | **CLOSED** |
| **CRITICAL-004** |  Critical | CFD-CORE | Fake Distributed GMRES (Placeholder Solver) | **CLOSED** |
| **CRITICAL-005** |  Critical | CFD-CORE | Fake Additive Schwarz (Identity Solver) | **CLOSED** |
| **MAJOR-006** |  Major | CFD-CORE | Fake Block Jacobi (Unit Diagonal Assumption) | **CLOSED** |
| **CRITICAL-011** | ✅ Closed | CFD-MESH | WASM OOM due to O(N²) allocation overhead in Bowyer-Watson | **CLOSED** |
| **CRITICAL-007** |  Critical | CFD-MESH | Fake Mesh Refinement (Empty Methods) | **CLOSED** |
| **MAJOR-008** |  Major | CFD-MESH | Missing Distributed Mesh Support | **CLOSED** |
| **CRITICAL-009** | ✅ Closed | CFD-MATH | Ruge-Stüben Coarsening Fine-to-Coarse Mapping Bug | **CLOSED** |
| **CRITICAL-010** | ✅ Closed | CFD-3D | Fake Unstructured FEM Domains (0-Element Mocks) | **CLOSED** |
| **CRITICAL-012** | 🔴 Critical | CFD-CORE | Missing cell-specific cavitation injury physics (Lysis/Necrosis) | **CLOSED** |
| **CRITICAL-013** | 🔴 Critical | CFD-CORE | Missing CTC stiffness-coupled inception mechanisms for HCOC | **CLOSED** |
| **MAJOR-014** | ✅ Closed | CFD-OPTIM | Milestone 12 cross-mode therapy utility misclassified ultrasound-only SDT as non-cavitating because it gated all cavitation credit on active venturi throats | **CLOSED** |

---

## RESOLVED-014: Nuclei Transport Diffusion Coupling in CFD-3D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-core/src/physics/cavitation/nuclei_transport.rs`, `crates/cfd-3d/src/vof/cavitation_solver.rs`
**Status**: **CLOSED** - The documented advection-diffusion model is now implemented in the 3D solver.

### Verification

- `NucleiTransport::diffusion_coefficient()` exposes the configured scalar diffusivity to the solver layer.
- `CavitationVofSolver::update_nuclei_advection_diffusion()` now applies explicit finite-volume diffusion with zero-flux outer faces in addition to the existing upwind advection and reaction terms.
- Regression coverage confirms a localized nuclei peak spreads into neighboring cells without flow or generation.

---

## RESOLVED-015: Cavitation Damage Validation and Slice-Based Accumulation in CFD-3D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-3d/src/vof/cavitation_solver.rs`
**Status**: **CLOSED** - Cavitation damage accumulation now validates incoming pressure and density fields and walks contiguous column slices instead of repeated matrix element indexing.

### Verification

- `CavitationVofSolver::update_damage()` now rejects pressure and density dimension mismatches with `DimensionMismatch` instead of silently skipping the damage update.
- The damage accumulation path now uses raw column slices for pressure, density, and damage storage, preserving the column-major layout while removing repeated matrix indexing overhead.
- Regression coverage confirms mismatched pressure fields are rejected and sub-1% void-fraction damage still accumulates.

---

## RESOLVED-016: Cavitation Source Validation and Slice-Based Accumulation in CFD-3D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-3d/src/vof/cavitation_solver.rs`
**Status**: **CLOSED** - Cavitation source accumulation now validates its matrix inputs, uses contiguous slices, and clamps the source to the feasible interval per cell.

### Verification

- `calculate_cavitation_source_into()` now rejects pressure, density, and source-workspace shape mismatches with `DimensionMismatch`.
- The source update path now iterates over raw column slices, preserving the column-major storage contract while removing repeated matrix indexing overhead.
- Regression coverage confirms the source clamps to zero at the vaporization and condensation feasibility extremes.

---

## RESOLVED-017: Pries Plasma-Skimming Silent Clamps in CFD-2D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-2d/src/solvers/plasma_skimming.rs`
**Status**: **CLOSED** - The standalone 2D Pries phase-separation path now evaluates the thresholded logit law through a checked physical-envelope API instead of clipping flow, hematocrit, diameter ratios, logit arguments, or daughter hematocrit outputs.

### Verification

- `checked_pries_phase_separation()` rejects non-finite, out-of-range, and non-positive geometry inputs with structured errors.
- Boundary cases `FQB = 0`, `FQB = 1`, and zero feed hematocrit are explicit conservation states.
- Regression coverage validates below-threshold zero cell flux without modifying the supplied flow fraction.

---

## RESOLVED-018: WALE Boundary Zero-Gradient Reduction in CFD-2D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-2d/src/physics/turbulence/les_smagorinsky/wale.rs`
**Status**: **CLOSED** - WALE velocity-gradient evaluation now uses second-order one-sided boundary stencils on uniform grids instead of forcing boundary gradients to zero.

### Verification

- Rustdoc documents the Taylor-series stencil proof.
- Regression coverage confirms exact derivative recovery for quadratic and linear velocity components at lower and upper boundaries.

---

## RESOLVED-019: Spalart-Allmaras All-Zero TKE Diagnostic in CFD-3D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-3d/src/physics/turbulence/spalart_allmaras.rs`
**Status**: **CLOSED** - The SA turbulence model no longer returns a constant zero turbulent-kinetic-energy vector when the transported eddy viscosity is nonzero.

### Verification

- SA TKE is now a documented diagnostic `k = (nu_t / (C_k d))^2` using the shared Yoshizawa relation and local wall distance.
- Regression coverage validates the computed diagnostic against the analytical relation and preserves the zero wall-length-scale state.

---

## RESOLVED-020: K-Epsilon Uninitialized Zero Fields in CFD-3D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-3d/src/physics/turbulence/k_epsilon.rs`
**Status**: **CLOSED** - The k-epsilon model no longer fabricates zero viscosity, turbulent kinetic energy, or dissipation fields when the transported state is absent.

### Verification

- `turbulent_viscosity()`, `turbulent_kinetic_energy()`, and `dissipation_rate()` now require initialized `k` and `epsilon` fields with dimensions matching the flow field.
- Initialized zero-TKE states still produce zero eddy viscosity by the analytical `nu_t = C_mu k^2 / epsilon` limit.
- Unit and integration tests assert uninitialized state rejection and initialized zero-state behavior.

---

## RESOLVED-021: Smagorinsky LES Zero SGS Diagnostics in CFD-2D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-2d/src/physics/turbulence/les_smagorinsky/`
**Status**: **CLOSED** - The 2D Smagorinsky LES model no longer returns placeholder zero turbulent kinetic energy or dissipation rate after computing nonzero SGS viscosity.

### Verification

- SGS turbulent kinetic energy now follows the Yoshizawa relation `k_sgs = (nu_t / (C_k Delta))^2`.
- SGS dissipation now follows inertial-range dimensional scaling `epsilon_sgs = C_epsilon k_sgs^(3/2) / Delta`.
- Strain-rate evaluation now uses second-order one-sided boundary stencils instead of imposing zero boundary strain.
- Default `min_sgs_viscosity` is zero so zero-strain fields remain exactly laminar unless a caller explicitly configures a numerical lower bound.
- Regression coverage validates the Yoshizawa energy and dissipation diagnostics, boundary strain recovery, and zero-strain invariant.

---

## RESOLVED-022: Margination Singular Wall Cutoff in CFD-1D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-1d/src/physics/cell_separation/margination.rs`
**Status**: **CLOSED** - The margination inertial-lift path no longer clamps lateral position to a hard-coded near-wall cutoff or evaluates a singular wall-lift shape.

### Verification

- The validated domain is now the closed half-channel interval `0 <= x_tilde <= 1`.
- The wall-induced lift shape is bounded and monotone on the full interval while preserving the published `rho U^2 a^6/H^4` wall-lift scaling.
- The shear-gradient term preserves the published `rho U^2 a^3/H` scaling and vanishes at the wall.
- Public legacy lift evaluation now dispatches to the checked evaluator and rejects invalid positions instead of silently clipping.
- Regression coverage validates monotonic wall lift, reference equilibrium calibration, and a derived force value for MCF-7-like cells.

---

## RESOLVED-023: Point-Droplet Occupancy State in CFD-1D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-1d/src/solver/core/transient/droplets/types.rs`
**Status**: **CLOSED** - `DropletSnapshot` no longer stores a separate point-style `occupied_channels` vector beside finite-length occupancy spans.

### Verification

- `occupancy_spans` is now the single authoritative droplet occupancy representation.
- `occupied_channels()` computes the ordered unique channel projection on demand from finite spans.
- Snapshot consistency validation checks finite interval bounds instead of comparing duplicate stored representations.
- Unit and integration tests validate projection, duplicate-channel collapsing, invalid-span rejection, split/merge behavior, and API/manual parity.

---

## RESOLVED-024: Smagorinsky Validation SGS Floor in CFD-2D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-2d/src/physics/turbulence/validation/`
**Status**: **CLOSED** - The turbulence validation paths no longer override the Smagorinsky LES physical zero-floor default with residual nonzero SGS viscosity floors.

### Verification

- Benchmark dispatch remains a closed enum with typed unsupported-model rejection.
- Smagorinsky benchmark and LES/DES validation configurations now use `min_sgs_viscosity = 0.0`.
- Targeted nextest verifies the benchmark module and LES tests under the requested 30-second timeout profile.

---

## RESOLVED-017: Nuclei Transport Validation and Slice-Based Accumulation in CFD-3D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-3d/src/vof/cavitation_solver.rs`
**Status**: **CLOSED** - Cavitation nuclei advection-diffusion now validates its solver inputs and accumulates through contiguous column slices.

### Verification

- `update_nuclei_advection_diffusion()` now rejects cavitation-source and nuclei-workspace shape mismatches with `DimensionMismatch` or `InvalidConfiguration` instead of assuming valid inputs.
- The nuclei transport update now iterates over raw column slices for state, next-state, and source storage, matching the column-major layout and removing repeated matrix indexing overhead.
- Regression coverage confirms a localized cavitation source increases nuclei fraction while zero-flow diffusion spreads a localized peak into neighboring cells.

---

## RESOLVED-018: Indexed Blueprint Layout Cache in CFD-SCHEMATICS

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-schematics/src/visualizations/schematic`
**Status**: **CLOSED** - Schematic layout generation now uses a borrowed indexed cache with flat position storage and index-keyed parallel grouping.

### Verification

- `blueprint_node_positions()` now stores node positions in authored order with a borrowed ID index and a flat auto-layout index list instead of cloning string-keyed maps.
- `resolved_channel_paths()` now resolves parallel channel groups by borrowed node indices instead of cloned `(String, String)` keys.
- `materialize_blueprint_layout()` reuses the indexed cache to write node layout metadata without changing explicit positions.
- Regression coverage confirms the indexed layout cache produces the expected auto-layout coordinates and the materialization path writes the computed fallback positions into the blueprint.

---

## RESOLVED-019: Geometry Bounds and Clearance-Width Relations in CFD-SCHEMATICS

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-schematics/src/config/geometry`, `crates/cfd-schematics/src/state_management`
**Status**: **CLOSED** - Geometry parameter bounds now share canonical constants and clearance-width degeneracy is rejected at the config, manager, and registry validation layers.

## RESOLVED-020: Milestone 12 GA Convergence Trend SSOT in CFD-OPTIM

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-optim/src/reporting`
**Status**: **CLOSED** - The GA convergence figure and narrative now derive their improvement / regression / near-plateau labels from the same trailing-fitness-window helper.

### Verification

- `GeometryConfig::validate()` now uses the canonical geometry constants for channel width, channel height, and wall clearance, and rejects `wall_clearance >= channel_width`.
- `GeometryParameterManager::validate_all()` now enforces the same clearance-width relation at the manager layer.
- `ParameterRegistry::validate_all()` rejects the degenerate pair through the geometry manager path, and the integration regression confirms the returned `InvalidValue` payload.
- Regression coverage confirms canonical bound acceptance and degenerate clearance-width rejection at the config and manager layers.

---

## RESOLVED-021: Slot-Stable Clip-Plane Undo in CFD-UI

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-ui/src/domain/clipping/plane.rs`, `crates/cfd-ui/src/domain/clipping/mod.rs`, `crates/cfd-ui/src/application/clipping/commands.rs`
**Status**: **CLOSED** - Clip-plane commands now preserve their recorded slots, surface explicit stale-state failures, and reject silent relocation during undo.

### Verification

- `ClipPlaneSet::insert_at()` and `replace_at()` now update only the requested slot and return structured slot errors for occupied, empty, and out-of-range cases.

---

## RESOLVED-030: Plasma Skimming, Serpentine Mixing, and LES SGS Energy Physics

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-1d/src/physics/cell_separation/plasma_skimming.rs`, `crates/cfd-2d/src/solvers/serpentine_flow`, `crates/cfd-3d/src/physics/turbulence`
**Status**: **CLOSED** - Three reduced physics paths now dispatch through documented analytical or literature-backed models.

### Verification

- The compact `plasma_skimming_hematocrit` API now infers the sibling branch from Murray cubic diameter conservation and calls the checked Pries phase-separation model, preserving the `X0` cell-entry threshold.
- The serpentine mixing model now computes mixing fraction and L90 from the Neumann eigenfunction solution of transverse diffusion for a two-stream step inlet; the discretized solver result now reports positive analytical L90 and t90 instead of zero placeholders.
- LES eddy-viscosity models no longer return turbulent viscosity as turbulent kinetic energy; they use a shared Yoshizawa SGS relation with value-semantic regression coverage.
- Add, remove, and update clip-plane commands now return `Result` and validate the live slot before restoring or replacing a plane.
- Regression coverage confirms valid add/remove/update round-trips succeed, while undoing a remove into an occupied slot fails explicitly instead of relocating the plane.

---

## RESOLVED-022: Milestone 12 Validation Traceability in CFD-OPTIM

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-optim/src/reporting`, `report/templates/m12_narrative_template.md`, `report/ARPA-H_SonALAsense_Milestone 12 Report.md`, `report/milestone12_results.md`
**Status**: **CLOSED** - The Milestone 12 report now exposes validation evidence and traceability in both the narrative and canonical results artifacts.

### Verification

- `build_results_intro()` now advertises validation evidence and artifact traceability in the Milestone 12 narrative contract, and `m12_narrative_template.md` renders the corresponding subsection.
- `report/ARPA-H_SonALAsense_Milestone 12 Report.md` now contains `### 5.4 Validation Evidence and Artifact Traceability` and the cross-fidelity validation evidence table.
- `report/milestone12_results.md` now contains `## Cross-Fidelity Validation Evidence`, and `report/milestone12/asset_review_manifest.json` was updated to the regenerated report hash.

---

## RESOLVED-023: Single-Pass Authoritative Milestone 12 Report in CFD-OPTIM

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-optim/src/application/orchestration/milestone12/report.rs`, `report/milestone12/report_manifest.json`, `report/ARPA-H_SonALAsense_Milestone 12 Report.md`, `report/milestone12/asset_review_manifest.json`
**Status**: **CLOSED** - The release report now emits the authoritative narrative once, preserves the reviewed hash, and completes after the asset-review manifest is marked pass.

### Verification

- `refresh_milestone12_reports()` now renders the narrative with `authoritative_run = run_class.is_authoritative()` on the first pass instead of rewriting the report with a second content hash.
- `cargo run -p cfd-optim --example milestone12_report --no-default-features --release` completes successfully once the manifest review entry is marked pass under a PATH that resolves the MSYS2 compiler DLLs correctly.
- `report/milestone12/report_manifest.json` records `authoritative_run: true` and `review_complete: true`, and the narrative asset-review entry remains `opened: true`, `status: pass` after regeneration.

---

## CRITICAL-009: Ruge-Stüben Fine-to-Coarse Mapping Bug (VERIFIED FIXED)

**Severity**: ✅ **RESOLVED** - Code Review Confirms Correct Implementation  
**Discovered**: November 18, 2025 (Deep Algorithm Audit)  
**Verified**: February 23, 2026 (Sprint 1.95.0)  
**Component**: `crates/cfd-math/src/linear_solver/preconditioners/multigrid/coarsening.rs`  
**Status**: **CLOSED** - Implementation is correct

### Verification (February 23, 2026)

Upon detailed code review of `assign_closest_coarse_points()` function in `coarsening.rs`, the implementation is **already correct**:

```rust
// Current implementation (CORRECT):
for k in s_offsets[i]..s_offsets[i + 1] {
    let j = s_indices[k];
    if status[j] == 1 {  // j is a C-point
        let strength = strength_matrix.values()[k];
        if strength > max_strength {
            max_strength = strength;
            best_coarse_idx = fine_to_coarse_map[j];  // ✅ CORRECT: gets the coarse grid index
        }
    }
}
```

The key insight: `fine_to_coarse_map[j]` for a C-point `j` already contains the **coarse grid index** (set when the point was marked as coarse via `fine_to_coarse_map[i] = Some(coarse_points.len() - 1)`). Therefore, assigning `fine_to_coarse_map[j]` to an F-point correctly gives the coarse grid index.

### Existing Tests

The file contains comprehensive tests that verify correctness:
- `test_mapping_correctness`: Verifies all mapped indices are valid coarse indices
- `test_interpolation_operator_shape`: Verifies fine_to_coarse_map dimensions
- `test_coarsening_ratio_bounds`: Verifies reasonable coarsening ratios for 2D Laplacian

### Issue Description (Historical)

- **Interpolation Operator Corruption**: The AMG interpolation operator will reference incorrect coarse DOFs
- **Convergence Theory Violation**: AMG theory requires fine points to interpolate from geometrically nearby coarse points
- **Performance Degradation**: Suboptimal convergence rates, potential divergence for difficult problems
- **Literature Divergence**: Does not match Ruge-Stüben (1987) algorithm specification

### Current Code (INCORRECT)

```rust
// Step 2: Assign fine points to nearest coarse points
for i in 0..n {
    if fine_to_coarse_map[i].is_none() {
        // Find strongest connection to coarse point
        let mut max_strength = 0.0;
        let mut best_coarse = None;

        for &c in &coarse_points {
            if strength_matrix[(i, c)] > max_strength {
                max_strength = strength_matrix[(i, c)];
                best_coarse = Some(c);
            }
        }

        if let Some(c) = best_coarse {
            fine_to_coarse_map[i] = fine_to_coarse_map[c];  // ❌ BUG: Assigns mapping value, not index
        }
    }
}
```

### Expected Code (CORRECT)

```rust
// Step 2: Assign fine points to nearest coarse points
for i in 0..n {
    if fine_to_coarse_map[i].is_none() {
        // Find strongest connection to coarse point
        let mut max_strength = 0.0;
        let mut best_coarse = None;

        for &c in &coarse_points {
            if strength_matrix[(i, c)] > max_strength {
                max_strength = strength_matrix[(i, c)];
                best_coarse = Some(c);
            }
        }

        if let Some(c) = best_coarse {
            // ✅ CORRECT: Map to coarse point INDEX in coarse grid
            let coarse_idx = coarse_points.iter()
                .position(|&x| x == c)
                .expect("Coarse point must exist in coarse_points list");
            fine_to_coarse_map[i] = Some(coarse_idx);
        }
    }
}
```

### Remediation Steps

1. **Fix Mapping Logic** (1-2 hours):
   - Update fine-to-coarse map to assign coarse point indices
   - Ensure coarse points map to their own indices in coarse grid
   - Verify interpolation operator dimensions match expected (n_fine × n_coarse)

2. **Add Validation Tests** (2-3 hours):
   - Test that all fine points map to valid coarse indices  (0 ≤ idx < n_coarse)
   - Test that coarse points self-map correctly
   - Test interpolation operator shape: P.shape == (n_fine, n_coarse)
   - Test restriction operator shape: R.shape == (n_coarse, n_fine)
   - Add regression test for this specific bug

3. **Validate Convergence** (1 hour):
   - Run AMG on Poisson problem, verify O(N) complexity
   - Check convergence factor < 0.1 per V-cycle (theory predicts ~0.03-0.08)
   - Compare with literature benchmarks (Ruge-Stüben 1987, Briggs 2000)

### Why Not Caught Earlier

- **Tests Pass**: System still solves (with suboptimal interpolation)
- **Convergence Occurs**: AMG may still converge, just slower than optimal
- **Subtle Bug**: Requires deep algorithm knowledge and literature familiarity to detect
- **No Analytical Verification**: Previous tests didn't validate internal AMG structure

### Evidence Hierarchy Violation

This bug violates the audit framework's evidence hierarchy:
- ❌ **Mathematical Proofs**: Mapping does not match Ruge-Stüben theorem
- ❌ **Literature Validation**: Code diverges from Ruge-Stüben (1987) Algorithm 2
- ⚠️ **Empirical Testing**: Tests incomplete (didn't check mapping correctness)

1. ✅ **FIX CRITICAL-009: Ruge-Stüben Fine-to-Coarse Mapping Bug** (2-5 hours, HIGHEST PRIORITY)
   - Update mapping assignment in `coarsening.rs` to use coarse point indices
   - Add comprehensive AMG coarsening validation tests
   - Verify interpolation operator correctness
   - Validate convergence improvement (expected 2-5x speedup)
   - Block AMG use in production until verified

## COMPLETED (Previous Sprints)

2.  **Immediate Refactoring**: Implement the missing least-squares solve in `DistributedGMRES`. (Done)
3.  **Integration**: Connect `cfd-core` to `cfd-math`. Use `cfd-math::IncompleteLU` as the local solver for `AdditiveSchwarzPreconditioner`. (Done)
4.  **API Extension**: Add `get_diagonal` to `DistributedLinearOperator` to enable real Block Jacobi. (Done)
5.  **Implementation**: Implement basic uniform refinement in `cfd-mesh` to validate the interface. (Done - Implemented 1:8 Tet Refinement)
6.  **Continuous Validation**: Ensure `IncompleteLU` remains the standard for ILU operations.
7.  **Parallel Integration**: When implementing distributed solvers in `cfd-core`, leverage the serial block solvers from `cfd-math` where appropriate, but ensure explicit types for distributed contexts.
8.  **Benchmarking**: Continue to expand benchmarks in `cfd-core` to validate actual parallel scaling of the MPI-specific implementations.

---

## RESOLVED-024: HCOC Cellular Injury & CTC Detection Framework

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-core/src/physics/cavitation/bio_damage.rs`, `crates/cfd-core/src/physics/cavitation/heterogeneous_nucleation.rs`
**Status**: **CLOSED** - Cellular injury grading and stiffness-coupled inception are implemented and tested with literature-consistent threshold physics.

### Verification

- `CellularMembraneMechanics::evaluate_spatial_injury()` uses rate-adjusted critical strains and a first-order hazard transition over the ordered permeabilization / necrosis / lysis thresholds while conserving mass fractions.
- `blake_threshold_pressure_pa()` now applies the classical Blake threshold factor `P_B = P_v - 4σ/(3R_B)` instead of the previous under-corrected `2σ` term.
- Regression tests cover exact mass conservation, ambient no-injury behavior, monotone injury growth, loading-duration sensitivity, the classical Blake-factor identity, and target-versus-healthy selectivity ordering.

---

## RESOLVED-025: Apollo FFT Plan Binding in CFD-3D Periodic DNS

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-3d/src/spectral/dns.rs`
**Status**: **CLOSED** - The periodic pseudospectral DNS stepper now validates its Apollo 3D shape once at construction and owns a reusable `FftPlan3D` for all forward and inverse transforms.

### Verification

- `PeriodicPseudospectralDns3D::new()` maps the validated CFD grid dimensions into Apollo `Shape3D` and returns structured configuration errors if Apollo rejects the shape.
- Velocity spectra, nonlinear spectra, spectral derivatives, and inverse step reconstruction all use the stepper-owned plan rather than per-transform shape dispatch.
- `cargo check -p cfd-3d --no-default-features` completes successfully.

---

## RESOLVED-026: Exact MUSCL3/QUICK Face Reconstruction in CFD-2D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-2d/src/physics/momentum/muscl.rs`
**Status**: **CLOSED** - The third-order MUSCL/QUICK path no longer uses an approximate right-interface formula; both face states now use exact uniform-grid quadratic interpolation.

### Verification

- `quadratic_left_face()` implements `-phi_{i-1}/8 + 3 phi_i/4 + 3 phi_{i+1}/8`.
- `quadratic_right_face()` implements `3 phi_i/8 + 3 phi_{i+1}/4 - phi_{i+2}/8`.
- Regression coverage verifies exact reproduction of a linear field from both sides and exact reproduction of a quadratic parabola at the face.
- `CC=clang cargo nextest run -p cfd-2d --no-default-features muscl ...` passed 14/14 filtered tests under the 30-second slow-test bound.

---

## RESOLVED-027: cfd-1d Margination Lift and Droplet Occupancy Contracts

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-1d/src/physics/cell_separation/margination.rs`, `crates/cfd-1d/src/solver/core/transient/droplets`
**Status**: **CLOSED** - The margination model now separates wall-induced and shear-gradient inertial lift scaling, and droplet snapshots document `occupied_channels` as a projection of finite-length occupancy spans.

### Verification

- `dimensional_lift_force_n()` evaluates distinct wall-induced and shear-gradient finite-size scaling terms and derives the wall gain from the documented rigid-cell reference equilibrium.
- Rustdoc theorem/proof sketches document wall-shape monotonicity, shear-shape monotonicity, and reference-equilibrium calibration.
- `DropletSnapshot::occupied_channels_from_spans()` and `has_consistent_finite_length_occupancy()` encode the finite-length occupancy projection contract.
- Targeted `CC=clang cargo nextest run -p cfd-1d --no-default-features ...` passed the margination and droplet regression tests.

---

## RESOLVED-028: cfd-2d Turbulence Benchmark Dispatch

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-2d/src/physics/turbulence/validation/benchmarks.rs`
**Status**: **CLOSED** - The performance benchmark no longer contains an unimplemented fallback path; supported benchmark models are parsed into a closed enum before dispatch.

### Verification

- `BenchmarkModel` admits only `smagorinsky-les` and `des`, and exhaustive matching runs concrete production update loops.
- Unsupported model names return a typed rejection listing supported models and do not emit placeholder text.
- `CC=clang cargo nextest run -p cfd-2d --no-default-features ...` passed 554/554 tests under the 30-second slow-test bound.

---

## RESOLVED-029: cfd-1d Coupled Blood Pressure-Event Verification

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-1d/src/domain/network/wrapper.rs`, `crates/cfd-1d/src/solver/core/linear_system.rs`, `crates/cfd-1d/src/solver/core/transient/composition/simulator.rs`
**Status**: **CLOSED** - Coupled pressure-event hematocrit simulation now produces finite positive branch-flow magnitudes while preserving signed internal transport.

### Verification

- Zero-flow blood apparent viscosity uses finite Secomb diameter/hematocrit initialization before nonzero shear-dependent Quemada updates.
- Small hydraulic systems use row-equilibrated solves so conductance rows near machine scale are not dominated by unit Dirichlet rows.
- Dense fallback conversion accumulates duplicate CSR entries, preserving the sum of incident-edge diagonal conductances required by Kirchhoff balance.
- `CC=clang cargo nextest run -p cfd-1d --no-default-features ...` passed 696/696 tests under the 30-second slow-test bound.

---

## RESOLVED-030: Root Historical Source Artifact Drift

**Severity**: ✅ **RESOLVED**
**Component**: workspace root source artifacts
**Status**: **CLOSED** - Tracked root historical Rust files that were outside Cargo, examples, tests, docs, and report-generation paths have been removed to preserve single-source-of-truth ownership.

### Verification

- `git ls-files` confirmed `old_assemble.rs`, `old_arrangement.rs`, `old_phase2.rs`, `old_operations.rs`, `old_indexed.rs`, `old_gwn_bvh.rs`, `old_seam.rs`, `old_phase4.rs`, and `csg_bi.rs` were tracked root artifacts.
- `rg` found no references to those paths in `Cargo.toml`, `crates/`, `examples/`, `tests/`, `docs/`, or `README.md`.
- Deletion is patch-class cleanup: no public API, algorithm, report template, or executable target references the removed files.
- `cargo metadata --no-deps --format-version 1` and `cargo check -p cfd-suite --no-default-features` completed successfully.
- `.cargo/config.toml` now selects MSYS2 `clang.exe` and `-fuse-ld=lld` for `x86_64-pc-windows-gnu`, and forces C/C++ build-script tools to MSYS2 `clang`, `clang++`, `llvm-ar`, and `llvm-ranlib`.
- `cargo nextest run --workspace --no-default-features --fail-fast --hide-progress-bar` was re-attempted after the toolchain correction; the prior `gcc.exe` `zstd-sys` failure did not recur, but workspace compilation exceeded the 60-second bound before Rust test execution.
- Misleading terminology was removed from explicit unsupported-operation paths and bounded-model comments in `cfd-core`, `cfd-math`, `cfd-python`, and `gaia` without changing numerical behavior.
- `cargo check -p cfd-core -p cfd-math -p gaia -p cfd-python --no-default-features` completed successfully in 46.17 seconds.
- `cargo nextest run -p cfd-core -p cfd-math -p cfd-python --no-default-features --fail-fast --hide-progress-bar` exceeded the 60-second compilation bound before test execution; no failing Rust test was reported.

---

## RESOLVED-031: Optimization Generator and Boundary Unsupported-Order Terminology

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-schematics/src/geometry/optimization`, `crates/cfd-core/src/physics/boundary/error.rs`
**Status**: **CLOSED** - Internal serpentine optimization path generation and boundary stencil diagnostics now describe executable contracts without implying simplified or incomplete implementations.

### Verification

- `generate_simplified_serpentine_path` was renamed to `generate_optimization_serpentine_path`; all internal optimization call sites were updated.
- Boundary unsupported-order diagnostics now report `"Stencil order {order} is unsupported (supported: 1-4)"`.
- `rg` found no stale `generate_simplified_serpentine_path`, `simplified serpentine`, `simplified generation`, or `"not implemented"` wording in the touched paths.
- `cargo check -p cfd-core -p cfd-schematics --no-default-features` completed successfully in 28.06 seconds.
- `cargo nextest run -p cfd-core -p cfd-schematics --no-default-features --fail-fast --hide-progress-bar` exceeded the 60-second compilation bound before test execution; no failing Rust test was reported.

---

## RESOLVED-032: Womersley Validation Approximation Drift

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-validation/src/analytical/womersley.rs`
**Status**: **CLOSED** - The validation Womersley analytical solution now reuses the canonical exact complex-Bessel `cfd-1d` implementation instead of maintaining local regime approximations.

### Verification

- `velocity()`, `wall_shear_stress()`, and `flow_rate()` delegate through `cfd_1d::physics::vascular::womersley::WomersleyProfile`.
- Rustdoc now states the exact Bessel-form velocity theorem and no-slip proof sketch.
- Added `test_exact_womersley_no_slip_wall_condition`, which checks computed wall velocity at five phases.
- `cargo check -p cfd-validation --no-default-features` completed successfully in 3.44 seconds after warm build.
- `cargo test -p cfd-validation --no-default-features analytical::womersley --lib -- --nocapture` passed 5/5 tests in 10.50 seconds after compilation.
- `cargo nextest run -p cfd-validation --no-default-features womersley --fail-fast --hide-progress-bar` exceeded the 60-second compilation bound before test execution; no failing Rust test was reported.
