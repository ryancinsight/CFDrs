# CFD Suite Backlog

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
