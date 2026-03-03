# Project Checklist

## Sprint 1.95.0: DEEP PHYSICS AUDIT & PIPELINE VALIDATION ✅ COMPLETE

**Status**: ✅ COMPLETE — R-P equation corrected, Smagorinsky test fixed, SDT pipeline validated (425K sweep), mesh-export pipeline verified (5 designs with STL/OpenFOAM/SVG/JSON).

### Physics Fixes
- [x] `cfd-core`: Rayleigh-Plesset `bubble_acceleration()` — added polytropic gas pressure `p_g0·(R_0/R)^{3κ}`, removed double surface tension
- [x] `cfd-3d`: Smagorinsky test — use `with_filter_width(cs, Δ, Δ, Δ)` matching grid spacing instead of `new(cs)` (filter_width=1.0)

### Code Quality Scan
- [x] No `unimplemented!()`, `todo!()`, or placeholder macros in any crate
- [x] All workspace tests pass (`cargo test --workspace --no-default-features --exclude cfd-python`)

### Pipeline Validation
- [x] `sdt_therapy` example: 425K candidate sweep, exit 0, top-5 CCT/CIFX designs (σ<0, FDA compliant)
- [x] `sdt_pipeline --features mesh-export`: 5 designs → STL + OpenFOAM polyMesh + SVG + JSON, exit 0
- [x] Report figures verified: all 10 referenced figures exist in `report/figures/`

---

## Sprint 1.90.0: CFD CRATES COMPREHENSIVE AUDIT & ENHANCEMENT ✅ COMPLETE

**Status**: ✅ COMPLETE — All 5 crates audited, theorem-documented, adversarially tested. All tests pass.

### Structural Changes
- [x] `cfd-math`: `performance_monitor.rs` → `diagnostics/performance_monitor.rs` (SRP)
- [x] `cfd-math`: `diagnostics/mod.rs` with Adaptive Threshold Optimality theorem + backward-compat shim
- [x] `cfd-core`: Removed empty orphan `management/neuronavigation/`
- [x] `cfd-1d`: `resistance/traits.rs` converted from misleading redirect comment to proper `pub use` re-export (SSOT + DIP)
- [x] `cfd-3d`: All 9 bare `.unwrap()` in `lib.rs` tests → `.expect("invariant message")` / `Result<()`
- [x] `cfd-3d`: `src/tests/adversarial_tests.rs` created and wired into module

### Theorem Documentation (28 formal theorems added)
- [x] `cfd-math`: Taylor Remainder, Spectral, SBP | Gauss-Legendre, Composite, Fubini | CSR SpMV, Cholesky, Parallel | Runge, Lebesgue, Cubic Spline | Dahlquist, RK4, CFL, RKC
- [x] `cfd-core/physics`: NS Existence (Leray), Well-Posedness (Lax), Frame Invariance (Noll), Buckingham Π
- [x] `cfd-1d/resistance`: Hagen-Poiseuille, Rectangular Duct, Dean Number, Kirchhoff's Laws, Knudsen Regime
- [x] `cfd-2d/solvers`: Lax-Richtmyer, Gauss Divergence, BGK+Chapman-Enskog, SIMPLE, PISO, Fixed-Point Contraction
- [x] `cfd-3d/vof`: PLIC Conservation, Accuracy O(h²), Young-Laplace, Maximum Principle
- [x] `cfd-3d/spectral`: Parseval, Aliasing/Nyquist (3/2 dealiasing), Gauss-Lobatto, Exponential Convergence

### Adversarial Tests (21 new tests)
- [x] `cfd-math`: ill-conditioned (κ=10^14), zero-RHS, dimension mismatch, n=500 SPD stress, NaN-safety, convergence history, sparse builder bounds
- [x] `cfd-math`: interpolation extrapolation rejection, single-point, duplicate nodes, quadratic exactness
- [x] `cfd-3d`: VOF conservation/bounds, Parseval round-trip, energy localization, aliasing, FEM config invariants, level set

### Verification Results
- [x] `cargo check` exit 0 on all 5 crates
- [x] `cargo test -p cfd-math` → all pass   |   `cargo test -p cfd-core` → all pass
- [x] `cargo test -p cfd-1d` → all pass   |   `cargo test -p cfd-2d` → all pass   |   `cargo test -p cfd-3d` → all pass

---

## Sprint 1.91.0: ARCHITECTURAL PURITY PHASE 2 ✅ COMPLETE
- [x] **WENO Deduplication**: `cfd-math/src/high_order/weno.rs` vs `src/spatial/weno.rs` (SSOT established, `spatial` deleted)
- [x] Update all `cfd_core::domain` → `cfd_core::geometry` references (REST-002) (Verified completed)
- [x] Remove simplified solvers from `cfd-core` (REST-004) (Replaced simplified blood/cavitation with exact math, removed obsolete files)
- [x] CSG Boolean exact predicates (cfd-mesh) (Verified exact GWN and arithmetic predicates exist)

---

## Sprint 1.86.0 Status: ✅ COMPLETE

- [x] Resolve Richardson extrapolation validation errors
- [x] Fix benchmarking framework thread-safety and accuracy
- [x] Verify multi-physics coupled validation tests
- [x] Ensure all benchmark suites are executable

### Validation & Benchmarking Fixes:
1. **Richardson Extrapolation** - Fixed grid ordering reversal and ratio calculation in `validation.rs`.
2. **MMS Tests** - Corrected sign error in Laplace equation solver.
3. **Memory Profiling** - Replaced `Relaxed` with `Acquire/Release` ordering in `memory.rs`.
4. **Scaling Analysis** - Fixed thread pool reinitialization panic and unified `BenchmarkConfig`.
5. **Multi-Physics** - Resolved MHD property test and coupled physics interaction failures.
6. **Reporting** - Fixed markdown report generation tests.

---

## Sprint 1.84.0: CFD-1D NUMERICAL AUDIT (IMMEDIATE)
**Status**: ✅ COMPLETE - Audited 1D CFD Physics and Numerical Methods

### CFD-1D Numerical Audit (December 20, 2025)
A deep audit of the 1D CFD implementation has identified and resolved several mathematical and physical correctness issues:

- ✅ **FIXED-1D-001**: RectangularChannelModel uses constant resistance for laminar flow.
- ✅ **FIXED-1D-002**: DarcyWeisbachModel corrected units and split into linear/quadratic coefficients.
- ✅ **FIXED-Gap**: MatrixAssembler units documentation corrected.
- ✅ **FIXED-Gap**: Mach number validation (Ma < 0.3) implemented.
- ✅ **FIXED-Gap**: Entrance length validation (L/Dh > 10) implemented.

### Sprint 1.84.0 Tasks
- [x] **Priority 1**: Fix RectangularChannelModel Math ✅
- [x] **Priority 2**: Fix DarcyWeisbachModel Math and Units ✅
- [x] **Priority 3**: Correct MatrixAssembler and Solver Documentation ✅
- [x] **Priority 4**: Add Physical Invariant Validation ✅
- [x] **Priority 5**: Verification ✅


---

## Sprint 1.87.0: ARCHITECTURAL PURITY (IMMEDIATE)
**Status**: 🚧 IN PROGRESS - Restructuring into Deep Vertical Tree

### Strategic Restructuring Tasks
- [x] **Phase 1: Audit & Cleanup**
    - [x] Delete `reynolds_stress.rs.backup` and `backend_example.rs` ✅
    - [x] Consolidate `cfd-math/src/vectorization/` into `cfd-math/src/simd/vectorization.rs` ✅
    - [x] Remove obsolete `vectorization` module from `cfd-math/src/lib.rs` ✅
    - [ ] Audit `cfd-core` for remaining "Potemkin" stubs
    - [ ] Verify `cfd-math` duplication (WENO)
- [ ] **Phase 2: Physics Consolidation (REST-001)**
    - [x] Verify `fluid`, `material`, `boundary`, `constants` are in `cfd-core/src/physics/` ✅
    - [ ] Update re-exports and documentation

**TODO(MEDIUM): Update physics module re-exports and documentation after consolidation. Dependencies: REST-001 completion. Reference: SOLID principles for clean module boundaries.**
- [ ] **Phase 3: Geometry Consolidation (REST-002)**
    - [x] Rename `domain` to `geometry` in `cfd-core` ✅
    - [x] Move `rhie_chow.rs` to `physics/fluid_dynamics` and remove `interpolation` ✅
    - [ ] Update all external references to `cfd_core::domain` to `cfd_core::geometry`

**TODO(MEDIUM): Update all external references from cfd_core::domain to cfd_core::geometry post-renaming. Dependencies: REST-002 completion. Mathematical foundation: Maintain API consistency for bounded contexts.**
- [ ] **Phase 4: Compute Consolidation (REST-003)**
    - [x] Move `gpu`, `mpi`, `simd` into `cfd-core/src/compute/` ✅
- [ ] **Phase 5: Solver Decoupling (REST-004)**
    - [ ] Remove simplified solvers from `cfd-core`
    - [ ] Migrate `NumericalMethodsService` to appropriate abstraction level
- [ ] **Phase 6: Duplicate Elimination (REST-005)**
    - [x] Eliminate SIMD duplication in `cfd-math` ✅
    - [ ] Eliminate WENO duplication in `cfd-math` (src/high_order/weno.rs vs src/spatial/weno.rs)

**TODO(HIGH): Complete duplicate elimination in cfd-math (WENO, SIMD already done). Dependencies: REST-005. Reference: DRY principle and performance optimization.**

**TODO(MEDIUM): Complete architectural restructuring to deep vertical tree architecture per gap_audit.md REST-001 through REST-005. Dependencies: Phase 1 must complete before Phase 2-6. Estimated effort: 15-20 hours. Mathematical foundation: SOLID principles and bounded contexts for numerical methods.**

---

## Audit Documents Generated

1. **SPRINT_AUDIT_2025_11_18.md** - Comprehensive stepwise audit (19,772 chars)
2. **gap_audit.md** - Updated with CRITICAL-009 detailed analysis
3. **ENHANCEMENT_PLAN_SPRINTS_1.83-1.90.md** - 8-sprint roadmap (24,584 chars)
4. **AUDIT_SUMMARY_2025_11_18.md** - Executive summary with action items

---

## Future Sprints (Post-1.83.0)

**Sprint 1.84.0**: Performance Optimization (6-8 hours)
- Cache AMG preconditioner (2-5x pressure solve speedup)
- Optional SIMD optimization for CG

**Sprints 1.85-1.87**: Test Coverage Expansion (15-20 hours)
- Target: 8.82% → 50% → 80%
- Focus: CFD-Math critical paths, CFD-2D pressure-velocity, integration tests

**Sprint 1.88.0: CFD-MESH ROBUSTNESS AUDIT & REFACTOR**
**Status**: ✅ COMPLETE - Audited Mesh Processing (CSG, BSP, Welding)

### CFD-Mesh Geometric Audit
A deep audit of the geometric algorithms identified heuristics requiring formal mathematical resolution:
- ✅ **AUDITED**: `csg/boolean.rs` relies on quantized heuristic containment.
- ✅ **AUDITED**: `csg/bsp.rs` splitting uses greedy heuristic with 1e-5 tolerance.
- ✅ **AUDITED**: `welding/snap.rs` lacks topology-preserving invariant checks.

### Tasks
- [ ] **Priority 1**: Refactor CSG Booleans to Exact Predicates (EMBER-style)
- [ ] **Priority 2**: Implement Plane-based BSP Tree exact classifications
- [ ] **Priority 3**: Upgrade SnappingGrid to topological-preserving welding

**Sprint 1.89.0+**: Strategic Enhancements
- GPU integration decision (implement or defer)
- MPI scaling validation (empirical benchmarks)
- Advanced AMG coarsening strategies (PMIS, HMIS, Falgout)

---

## Production Readiness

### ✅ APPROVED (Non-AMG)
- GMRES, BiCGSTAB, CG solvers (mathematically verified)
- Pressure correction SIMPLE/PISO (Patankar-correct)
- Momentum solvers and turbulence models

### ❌ BLOCKED
- AMG Preconditioner (until CRITICAL-009 fixed - Sprint 1.83.0)

### ⚠️ GAPS (Non-Blocking)
- Test coverage (improvement plan in place)
- Dead code cleanup (feature-gating recommended)

---

**Next Action**: Begin Sprint 1.83.0 CRITICAL-009 remediation (estimated 2-5 hours to fix)

**Timeline to Production-Ready AMG**: 8-12 hours (Sprint 1.83.0 completion)

**TODO(HIGH): Implement comprehensive performance benchmarking framework - strong/weak scaling analysis, regression detection, communication profiling per ASME V&V 20-2009. Dependencies: AMG convergence validation. Estimated effort: 8-12 hours.**
