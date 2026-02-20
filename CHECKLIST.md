# Project Checklist

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

**Sprint 1.88.0+**: Strategic Enhancements
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
