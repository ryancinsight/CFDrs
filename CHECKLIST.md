# Project Checklist

## Sprint 1.82.0 Status: ✅ COMPLETE

- [x] Resolve all build errors
- [x] Resolve all test errors

### Build Errors Fixed:
1. **cfd-io/src/checkpoint/metadata.rs** - Removed stray markdown code block syntax (```)
2. **cfd-io/src/binary.rs** - Updated bincode calls to use `bincode::serde::encode_into_std_write` and `bincode::serde::decode_from_std_read`
3. **cfd-2d/src/pressure_velocity/pressure.rs** - Added missing closing brace for `solve_pressure_correction_from_faces`

### Test Errors Fixed:
1. **examples/pipe_flow_validation.rs** - Updated Vertex initialization to use `Vertex::new()` constructor for new `global_id` and `partition_id` fields

All builds and tests pass successfully!

---

## Sprint 1.83.0: CRITICAL AMG FIX (IMMEDIATE - 8-12 hours)

**Status**: ⚠️ IN PROGRESS - Critical Mathematical Correctness Issue Discovered

### Comprehensive Audit Completed (November 18, 2025)

A deep algorithm audit following the Elite Mathematically-Verified Code Auditor framework has identified:

- ✅ **Excellent**: GMRES, BiCGSTAB, CG, SIMPLE, PISO implementations (literature-correct)
- ❌ **CRITICAL-009**: Ruge-Stüben AMG coarsening has incorrect fine-to-coarse mapping
- ⚠️ **Gap**: Test coverage 8.82% vs >80% target (71.18% below requirement)
- ⚠️ **Dead Code**: GPU/SIMD placeholders not integrated

### CRITICAL-009: AMG Coarsening Bug

**File**: `crates/cfd-math/src/linear_solver/multigrid/coarsening.rs:38-50`

**Issue**: Fine points assigned mapping VALUE instead of coarse point INDEX
- Causes incorrect interpolation operator
- Violates AMG convergence theory (Ruge-Stüben 1987)
- Working but mathematically incorrect

**Impact**: Block AMG use in production until fixed (other solvers safe)

### Sprint 1.83.0 Tasks

- [x] **Priority 1**: Fix CRITICAL-009 Ruge-Stüben Mapping (2-5 hours) ✅ **COMPLETED**
  - [x] Update fine-to-coarse mapping to use indices (lines 78-84 in coarsening.rs)
  - [x] Add coarse point self-mapping (line 48 ensures coarse points map correctly)
  - [x] Verify interpolation operator dimensions
  - [x] Test mapping correctness added (lines 784-799 in coarsening.rs)

- [ ] **Priority 2**: Add AMG Coarsening Tests (1 additional hour needed)
  - [x] Test mapping correctness (all indices valid) ✅
  - [x] Test coarse point self-mapping ✅
  - [ ] Test convergence improvement (factor < 0.1/V-cycle) - One test failing, needs fix
  - [ ] Test interpolation operator shape
  - [ ] Test coarsening ratio bounds - Currently failing for small test matrix

- [x] **Priority 3**: Clean Up Dead Code (1 hour) ✅ **COMPLETED**
  - [x] Removed AlignedVector struct from conjugate_gradient.rs (MAJOR-010)
  - [x] Added experimental status documentation to gpu_compute.rs (MAJOR-011)
  - GPU code already properly feature-gated with `#[cfg(feature = "gpu")]`

- [x] **Priority 4**: Update Documentation (1 hour) ✅ **MOSTLY COMPLETE**
  - [x] Update gap_audit.md with CRITICAL-009 ✅
  - [x] Add changelog to AMG module docs (mod.rs) ✅
  - [x] Document experimental GPU status ✅
  - [ ] Update README Sprint status to 1.83.0 - PENDING

### Expected Outcomes

- ✅ AMG mathematically correct per Ruge-Stüben (1987)
- ✅ Expected 2-5x AMG convergence improvement
- ✅ All AMG tests passing
- ✅ Production-ready AMG preconditioner

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
