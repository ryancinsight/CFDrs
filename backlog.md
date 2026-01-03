# Project Backlog

## Strategy
- **Prioritization**: Math Correctness > Critical Bugs > Test Coverage > Features > Optimization
- **Workflow**: gap_audit -> backlog -> checklist -> implementation

## Sprint 1.84.0: CFD-1D Numerical Audit (Immediate)

### High Priority
- [x] **1D-001**: Fix RectangularChannelModel laminar resistance (must be constant).
- [x] **1D-002**: Fix DarcyWeisbachModel units and non-linear consistency.
- [x] **1D-003**: Correct MatrixAssembler units documentation.
- [x] **1D-004**: Add Mach number and entrance length validation.
- [x] **1D-005**: Integrate invariant validation into simulation loop.

## Sprint 1.86.0: CFD Validation & Benchmarking Audit (Immediate)

### High Priority
- [x] **VAL-001**: Fix Richardson extrapolation grid ordering and ratio calculation.
- [x] **VAL-002**: Correct Laplace equation solver sign error in MMS tests.
- [x] **BENCH-001**: Fix atomic ordering in memory profiling (Acquire/Release).
- [x] **BENCH-002**: Improve memory metric display accuracy (peak vs total).
- [x] **BENCH-003**: Fix scaling analysis thread pool initialization panic.
- [x] **BENCH-004**: Align benchmarking suite with unified BenchmarkConfig.
- [x] **TEST-001**: Resolve multi-physics coupled interaction test failures.
- [x] **TEST-002**: Fix markdown report generation test failures.

## Sprint 1.87.0+ Items (Planned)

#### Sprint 1.84.0: Performance Optimization
- [ ] Cache AMG preconditioner (2-5x pressure solve speedup)
- [ ] Optional SIMD optimization for CG

#### Sprint 1.85-1.87: Test Coverage Expansion
- [ ] Target: 8.82% → 50% → 80%
- [ ] Focus: CFD-Math critical paths, CFD-2D pressure-velocity, integration tests

#### Sprint 1.88.0+: Strategic Enhancements
- [ ] GPU integration decision (implement or defer)
- [ ] MPI scaling validation (empirical benchmarks)
