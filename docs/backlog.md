# CFD Suite - Technical Backlog (SSOT)

## Sprint 1.28.0-CONVERGENCE-REMEDIATION

### Critical Priority (P0) - BLOCKING CONVERGENCE
- [ ] **FIX-COMPILATION**: 17 errors in pipe_flow_1d example
  - **Impact**: Examples non-functional, API mismatches
  - **Dependencies**: None
  - **Risk**: HIGH - Blocks validation workflows
  - **ETA**: 2h
  - **Acceptance**: `cargo test --example pipe_flow_1d --no-run` succeeds

- [ ] **REDUCE-CLIPPY**: 756 → <100 warnings per standards
  - **Impact**: Static analysis debt per IEEE TSE 2022 safety requirements
  - **Dependencies**: Fix compilation first
  - **Risk**: MEDIUM - Technical debt accumulation
  - **ETA**: 4h
  - **Acceptance**: `cargo clippy --workspace` shows <100 warnings

### High Priority (P1) - INFRASTRUCTURE
- [ ] **DOC-STRUCTURE**: Reorganize per standards (Phase 1 requirement)
  - **Tasks**:
    - Move CHECKLIST.md → docs/checklist.md
    - Move PRD.md → docs/prd.md  
    - Establish docs/ as canonical location
  - **Dependencies**: None
  - **Risk**: LOW - Organizational change
  - **ETA**: 1h
  - **Acceptance**: All docs in docs/ directory

- [ ] **MODULE-AUDIT**: Validate <400 lines per module (Rust forums standard)
  - **Impact**: SOC/modularity/extensibility per SOLID principles
  - **Dependencies**: Doc structure complete
  - **Risk**: MEDIUM - May require refactoring
  - **ETA**: 3h
  - **Acceptance**: All modules <400 lines, proper domain separation

### Medium Priority (P2) - VALIDATION
- [ ] **TEST-RUNTIME**: Ensure <30s per docs/checklist.md requirement
  - **Impact**: CI/CD efficiency, developer productivity
  - **Dependencies**: Compilation fixes
  - **Risk**: LOW - Test optimization
  - **ETA**: 2h
  - **Acceptance**: `cargo test --workspace` completes <30s

- [ ] **PHYSICS-VALIDATION**: Momentum solver accuracy per Chapman-Enskog
  - **Impact**: CFD correctness per literature standards
  - **Dependencies**: Test fixes
  - **Risk**: MEDIUM - May reveal solver issues
  - **ETA**: 4h
  - **Acceptance**: Literature validation <2% error

### Low Priority (P3) - OPTIMIZATION (DEFERRED POST-CONVERGENCE)
- [ ] **BENCHMARKS**: Criterion integration (deferred per standards)
- [ ] **NO_STD**: Embedded CFD support (deferred per standards)
- [ ] **SIMD-OPTIMIZATION**: AVX2/NEON tuning (deferred per standards)

## Technical Debt Inventory

### Current State (Sprint 1.27.0)
- ✅ **Build Quality**: Zero compilation warnings achieved
- ✅ **Solver Physics**: Momentum equation properly implemented  
- ✅ **Example Infrastructure**: Partially restored (spectral_3d fixed)
- ⚠️ **Static Analysis**: 756 clippy warnings (vs target <100)
- ❌ **Documentation**: Missing SSOT artifacts (backlog.md, proper structure)

### Risk Assessment
- **CRITICAL**: API incompatibilities blocking examples
- **HIGH**: Clippy warning accumulation (safety implications per IEEE)
- **MEDIUM**: Module size violations (maintenance burden)
- **LOW**: Test performance (CI efficiency)

## Sprint Retrospective Framework

### Definition of Done
1. All P0 items completed and validated
2. Build/test/clippy metrics within standards
3. Documentation structure per requirements
4. Technical debt reduced, not increased

### Success Metrics
- Compilation errors: 17 → 0
- Clippy warnings: 756 → <100  
- Module violations: TBD → 0
- Test runtime: TBD → <30s

### Iteration Plan
- **Sprint 1.28.0**: Convergence remediation (this sprint)
- **Sprint 1.29.0**: Performance validation vs literature
- **Sprint 1.30.0**: Architecture optimization (if converged)

## Dependencies & Constraints

### External Dependencies
- Rust toolchain: 1.70+ (MSRV)
- Clippy: Built-in static analysis
- Criterion: Benchmarking (deferred)

### Internal Dependencies
- cfd-core: Abstractions layer
- cfd-validation: Literature benchmarks
- Examples: API demonstration

### Resource Constraints
- Development time: 8-12h per sprint
- CI/CD budget: <5min build/test cycle
- Memory: Efficient patterns required

## Notes
- This backlog serves as SSOT per requirements
- Updated every sprint per 3-sprint ADR/SRS cycle
- Priorities aligned with convergence requirements
- Defers optimization until core completion per standards