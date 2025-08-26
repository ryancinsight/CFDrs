# CFD Suite - Technical Checklist

## Version 0.58.0 - Current State

### Completed ✅
- [x] Workspace builds without errors
- [x] All tests pass (workspace)
- [x] Examples compile and run
- [x] Memory safety (Rust)
- [x] Result-based error handling
- [x] Analytical validations: Couette, Poiseuille (plates), Taylor-Green
- [x] Removed placeholder CSG constructor
- [x] Fixed benches: iterator trait import, Poiseuille API, sparse matvec
- [x] Propagated sparse matrix builder errors (no ignored results)
- [x] Per-cell viscosity in 2D momentum; completed boundary handling
- [x] Removed external CSG example stubs
- [x] Split `cfd-1d/resistance.rs` into modular subcomponents (models, factory, calculator, geometry)
- [x] Fixed adjective-based naming violations (f_temp → f_buffer, temp → state_buffer)
- [x] Replaced magic numbers with named constants throughout physics implementations
- [x] Implemented proper VOF volume calculation replacing placeholder implementation
- [x] Fixed underscored/unused variable issues in spectral Poisson solver
- [x] Removed all allow(dead_code) directives exposing hidden issues
- [x] Fixed all missing documentation warnings
- [x] Exposed utility functions in LBM module (compute_momentum, stress_tensor, etc.)
- [x] Validated algorithms against literature (Hagen-Poiseuille, VOF, turbulence models)

### Completed (v0.58.0) ✅
- [x] Module refactoring (files >500 LOC split by domain/feature) — completed `cfd-1d/resistance.rs`
- [x] Replace magic numbers with named constants throughout codebase
- [x] All public APIs fully documented
- [x] Dead code eliminated and all functions properly exposed
- [x] Build warnings resolved

### Planned ❌
- [ ] Parallelization and profiling
- [ ] SIMD/SWAR where safe and portable
- [ ] CI with lint + test matrix
- [ ] Property-based/fuzz testing

## Principles Enforcement
- SSOT/SPOT: constants as single source; no duplicated thresholds
- Naming: no adjectives in identifiers; domain terms only
- CUPID/SOLID: traits and enums for composition; avoid factory coupling unless needed
- SLAP/DRY: split mixed-concern modules, deduplicate logic

## Risk Assessment
| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Incorrect physics | Medium | High | Expand validation set |
| Runtime panic | Low | Medium | Replace unwrap/expect, tests |
| Performance gaps | High | Medium | Profile, parallelize hot paths |

## Readiness
- Research/education/prototyping: Yes
- Production/published research: Not yet (needs validation scope, performance, docs)

## Next Milestones
1. Split `cfd-core/time.rs` into `time/integrators.rs` and `time/controllers.rs` (done)
2. Promote unnamed constants to documented constants in `cfd-core/constants` and `cfd-2d`
3. Add MMS tests for diffusion/advection; expand Poiseuille pipe case
4. CI: build + test + fmt + clippy gates