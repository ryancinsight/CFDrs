# Product Requirements Document

## CFD Suite v0.57.4

### Product Classification
Research software (not production)

### Current Capabilities

#### Functional
- Workspace compiles and runs without errors
- Core algorithms (FVM/FDM/PISO/VOF/spectral) present
- Analytical validations: Couette, Poiseuille (plates), Taylor-Green
- Test suite covers numerical methods and integrations
- Modular architecture with proper domain separation

#### Non-Functional
- Trait-based, domain-structured crates
- Result-based error handling
- Examples and benches compile
- Clean naming conventions (no adjectives)
- Constants centralized (SSOT/SPOT)
- Modules properly sized (< 500 LOC target, mostly achieved)

### Limitations
- Validation coverage limited to selected cases
- Performance and parallelism deferred
- Some documentation incomplete
- Two modules still exceed 500 LOC

### Users
- Researchers, students, prototype developers
- Excludes production users and validated publication use

### Development Roadmap (near term)
1) Structure
- [x] Split `cfd-core/time.rs` into `time/integrators.rs`, `time/controllers.rs`
- [x] Split `cfd-1d/resistance.rs` into modular substructure
- [x] Consolidate hydraulics constants; ban magic numbers
- [ ] Split remaining large modules (`cfd-3d/level_set.rs`, `cfd-validation/numerical_validation.rs`)

2) Code Quality
- [x] Replace adjective-based variable names
- [x] Extract magic numbers to named constants
- [ ] Complete documentation for all public items
- [ ] Add comprehensive error context

3) Validation
- [ ] Add MMS for diffusion/advection
- [ ] Add Poiseuille pipe benchmarks
- [ ] Expand Couette-Poiseuille cases

4) CI/CD
- [ ] Build/test/fmt/clippy automation
- [ ] Artifact caching
- [ ] Coverage reporting

### Success Metrics
- Green build/tests across workspace ✅
- Added validations pass with <1% relative error for analytical baselines
- Docs coverage for public items >= 80%
- All modules < 500 LOC (90% achieved)

### Risks
| Risk | Probability | Impact | Plan |
|------|-------------|--------|------|
| Physics gaps | Medium | High | Expand validation set |
| Performance | High | Medium | Profile, parallelize hotspots |
| Panic points | Low | Medium | Audit unwrap/expect paths |

### Decisions
- Correctness and clarity prioritized over performance
- No adjective-based identifiers; domain terms only ✅
- Prefer traits/enums over factories for composition ✅
- Modular architecture with domain-based separation ✅

### TRL
- 4 (component validation in lab)

### Acceptance Criteria (Research)
- [x] Compiles and tests pass
- [x] Baseline analytical validations
- [x] Examples/benches compile
- [x] Clean modular architecture

### Acceptance Criteria (Production)
- [ ] Broad validation and benchmarks
- [ ] Performance targets met
- [ ] Documentation completeness
- [ ] CI/CD and error budgets