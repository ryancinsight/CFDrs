# Product Requirements Document

## CFD Suite v0.57.4

### Product Classification
Research software (not production)

### Current Capabilities

#### Functional
- Workspace compiles and runs
- Core algorithms (FVM/FDM/PISO/VOF/spectral) present
- Analytical validations: Couette, Poiseuille (plates), Taylor-Green
- Test suite covers numerical methods and integrations

#### Non-Functional
- Trait-based, domain-structured crates
- Result-based error handling
- Examples and benches compile

### Limitations
- Validation coverage limited to selected cases
- Performance and parallelism deferred
- Docs incomplete; warnings for missing docs on public items

### Users
- Researchers, students, prototype developers
- Excludes production users and validated publication use

### Development Roadmap (near term)
1) Structure
- [x] Split `cfd-core/time.rs` into `time/integrators.rs`, `time/controllers.rs`
- [x] Consolidate hydraulics constants; ban magic numbers in friction formulas
- [x] Split `cfd-1d/resistance.rs` into domain-based modules
- [x] Remove duplicate/misnamed directories (crases)
- [x] Fix all adjective-based naming violations

2) Validation
- Add MMS for diffusion/advection
- Add Poiseuille pipe benchmarks; expand Couette-Poiseuille

3) Quality
- Reduce warnings (unused imports/vars)
- Document public constants and fields
- Remove placeholder literature modules from public API until validated

4) CI
- Build/test/fmt/clippy; artifact caching

### Success Metrics
- Green build/tests across workspace
- Added validations pass with <1% relative error for analytical baselines
- Docs coverage for public items >= 80%

### Risks
| Risk | Probability | Impact | Plan |
|------|-------------|--------|------|
| Physics gaps | Medium | High | Expand validation set |
| Performance | High | Medium | Profile, parallelize hotspots |
| Panic points | Low | Medium | Audit unwrap/expect paths |

### Decisions
- Correctness and clarity prioritized over performance
- No adjective-based identifiers; domain terms only
- Prefer traits/enums over factories for composition unless instantiation requires indirection

### TRL
- 4 (component validation in lab)

### Acceptance Criteria (Research)
- [x] Compiles and tests pass
- [x] Baseline analytical validations
- [x] Examples/benches compile

### Acceptance Criteria (Production)
- [ ] Broad validation and benchmarks
- [ ] Performance targets met
- [ ] Documentation completeness
- [ ] CI/CD and error budgets