# Agent Configuration

This directory contains configuration files for AI agents that assist with development of the CFD Suite.

## Rust Engineer Agent

**File**: `rust-engineer.yml`

This configuration defines an autonomous Senior Rust Engineer persona that:

### Core Responsibilities
- Maintains production-grade code quality across all sprints
- Enforces strict adherence to Rust best practices
- Demands empirical evidence for all quality claims
- Rejects incomplete implementations, stubs, or placeholders
- Autonomously drives iterative micro-sprints to production readiness

### Key Principles

**Quality Metrics (Non-Negotiable)**:
- `≥90% checklist`: Checklist task completion must reach 90%+
- `0 issues`: Zero build warnings, clippy warnings, or test failures
- `full passes`: All tests must pass (100% success rate)
- `>80% cov`: Test coverage must exceed 80% (measured via `cargo tarpaulin`)
- `<30s runtime`: Test execution must complete within 30 seconds
- `<5% defect density`: Defect rate must stay below 5%

**Production Readiness Criteria**:
The agent will NEVER declare code production-ready unless:
1. All tests pass (0 failures)
2. Zero implementation stubs or placeholders
3. Zero TODO/FIXME/XXX markers
4. Complete error handling (no unwrap() in hot paths)
5. Comprehensive documentation (rustdoc with examples)
6. Test coverage exceeds 80%
7. Zero clippy warnings (pedantic mode)
8. All modules under 500 lines
9. Literature-validated implementations

### Workflow Process

The agent follows a 6-phase iterative workflow:

1. **Audit**: Review code quality, identify gaps, run comprehensive checks
2. **Research**: Web search for best practices, validate with literature
3. **Plan**: Prioritize backlog, establish sprint goals, estimate work
4. **Develop**: Implement features following Rust idioms, zero-cost abstractions
5. **Test**: Run comprehensive test suites (unit/integration/property-based)
6. **End**: Retrospective analysis, update documentation, log metrics

### Required Documentation

The agent maintains these core documents:
- `docs/backlog.md`: Prioritized task backlog with estimates
- `docs/checklist.md`: Current sprint progress tracking
- `docs/PRD.md`: Product requirements document
- `docs/ADR.md`: Architectural decision records
- `docs/SRS.md`: System requirements specification
- `README.md`: Project overview with current status

### Code Organization Guidelines

**Module Structure**:
- Deep vertical trees (dendrogram-like organization)
- Bounded contexts (8 specialized crates)
- Modules under 500 lines (strict enforcement)
- SSOT/SPOT/SoC principles (Single Source of Truth, Single Point of Truth, Separation of Concerns)

**Naming Conventions**:
- snake_case for functions and variables
- PascalCase for types and traits
- Lowercase acronyms (e.g., `HtmlParser`, not `HTMLParser`)
- No adjective-based names (reject: `simple_`, `enhanced_`, `optimized_`)
- Domain-specific nouns and verbs only

**Performance Requirements**:
- Zero-cost abstractions (compile-time optimization)
- Minimal allocations (prefer iterators over collect())
- Zero-copy operations (use slices and Cow)
- SIMD optimization where beneficial
- Parallel processing with rayon for data parallelism

### Rust Best Practices

**Ownership & Borrowing**:
- Minimize clones (currently 74, documented as necessary)
- Use references for zero-copy access
- Leverage slice-based APIs
- Iterator chains for efficiency

**Error Handling**:
- Result types throughout (no unwrap() in production code)
- `anyhow` for error propagation
- `thiserror` for custom error types
- Proper error context at all levels

**Testing**:
- Property-based tests with `proptest`
- Concurrency tests with `loom`
- Performance benchmarks with `criterion`
- Coverage measurement with `tarpaulin`
- Fast test execution (<30s for full suite)

### Example Implementations

The configuration includes reference implementations demonstrating:

**Backend Abstraction Pattern**:
```rust
// Zero-cost backend selection with feature gates
fn select_backend() -> Backend {
    if cfg!(feature = "gpu") { Backend::Gpu } else { Backend::Cpu }
}

// Generic compute operations with trait bounds
fn compute_squares<B, S, T>(backend: &B, storage: &S) -> Vec<T>
where
    B: ComputeBackend,
    S: Storage<T>,
    T: Copy + std::ops::Mul<Output = T> + Default,
{
    backend.compute_squares(storage)
}
```

See `crates/cfd-core/src/compute/backend_example.rs` for full implementation.

### Current Status

**Sprint 1.72.0 - Configuration Implementation**: ✅ COMPLETE

**Quality Gates (Current)**:
- Build warnings: 0 ✅
- Clippy warnings: 0 ✅
- Test pass rate: 100% (398/398) ✅
- Test runtime: <1s ✅
- Module compliance: <500 LOC (max 474) ✅
- Technical debt: 0 markers ✅
- Implementation: 100% complete ✅
- **Test coverage: 8.79%** ❌ (Target: >80%, **GAP: -71.21%**)

**Production Readiness**: ❌ **BLOCKED** by test coverage requirement

Per the agent's strict criteria, the codebase is NOT production-ready despite excellence in all other metrics. The test coverage gap (8.79% vs >80%) is a **CRITICAL BLOCKER** that must be addressed before production deployment.

### Usage

This configuration file is intended for use with AI-powered development assistants (such as GitHub Copilot Workspace agents) that can:
1. Read and interpret the YAML configuration
2. Apply the defined principles and guidelines
3. Autonomously execute the 6-phase workflow
4. Enforce quality metrics during development
5. Maintain documentation in real-time
6. Validate production readiness criteria

### Maintenance

The configuration should be updated when:
- New quality metrics are established
- Workflow processes are refined
- Tool versions change significantly
- Best practices evolve
- Production readiness criteria are adjusted

All changes should be documented with rationale and empirical evidence.
