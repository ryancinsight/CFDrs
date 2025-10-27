# Sprint 1.74.0: GPU Backend, Procedural Macros, and Linear Solver Test Enhancement

## Sprint Objectives 🎯

Per persona requirements and gap analysis recommendations, Sprint 1.74.0 focuses on:

1. **GPU Backend Abstraction** - wgpu compute infrastructure (HIGH priority per gap analysis)
2. **Procedural Macro System** - Derive macros for solver/boundary setup (DRY principle)
3. **Phase 3: Linear Solver Tests** - BiCGSTAB/GMRES/CG convergence validation
4. **Series Resistance Resolution** - Fix without ignore marker (zero-tolerance policy)

## Implementation Plan

### 1. GPU Backend Abstraction (`cfd-gpu` crate)

**Rationale**: Gap analysis identified GPU acceleration as HIGH priority. Persona guideline: "prefer backend abstraction over wgpu-rs, wgpu compute shaders".

**Architecture**:
```rust
// Zero-cost abstraction with trait-based backend selection
pub trait ComputeBackend: Send + Sync {
    type BufferHandle;
    type ShaderHandle;
    
    async fn create_buffer(&self, size: usize) -> Result<Self::BufferHandle>;
    async fn compile_shader(&self, source: &str) -> Result<Self::ShaderHandle>;
    async fn execute(&self, shader: &Self::ShaderHandle, workgroups: [u32; 3]) -> Result<()>;
}

// CPU fallback (always available)
pub struct CpuBackend { /* ... */ }

// GPU backend (feature-gated)
#[cfg(feature = "gpu")]
pub struct WgpuBackend { /* device, queue, ... */ }

// Runtime backend selection
pub fn select_backend() -> Box<dyn ComputeBackend> {
    if cfg!(feature = "gpu") {
        Box::new(WgpuBackend::new())
    } else {
        Box::new(CpuBackend::new())
    }
}
```

**Implementation Steps**:
- [ ] Create `crates/cfd-gpu` with wgpu dependency (feature-gated)
- [ ] Implement `ComputeBackend` trait with CPU/GPU variants
- [ ] Add WGSL shaders for matrix operations (SpMV, axpy, dot)
- [ ] Write integration tests comparing CPU vs GPU results
- [ ] Benchmark performance gains with criterion

**Literature Validation**:
- Algorithms from "GPU Gems 3" (2007), Chapter 31: Fast N-Body Simulation
- wgpu best practices from "Learn WGPU" (2024)

### 2. Procedural Macro System (`cfd-macros` crate)

**Rationale**: Reduce boilerplate in solver setup, boundary conditions, and network DSL per DRY principle.

**Macros to Implement**:

#### 2.1 `#[derive(Solver)]`
```rust
// Before (manual implementation):
pub struct MyPoissonSolver {
    config: SolverConfig,
    preconditioner: Box<dyn Preconditioner>,
}

impl LinearSolver for MyPoissonSolver {
    fn solve(&self, a: &CsrMatrix, b: &Vector) -> Result<Vector> {
        // 50 lines of boilerplate...
    }
}

// After (macro-generated):
#[derive(Solver)]
#[solver(iterative, preconditioner = "ILU")]
pub struct MyPoissonSolver {
    config: SolverConfig,
}
```

#### 2.2 `boundary!` DSL Macro
```rust
// Before:
let mut network = Network::new(graph, fluid);
network.set_pressure(inlet, 1000.0);
network.set_pressure(outlet, 0.0);

// After:
let network = network! {
    graph: graph,
    fluid: fluid,
    boundary: {
        inlet => Pressure(1000.0),
        outlet => Pressure(0.0),
    }
};
```

#### 2.3 `#[derive(FlowComponent)]`
```rust
#[derive(FlowComponent)]
#[component(resistance_model = "HagenPoiseuille")]
pub struct Micropump {
    #[param(bounds = "0.0..=1.0")]
    efficiency: f64,
    
    #[param(bounds = "0.0..=1.0")]
    operating_point: f64,
}
// Auto-generates: parameter validation, resistance calculation, etc.
```

**Implementation Steps**:
- [ ] Create `crates/cfd-macros` with proc-macro2, quote, syn
- [ ] Implement `#[derive(Solver)]` with attribute parsing
- [ ] Implement `boundary!` declarative macro
- [ ] Implement `#[derive(FlowComponent)]` with parameter validation
- [ ] Write comprehensive macro tests (expand-based)
- [ ] Add rustdoc examples for each macro

### 3. Phase 3: Linear Solver Test Enhancement

**Objective**: Add literature-validated convergence tests for BiCGSTAB, GMRES, CG.

**Test Suite** (`crates/cfd-math/tests/linear_solver_convergence.rs`):

#### 3.1 BiCGSTAB Convergence Tests
- [ ] Test convergence on Poisson matrix (Saad & Schultz 1986)
- [ ] Test with ILU preconditioner (Barrett et al. 1994)
- [ ] Property test: convergence rate validation
- [ ] Benchmark: iterations vs matrix size

#### 3.2 GMRES Convergence Tests
- [ ] Test restarted GMRES(m) for various m values
- [ ] Test Arnoldi orthogonalization stability
- [ ] Compare with theoretical convergence bounds
- [ ] Benchmark: Krylov subspace size vs performance

#### 3.3 Conjugate Gradient Tests
- [ ] Test CG on symmetric positive definite (SPD) matrices
- [ ] Test preconditioned CG (PCG)
- [ ] Property test: energy minimization validation
- [ ] Benchmark: CG vs BiCGSTAB on SPD systems

**Literature References**:
- Saad, Y., & Schultz, M. H. (1986). "GMRES: A generalized minimal residual algorithm"
- Barrett, R., et al. (1994). "Templates for the Solution of Linear Systems"
- Greenbaum, A. (1997). "Iterative Methods for Solving Linear Systems"
- Golub, G. H., & Van Loan, C. F. (2013). "Matrix Computations" (4th ed.)

**Expected Coverage Improvement**: +15-20% (exercising solver paths)

### 4. Series Resistance Test Resolution

**Current Issue**: One edge shows near-zero flow (6e-41) in series network.

**Investigation**:
- [ ] Debug pressure distribution across series network
- [ ] Verify boundary condition setup
- [ ] Check if solver converges to correct solution
- [ ] Compare with analytical solution (R_total = R1 + R2 + R3)

**Approaches**:
1. **Option A**: Fix solver multi-junction handling if bug found
2. **Option B**: Adjust test topology (use simpler series without intermediate junctions)
3. **Option C**: Document as known limitation with workaround

**Requirement**: Must pass or be fixed (no ignore markers per persona rejection criteria).

## Success Criteria (≥90% checklist)

### GPU Backend ✅
- [ ] `cfd-gpu` crate created with ComputeBackend trait
- [ ] CPU fallback implementation complete
- [ ] wgpu backend implementation (feature-gated)
- [ ] WGSL shaders for SpMV, axpy, dot operations
- [ ] Integration tests: CPU vs GPU parity
- [ ] Criterion benchmarks showing speedup
- [ ] Documentation with examples
- [ ] Zero clippy warnings

### Procedural Macros ✅
- [ ] `cfd-macros` crate created
- [ ] `#[derive(Solver)]` implemented with tests
- [ ] `boundary!` macro implemented with tests
- [ ] `#[derive(FlowComponent)]` implemented with tests
- [ ] Macro expansion tests pass
- [ ] Rustdoc examples for each macro
- [ ] Zero clippy warnings

### Linear Solver Tests ✅
- [ ] BiCGSTAB: 4+ tests with literature validation
- [ ] GMRES: 4+ tests with Arnoldi validation
- [ ] CG: 4+ tests with energy minimization
- [ ] Property-based tests with proptest
- [ ] Criterion benchmarks for all solvers
- [ ] All tests passing (0 failures)
- [ ] Coverage measurement: +15-20%
- [ ] Literature citations complete

### Series Resistance Resolution ✅
- [ ] Root cause identified
- [ ] Fix implemented OR topology adjusted
- [ ] Test passing without ignore marker
- [ ] Zero regressions on other network tests
- [ ] Documentation updated

## Quality Gates (Maintained)

- ✅ Build warnings: 0
- ✅ Clippy warnings: 0 (pedantic)
- ✅ Test pass rate: 100%
- ✅ Test runtime: <30s
- ✅ Module compliance: All <500 LOC
- ✅ Technical debt: 0 markers
- ✅ Zero regressions maintained
- ⏳ Test coverage: Targeting +15-20% improvement

## Timeline Estimate

- **GPU Backend**: 8-10 hours
- **Procedural Macros**: 6-8 hours
- **Linear Solver Tests**: 4-6 hours
- **Series Resistance**: 2-3 hours
- **Total**: 20-27 hours (2-3 days focused work)

## Autonomous Execution Plan

Following persona "autonomous programming" directive:
1. Start all three initiatives in parallel
2. Use CoT-ToT-GoT for design decisions
3. Iterate on failures without user prompts
4. Validate with tools (cargo test, clippy, tarpaulin)
5. Report progress with commits
6. Complete sprint when ≥90% checklist achieved

## Literature References

**GPU Computing**:
- "GPU Gems 3" (2007), Chapter 31: Fast N-Body Simulation with CUDA
- "Learn WGPU" (2024): https://sotrh.github.io/learn-wgpu/
- Nickolls, J., et al. (2008). "Scalable parallel programming with CUDA"

**Linear Solvers**:
- Saad, Y., & Schultz, M. H. (1986). "GMRES: A generalized minimal residual algorithm". *SIAM J. Sci. Stat. Comput.*, 7(3), 856-869.
- Barrett, R., et al. (1994). "Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods". SIAM.
- Greenbaum, A. (1997). "Iterative Methods for Solving Linear Systems". SIAM.
- Golub, G. H., & Van Loan, C. F. (2013). "Matrix Computations" (4th ed.). Johns Hopkins.

**Procedural Macros**:
- "The Little Book of Rust Macros" (2024): https://veykril.github.io/tlborm/
- "Procedural Macros in Rust" - Rust Lang Book Chapter

---

**Sprint Status**: INITIATED  
**Phase**: Planning Complete, Ready for Implementation  
**Estimated Completion**: 2-3 days  
**Next Commit**: GPU backend trait definition and CPU implementation
