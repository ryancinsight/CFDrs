# CFD Suite Architectural Restructuring Plan

## Critical Violations Requiring Immediate Action

### 1. Monolithic Module Violations (GRASP/SOC)

#### cfd-core/src/fluid.rs (437 lines)
**Current State**: Conflates fluid properties, non-Newtonian models, temperature dependencies, and predefined fluids
**Required Action**: Split into domain-oriented modules:
```
fluid/
├── mod.rs           # Core FluidProperties trait
├── properties.rs    # Basic property calculations
├── newtonian.rs     # Newtonian fluid implementations
├── non_newtonian.rs # Power-law, Bingham, etc.
├── temperature.rs   # Temperature-dependent models
├── database.rs      # Predefined fluids (water, air, etc.)
└── validation.rs    # Property validation and bounds checking
```

#### cfd-2d/src/fields.rs (370 lines)
**Current State**: Mixes field storage, access patterns, and simulation state
**Required Action**: Decompose into:
```
fields/
├── mod.rs        # Field traits and interfaces
├── storage.rs    # Field2D storage implementation
├── access.rs     # Zero-copy access patterns
├── state.rs      # SimulationFields aggregate
└── operators.rs  # Field operations (add, multiply, etc.)
```

### 2. Zero-Copy Violations (17 remaining clones)

**Location Analysis**:
- cfd-2d/src/schemes/time_integration.rs: 4 clones
- cfd-1d/src/solver/mod.rs: 3 clones
- cfd-2d/src/solvers/accelerated.rs: 2 clones
- cfd-2d/src/pressure_velocity/rhie_chow.rs: 2 clones

**Resolution Strategy**:
- Replace Vec clones with slice references
- Use Cow<'_, [T]> for conditional ownership
- Implement view-based field access
- Utilize iterator chains instead of collecting

### 3. Stub Implementation Violations

**Critical Stubs Identified**:
- GPU compute kernels return Ok(()) without computation
- Validation module exists but performs no validation
- 3D module structure present but non-functional

**Required Actions**:
- DELETE cfd-3d module entirely until 2D is validated
- DELETE GPU modules until CPU implementation is proven
- IMPLEMENT actual validation against literature benchmarks

### 4. Naming Convention Violations

**Identified Issues**:
- "Simple" in SIMPLE solver (acceptable - algorithm name)
- "Advanced" in comments (must be removed)
- "Enhanced" in old commit messages (historical, ignore)

### 5. Literature Validation Gaps

**Unvalidated Implementations**:
- QUICK scheme lacks Leonard (1979) validation
- Rhie-Chow interpolation incomplete
- VOF method non-functional despite claims
- Level Set exists as empty structure

**Validation Requirements**:
- Each numerical method must reference specific equation numbers
- Test cases must reproduce published benchmarks
- Error bounds must match literature predictions

## Implementation Priority

### Phase 1: Structural Decomposition (Immediate)
1. Split fluid.rs into modular structure
2. Decompose fields.rs into access patterns
3. Delete non-functional 3D and GPU modules

### Phase 2: Zero-Copy Implementation (Week 1)
1. Replace time integration clones with views
2. Implement Cow-based field ownership
3. Use iterator chains for transformations

### Phase 3: Stub Elimination (Week 2)
1. Implement actual validation tests
2. Complete Rhie-Chow interpolation
3. Fix pressure-velocity coupling

### Phase 4: Literature Validation (Week 3)
1. Validate QUICK against Leonard (1979)
2. Verify k-ε against DNS data
3. Benchmark lid-driven cavity against Ghia et al. (1982)

## Success Metrics

- Zero modules > 300 lines
- Zero clone operations in hot paths
- Zero Ok(()) stubs
- 100% literature-validated algorithms
- Zero unused variables
- Zero missing documentation

## Architectural Principles Enforcement

### SSOT (Single Source of Truth)
- All constants in dedicated modules
- No duplicate implementations
- Single validation framework

### CUPID (Composable, Unix philosophy, Predictable, Idiomatic, Domain-based)
- Plugin-based solver architecture
- Trait-based extensibility
- Domain-oriented module structure

### Zero-Cost Abstractions
- Iterator-based processing
- Compile-time dispatch
- No runtime overhead from abstractions

## Validation Against Rust Best Practices

Per The Rust Book Chapter 7:
- Modules should be small and focused
- Use mod.rs for module interfaces
- Keep implementation details private
- Expose clean public APIs

## Next Immediate Action

Execute structural decomposition of fluid.rs NOW.