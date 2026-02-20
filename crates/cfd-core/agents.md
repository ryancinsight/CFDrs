# cfd-core — Agent Reference

> **Role**: Foundation crate. Every other crate in the workspace depends on `cfd-core`.  
> **No external solver dependencies** — this crate must remain at the bottom of the DAG.

---

## Purpose

`cfd-core` provides:
- **Physics models**: fluid properties, boundary conditions, physical values (newtypes), cavitation, haemolysis, material interfaces
- **Compute infrastructure**: CPU/GPU/MPI dispatch, SIMD kernels, solver abstraction, time integration
- **Geometry domain**: `Domain`, structured grid shapes, mesh connectivity for FEM/FVM consumers
- **Plugin system**: `SimulationPlugin` registry, aggregates (problem, state, parameters, metadata)
- **Prelude**: curated imports (`Fluid`, `Pressure`, `Velocity`, `BoundaryCondition`, …) for downstream crates

---

## Module Structure

```
src/
  lib.rs                   prelude + top-level pub mods
  abstractions/
    problem.rs             Problem statement abstraction
    state.rs               SimulationState
  compute/
    backend.rs             UnifiedCompute — CPU/GPU/Hybrid selection
    cpu.rs                 CPU compute paths
    dispatch.rs            Backend dispatcher
    gpu/
      mod.rs
      buffer.rs            WGPU buffer management
      constants.rs         WGSL constant definitions
      field_ops.rs         GPU field operations
      kernels/             Advection, diffusion, laplacian, pressure, velocity, turbulence, arithmetic
      pipeline.rs          Compute pipeline setup
      poisson_solver.rs    GPU Poisson solve
      shaders.rs           WGSL shader source
      turbulence_compute.rs GPU LES/DES kernels
      validation_tests.rs
    mpi/
      communicator.rs      MPI communicator wrapper
      decomposition.rs     Cartesian domain decomposition
      distributed_grid.rs  Distributed structured grid
      distributed_solvers.rs Parallel CG/GMRES
      ghost_cells.rs       Halo exchange
      performance_validation.rs Scaling benchmarks
    simd.rs                Architecture-conditional SIMD dispatch
    simd/
      x86.rs               AVX2/SSE4.1 kernels
      aarch64.rs           NEON kernels
    solver/
      config.rs            SolverConfig — tolerance, max_iter, backend
      convergence.rs       ConvergenceCriteria, ConvergenceMonitor
      direct.rs            Direct solver wrappers
      iterative.rs         IterativeSolverResult
      mod.rs
      monitor.rs           SolverMonitor — residual history
      traits.rs            Solver, IterativeLinearSolver
    tests.rs
    time/
      controllers.rs       TimeStepController (CFL-based)
      integrators.rs       TimeIntegrator trait
      mod.rs
    traits.rs              ComputeBackend trait
  error.rs                 Error enum + Result alias
  error/enhanced.rs        Structured error context
  geometry/
    mesh/                  Connectivity, operations, quality, refinement, service, statistics
    mod.rs                 Domain, DomainBounds
    shapes.rs              Geometric primitives for domain construction
  management/
    aggregates/            Problem, simulation, state, metadata, parameters, configuration aggregates
    broadcast.rs           Event broadcasting
    conversion.rs          Type conversion utilities (pub-re-exported)
    factory.rs             SimulationFactory
    mod.rs
    plugin/
      dependency.rs        Plugin dependency resolution
      health.rs            Plugin health checking
      mod.rs
      registry.rs          SimulationPluginRegistry
      storage.rs           Plugin storage
      traits.rs            Plugin, SimulationPlugin traits
  physics/
    boundary/              BoundaryCondition, WallType, manager, inlet/outlet/wall/ghost cells
    cavitation/            RayleighPlesset, cavitation number, venturi models, damage
    constants/             Physical and mathematical constants
    fluid/                 Fluid, ConstantPropertyFluid, blood, non-Newtonian, temperature dependence
    fluid_dynamics/        FlowField, FlowOperations, RANSModel, TurbulenceModel, RhieChow, FlowClassifier
    hemolysis.rs           Giersiepen–Wurzinger haemolysis correlation, BloodTrauma, PlateletActivation
    material/              SolidProperties, MaterialDatabase, interface coupling
    values/                Pressure, Velocity, Temperature, ReynoldsNumber (newtypes)
```

---

## Key Public API (`prelude`)

```rust
use cfd_core::prelude::*;

// Physical values (newtypes with dimensional safety)
let p = Pressure::new(101_325.0);      // Pa
let v = Velocity::new(0.5);            // m/s
let t = Temperature::new(310.15);      // K

// Fluid model
let fluid = ConstantPropertyFluid::<f64>::water();
let blood  = ConstantPropertyFluid::<f64>::blood_casson();

// Boundary conditions
let bc = BoundaryCondition::inlet_velocity(v);
let wall = WallType::NoSlip;

// Solver configuration
let conf = SolverConfig { max_iter: 1000, tolerance: 1e-8, ..Default::default() };

// Plugin system
struct MyPlugin;
impl SimulationPlugin for MyPlugin { /* … */ }
```

---

## Physics Subsystems

### Fluid Properties

| Type | Model |
|------|-------|
| `ConstantPropertyFluid` | Newtonian; `.water()`, `.blood_casson()`, `.air()` constructors |
| `NonNewtonianFluid` | Power-law, Bingham, Carreau-Yasuda |
| `TemperatureDependentFluid` | Viscosity as polynomial in T |
| `FluidDatabase` | Lookup by name |

### Cavitation (`physics/cavitation/`)

**Rayleigh-Plesset equation** (*Theorem*):
```
R·R̈ + (3/2)·Ṙ² = (1/ρ)·[pᵢ(t) - p∞(t) - 4μṘ/R - 2σ/R]
```
Stages: `Inception → Growth → Collapse`. `CavitationNumber = (p∞ - pᵥ) / (½ρV²)`.  
`VenturiCavitation` models throat pressure drop and bubble collapse radius.

### Haemolysis (`physics/hemolysis.rs`)

**Giersiepen–Wurzinger correlation** (*Theorem*):
```
HI(τ, t) = 3.62 × 10⁻⁷ · τ^2.416 · t^0.785
```
`HemolysisCalculator` accumulates index along pathlines.  
`BloodTraumaSeverity`: Negligible / Moderate / Severe / Critical.  
FDA threshold: `≤ 150 Pa` sustained wall shear for blood-contacting devices.

### Boundary Conditions (`physics/boundary/`)

| Type | Description |
|------|-------------|
| `Inlet` | Fixed velocity or mass-flow |
| `Outlet` | Zero-gradient or fixed pressure |
| `Wall` | No-slip or partial slip (`WallType`) |
| `Symmetry` | Zero-normal-gradient |
| `Periodic` | Matched-pair coupling |
| `TimeDependentBC` | Womersley pulsatile inlet |
| `GhostCellBC` | MPI ghost-cell population |

---

## Compute Subsystems

### GPU (`compute/gpu/`)

- Backend: `wgpu` (Vulkan / Metal / DX12 via WebGPU)
- Shader language: WGSL (`shaders.rs`)
- Kernels: advection, diffusion, pressure-Poisson, velocity update, Laplacian, turbulence (Smagorinsky LES)
- Pipeline: `ComputePipeline` manages bind groups and workgroup dispatch
- CPU fallback: `UnifiedCompute::select_backend()` falls back when GPU unavailable

### MPI (`compute/mpi/`)  *(feature = `mpi`)*

- `DomainDecomposition`: Cartesian 2D/3D partition with load balancing
- `GhostCellExchange`: halo communication for all field types
- `DistributedLinearSolver`: parallel GMRES / BiCGSTAB
- `AdaptiveMeshRefinement` with repartitioning
- `PerformanceValidation`: strong/weak scaling benchmarks

### SIMD (`compute/simd/`)

- `arch_detect` → selects AVX2 / SSE4.1 / NEON / SWAR at runtime
- **Note**: SpMV SIMD performance benchmark and decision are owned by `cfd-math`. See `crates/cfd-math/agents.md § SIMD Status`.

---

## Feature Flags

The workspace `agents.md` is the single source of truth for the feature flag registry.
Flags that affect this crate:

| Flag | Default | Effect on `cfd-core` |
|------|---------|----------------------|
| `gpu` | ✅ on | Enables `wgpu` device + all GPU compute kernels |
| `mpi` | off | Enables MPI communicator stack + `DistributedLinearSolver` |
| `simd` | off | Enables explicit SIMD dispatch paths (see `cfd-math` for performance rationale) |

> `csg`, `millifluidic`, and other flags do not affect `cfd-core`.

---

## Key Mathematical Theorems

| Module | Theorem |
|--------|---------|
| `physics/cavitation/rayleigh_plesset.rs` | Rayleigh-Plesset bubble dynamics |
| `physics/hemolysis.rs` | Giersiepen–Wurzinger haemolysis index |
| `compute/solver/convergence.rs` | CG converges in ≤n steps (SPD matrix) |
| `compute/gpu/kernels/diffusion.rs` | Discrete maximum principle for diffusion |

---

## Prohibited Patterns

- No circular imports back from `cfd-math`, `cfd-mesh`, `cfd-2d`, etc.
- `physics/` subsystems must not import from `compute/` (domain must stay pure)
- GPU kernel code must remain in `compute/gpu/kernels/` — no physics logic in shaders
- `management/plugin/` is the only layer allowed to import from both `physics/` and `compute/`
