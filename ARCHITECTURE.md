# Architecture

## Workspace structure

CFDrs is a Cargo workspace with 10 crates:

| Crate | Role |
|---|---|
| `cfd-core` | Core abstractions, fluid properties, boundary conditions, solvers, **canonical error types** |
| `cfd-math` | Numerical methods, linear solvers, SIMD operations |
| `cfd-io` | File I/O (VTK, HDF5, CSV), parallel I/O, checkpointing |
| `cfd-1d` | 1D pipe networks, microfluidics simulation |
| `cfd-2d` | 2D solvers, SIMPLE/PISO algorithms, LBM foundations |
| `cfd-3d` | 3D FEM, spectral methods, multiphase foundations |
| `cfd-schematics` | Microfluidic schematic design, geometry, configuration, visualization |
| `cfd-optim` | Design optimization, genetic algorithms, blueprint candidates |
| `cfd-validation` | Convergence studies, error metrics, benchmarks |
| `cfd-python` | Python bindings |

## Error handling strategy

All domain-specific error types live in **`cfd-core::error`** as a single canonical location. Other crates either re-export from cfd-core or define thin local types with `From` conversions into the central `Error` enum.

### `cfd_core::error::Error`

The unified error enum with 28 variants covering the full domain:

```
InvalidInput, InvalidConfiguration, Numerical(NumericalErrorKind),
GpuCompute, Convergence(ConvergenceErrorKind), Plugin(PluginErrorKind),
Solver, ConversionError, Io, Serialization, UnsupportedOperation,
PhysicsViolation, DimensionMismatch, IndexOutOfBounds, Boundary(BoundaryErrorKind),
Geometry(GeometryErrorKind), Configuration(ConfigurationErrorKind),
Visualization(VisualizationErrorKind), Strategy(StrategyErrorKind),
Parameter(ParameterErrorKind), Validation(ValidationErrorKind),
Registry(RegistryErrorKind), Constraint(ConstraintErrorKind),
Dependency(DependencyErrorKind), Adaptation(AdaptationErrorKind),
ResistanceCalculation(ResistanceCalculationErrorKind),
NotImplemented, WithContext
```

### Kind enums (15 total)

Each domain area has a `*ErrorKind` enum with:

- Named variants with structured fields (not just `String`)
- `Display` impl for human-readable messages
- `std::error::Error` impl for interop
- Inherent convenience constructors (e.g., `GeometryErrorKind::invalid_box_dimensions(w, h)`)
- `From<Kind> for Error` impl for automatic `?` conversion

The 15 Kind enums:

| Kind | Domain |
|---|---|
| `NumericalErrorKind` | Solver numerics (division by zero, singular matrix, etc.) |
| `ConvergenceErrorKind` | Iterative solver convergence |
| `BoundaryErrorKind` | Boundary condition setup |
| `PluginErrorKind` | Plugin lifecycle |
| `GeometryErrorKind` | Geometry generation and validation |
| `ConfigurationErrorKind` | Config validation |
| `VisualizationErrorKind` | Rendering and output |
| `StrategyErrorKind` | Channel routing strategies |
| `ParameterErrorKind` | Parameter management (not found, read-only, etc.) |
| `ValidationErrorKind` | Validation rules |
| `RegistryErrorKind` | Parameter manager registry |
| `ConstraintErrorKind` | Range/set/custom constraints |
| `DependencyErrorKind` | Parameter dependency graphs |
| `AdaptationErrorKind` | Adaptive parameter behavior |
| `ResistanceCalculationErrorKind` | 1D resistance models |

### Cross-kind conversions

Two lossy `From` impls support state management workflows where constraint/validation errors propagate into parameter errors:

```rust
From<ConstraintErrorKind> for ParameterErrorKind  // → InvalidValue
From<ValidationErrorKind> for ParameterErrorKind   // → InvalidValue
```

### Crate-specific patterns

**cfd-schematics** — fully consolidated. All error types are re-exports or type aliases from cfd-core:

```rust
pub type GeometryError = GeometryErrorKind;
pub type ConfigurationError = ConfigurationErrorKind;
pub type VisualizationError = VisualizationErrorKind;
pub type StrategyError = StrategyErrorKind;
pub type SchemeError = Error;
pub type AdaptationError = AdaptationErrorKind;
// ... domain-specific Result types use Kind types directly
pub type GeometryResult<T> = Result<T, GeometryErrorKind>;
pub type ConfigurationResult<T> = Result<T, ConfigurationErrorKind>;
```

Callers use inherent constructors via type aliases (e.g., `ConfigurationError::invalid_arc_config(...)`), which resolve to `ConfigurationErrorKind::invalid_arc_config(...)`. The `From<Kind> for Error` impls handle `?` conversion in functions returning `Result<T, Error>`.

**cfd-1d** — `ResistanceCalculationError` consolidated into cfd-core. The old `solver/analysis/error.rs` was deleted; callers import via `use cfd_core::error::ResistanceCalculationErrorKind as ResistanceCalculationError`.

**cfd-optim** — keeps `OptimError` local (2 of 4 variants are optimizer-specific). A `From<OptimError> for Error` impl enables cross-crate `?` propagation:

```rust
impl From<OptimError> for cfd_core::error::Error {
    fn from(err: OptimError) -> Self {
        let msg = err.to_string();
        match err {
            OptimError::EmptyCandidates
            | OptimError::InsufficientCandidates { .. }
            | OptimError::InvalidParameter(_) => Self::InvalidInput(msg),
            OptimError::PhysicsError { .. } => Self::Solver(msg),
        }
    }
}
```

### Workspace-wide audit (confirmed)

All 10 crates use `cfd_core::error` as the single source of truth for errors. Zero local `Error` enums exist outside cfd-core and cfd-optim:

| Crate | Local error types | Status |
|---|---|---|
| `cfd-core` | Canonical `Error` + 15 Kind enums | **SSOT** |
| `cfd-math` | 0 | Uses `cfd_core::error` directly |
| `cfd-io` | 0 | Uses `cfd_core::error` directly |
| `cfd-1d` | 0 | `ResistanceCalculationError` consolidated into cfd-core |
| `cfd-2d` | 0 | Uses `cfd_core::error` directly |
| `cfd-3d` | 0 | Uses `cfd_core::error` directly |
| `cfd-schematics` | 0 | All former local types replaced with type aliases to Kind types |
| `cfd-optim` | `OptimError` (local) | `From<OptimError> for Error` bridge for cross-crate `?` propagation |
| `cfd-validation` | 0 | Uses `cfd_core::error` directly |
| `cfd-python` | 0 | Uses `cfd_core::error` directly |

### Guidelines for new error types

1. **Add the Kind enum to `cfd-core/src/error.rs`** with structured variants, `Display`, `std::error::Error`, inherent constructors, and `From<Kind> for Error`.
2. **Re-export via type alias** in the consuming crate if the old name is needed for backward compatibility.
3. **Keep local only if** the error type has variants with no meaningful core equivalent (like `OptimError::EmptyCandidates`). Add a `From` impl instead.
4. **Never define a local `Error` enum** in sub-crates that duplicates variants already in cfd-core.
