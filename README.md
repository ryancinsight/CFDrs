# CFD Suite - Rust Implementation

**Version 0.57.3** - Research Software

## Status

- Builds and tests pass across workspace (benches/examples compile)
- Analytical validations included for Couette, Poiseuille (plates), Taylor-Green
- Domain-structured crates present; further module splits planned

## Verified Functionality
- ✅ Build succeeds (workspace)
- ✅ Tests pass (workspace)
- ✅ Examples compile and run
- ✅ Memory safe (Rust guarantees)
- ✅ Result-based error handling

## Technical Debt (tracked)
- Modules over 500 LOC to split by feature (e.g., `cfd-1d/resistance.rs`, `cfd-3d/level_set.rs`)
- Missing documentation warnings for constants and fields
- Validation scope to expand beyond initial cases
- Parallelization and performance not addressed

## Architecture
```
cfd-suite/
├── cfd-core/       # Core abstractions, plugin system, time (integrators/, controllers/)
├── cfd-math/       # Numerical methods, sparse CSR, solvers
├── cfd-mesh/       # Mesh, grid, quality, (CSG gated by feature)
├── cfd-1d/         # 1D networks and resistance models
├── cfd-2d/         # 2D fields, discretization, solvers
├── cfd-3d/         # 3D spectral, VOF, level set
├── cfd-io/         # I/O
└── cfd-validation/ # Analytical solutions, benchmarks, convergence tools
```

## Usage
```
cargo build --workspace
cargo test --workspace --all-targets
cargo run --example pipe_flow_1d --release
```

## Design Principles
- SSOT/SPOT: constants centralized; replace magic numbers with named constants
- SOLID/CUPID/GRASP/SLAP/DRY/CLEAN: traits for composition; avoid mixing concerns
- No adjective-bearing identifiers in APIs (ban subjective names); use domain terms

## Validation
- Analytical: Couette, Poiseuille (plates), Taylor-Green initial/decay
- Numerical: linear solver convergence tests and criteria
- Literature placeholder entries removed from public API until validated
- Add manufactured solutions and benchmark comparisons next

## Limits (non-exhaustive)
- Limited validation coverage beyond listed cases
- Performance and parallelism not targeted yet
- Missing docs for some public items

## TRL
- TRL 4 (component validation in lab environment)

## License
MIT OR Apache-2.0