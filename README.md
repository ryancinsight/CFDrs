# CFD Simulation Suite

A modular Computational Fluid Dynamics (CFD) simulation framework in Rust, emphasizing clean architecture, performance optimization, and extensibility.

## Architecture

The suite is organized into 8 specialized crates:

- **cfd-core**: Core abstractions, fluid properties, boundary conditions
- **cfd-math**: Numerical methods, linear solvers, SIMD operations
- **cfd-mesh**: Mesh generation, topology, quality metrics
- **cfd-io**: File I/O, VTK support, checkpointing
- **cfd-1d**: 1D pipe networks, microfluidics
- **cfd-2d**: 2D solvers, SIMPLE/PISO algorithms
- **cfd-3d**: 3D FEM, spectral methods, multiphase foundations
- **cfd-validation**: Convergence studies, error metrics, benchmarks

## Features

### Implemented
- ✅ SIMPLE/PISO pressure-velocity coupling
- ✅ Discretization schemes (Central, Upwind, Power Law, QUICK)
- ✅ Linear solvers (CG, BiCGSTAB)
- ✅ Time integration (Euler, RK4)
- ✅ SIMD vectorization (AVX2/SSE/NEON with SWAR fallback)
- ✅ Structured mesh support
- ✅ VTK file I/O

### In Progress
- ⚠️ GPU acceleration (wgpu infrastructure present)
- ⚠️ Turbulence models (k-ε, k-ω SST partial)
- ⚠️ Multiphase methods (VOF, Level Set foundations)
- ⚠️ Unstructured mesh support

## Quick Start

```rust
use cfd_suite::prelude::*;
use cfd_suite::core::Result;

fn main() -> Result<()> {
    // Create fluid properties
    let fluid = Fluid::<f64>::water()?;
    
    // Set up 2D grid
    let grid = StructuredGrid2D::<f64>::new(
        100, 100,  // nx, ny
        0.0, 1.0,  // x bounds
        0.0, 1.0   // y bounds
    )?;
    
    // Configure solver
    let config = SolverConfig::default()
        .with_time_step(0.001)
        .with_max_iterations(1000);
    
    // Run simulation
    // ... solver implementation
    
    Ok(())
}
```

## Building

### Requirements
- Rust 1.89+ (2025 edition)
- Optional: HDF5 libraries for HDF5 support

### Basic Build
```bash
cargo build --release
```

### With GPU Support
```bash
cargo build --release --features gpu
```

### With All Features
```bash
# Note: Requires HDF5 system libraries
cargo build --release --all-features
```

## Design Principles

The codebase follows these principles:
- **SOLID**: Single responsibility, clean interfaces
- **CUPID**: Composable, Unix philosophy, predictable, idiomatic, domain-based
- **SSOT**: Single Source of Truth for all implementations
- **Zero-Copy**: Minimizing allocations (work in progress)

## Current Status: ALPHA

This is an alpha release undergoing active development:
- Core architecture refactored for modularity
- SIMD implementation unified
- Some features incomplete or untested
- API subject to change

## Contributing

Contributions welcome! Please ensure:
- Code follows Rust idioms
- Modules stay under 500 lines
- No redundant implementations
- Comprehensive tests for new features

## License

MIT OR Apache-2.0

## References

- Patankar, S.V. (1980). Numerical Heat Transfer and Fluid Flow
- Versteeg, H.K. & Malalasekera, W. (2007). An Introduction to Computational Fluid Dynamics
- Leonard, B.P. (1979). A stable and accurate convective modelling procedure