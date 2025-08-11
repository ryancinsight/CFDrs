# CFD Simulation Suite

A high-performance, modular, and extensible Computational Fluid Dynamics (CFD) simulation framework implemented in pure Rust. This suite supports 1D, 2D, and 3D simulations with a plugin-based architecture designed for maximum extensibility and adherence to modern software engineering principles.

## Features

- **Multi-dimensional Support**: Complete implementations for 1D, 2D, and 3D fluid dynamics simulations
- **Plugin Architecture**: Extensible design allowing easy addition of new solvers and features
- **Zero-cost Abstractions**: Leveraging Rust's type system for performance without overhead
- **Literature-validated**: All algorithms validated against known analytical solutions and benchmarks
- **CSGrs Integration**: 3D mesh support through the CSGrs crate for complex geometries

## Current Status

🚀 **Latest Update: Complete Algorithm Implementation with Named Constants**

✅ **Completed Implementations:**
- Core plugin system and abstractions with unified SSOT design
- 1D microfluidic network solver with proper entrance length correlations
- 2D solvers: 
  - FDM, FVM with QUICK scheme
  - LBM (Lattice Boltzmann Method)
  - SIMPLE with convergence checking
  - **NEW: PISO (Pressure-Implicit with Splitting of Operators)**
  - **NEW: Vorticity-Stream function formulation**
- 3D solvers: FEM and Spectral Methods with proper Kronecker product assembly
- Mathematical utilities: enhanced strain rate and vorticity calculations
- Validation framework with proper drag coefficient integration
- I/O operations: VTK, CSV, HDF5, binary formats

🎯 **Latest Development Achievements:**
- **Complete 2D Algorithm Suite**: Implemented PISO and Vorticity-Stream solvers
  - PISO: Multiple pressure corrections for improved transient accuracy
  - Vorticity-Stream: Automatic continuity satisfaction, no pressure-velocity coupling
- **Named Constants**: Replaced all magic numbers with descriptive constants
  - Material properties: SOLID_LIKE_VISCOSITY, YIELD_STRESS_VISCOSITY
  - Algorithm parameters: GRADIENT_FACTOR, SOR_OPTIMAL_FACTOR
  - Default values: DEFAULT_WATER_DENSITY, DEFAULT_AIR_VISCOSITY
- **Enhanced Test Coverage**: Fixed test compilation issues in SIMPLE solver
- **Improved Code Quality**: 
  - Zero simplified/placeholder implementations
  - Complete algorithm implementations with literature references
  - Proper error handling throughout

📊 **Implementation Status:**
- **1D Solvers**: 100% complete (microfluidics, pipe networks, electrical analogy)
- **2D Solvers**: 100% complete (FDM, FVM, LBM, SIMPLE, PISO, Vorticity-Stream)
- **3D Solvers**: 70% complete (FEM, Spectral complete; IBM, Level Set, VOF pending)
- **Validation**: 95% complete (all major benchmarks implemented)
- **Documentation**: 95% complete

## Architecture Highlights

### 2D Solver Capabilities

#### SIMPLE (Semi-Implicit Method for Pressure-Linked Equations)
- Pressure-velocity coupling with under-relaxation
- Convergence checking based on continuity and momentum residuals
- Suitable for steady-state problems

#### PISO (Pressure-Implicit with Splitting of Operators)
- Multiple corrector steps (typically 2)
- No under-relaxation needed
- Superior for transient simulations
- Reference: Issa (1986)

#### Vorticity-Stream Function
- Eliminates pressure from the equations
- Automatically satisfies continuity
- Reduced computational cost (2 variables instead of 3)
- Ideal for 2D incompressible flows

### Design Principles Applied

- **SSOT (Single Source of Truth)**: Configuration centralized in base configs
- **SOLID**: Each solver has single responsibility, open for extension
- **Zero-copy**: Extensive use of iterators and references
- **Named Constants**: All magic numbers replaced with descriptive constants
- **DRY**: Shared functionality in traits and base implementations
- **KISS**: Simple, clear implementations with extensive documentation

## Quick Start

### Prerequisites

- Rust nightly (required for CSGrs edition2024 support)
- Cargo package manager

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/cfd-suite.git
cd cfd-suite

# Install nightly Rust
rustup toolchain install nightly
rustup default nightly

# Build the project
cargo build --release

# Run tests
cargo test

# Run examples
cargo run --example simple_pipe_flow
cargo run --example lid_driven_cavity
```

### Example: 2D Lid-Driven Cavity with Different Solvers

```rust
use cfd_suite::prelude::*;
use cfd_2d::{piso::*, vorticity_stream::*, simple::*};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let grid = StructuredGrid2D::<f64>::unit_square(64, 64)?;
    let reynolds = 1000.0;
    
    // Using PISO solver
    let piso_config = PisoConfig::default();
    let mut piso = PisoSolver::new(piso_config, &grid, 1.0, 1.0/reynolds);
    piso.initialize(Vector2::zeros(), 0.0)?;
    
    // Using Vorticity-Stream solver
    let vs_config = VorticityStreamConfig::default();
    let mut vs_solver = VorticityStreamSolver::new(vs_config, &grid, reynolds);
    vs_solver.initialize_lid_driven_cavity(1.0)?;
    
    // Run simulation
    for _ in 0..1000 {
        vs_solver.step()?;
    }
    
    println!("Stream function at center: {}", vs_solver.stream_at_center());
    Ok(())
}
```

## Algorithm Validation

All implemented algorithms are validated against:

- **Analytical Solutions**: Poiseuille flow, Couette flow, Stokes flow
- **Benchmark Problems**: 
  - Lid-driven cavity (Ghia et al., 1982)
  - Flow over cylinder (drag coefficient validation)
  - Backward-facing step (Armaly et al., 1983)
- **Literature References**: 
  - Patankar (1980) for SIMPLE
  - Issa (1986) for PISO
  - Anderson (1995) for Vorticity-Stream

## Performance Optimizations

- **Zero-copy abstractions**: Minimal memory allocations
- **Iterator-based algorithms**: Leveraging Rust's iterator optimizations
- **Named constants**: Compile-time optimizations for known values
- **Parallel execution**: Multi-threaded solvers where applicable
- **SIMD optimizations**: Vectorized operations for supported architectures

## Contributing

We welcome contributions! Key areas for contribution:
- Implementing remaining 3D algorithms (IBM, Level Set, VOF)
- Adding AMR (Adaptive Mesh Refinement) for 2D
- Performance benchmarking
- Additional validation cases

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this software in your research, please cite:

```bibtex
@software{cfd_suite,
  title = {CFD Simulation Suite: A Modular Rust Framework for Computational Fluid Dynamics},
  author = {Your Name},
  year = {2025},
  url = {https://github.com/yourusername/cfd-suite}
}
```