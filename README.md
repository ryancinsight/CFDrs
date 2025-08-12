# CFD Simulation Suite

A Rust-based computational fluid dynamics (CFD) simulation framework for 1D, 2D, and 3D problems. This is an educational and research project that implements various CFD algorithms and numerical methods.

## Project Status

⚠️ **This project is under active development and is not production-ready.**

### What's Implemented

- **1D Solvers**: Basic pipe flow and network analysis
- **2D Solvers**: 
  - SIMPLE algorithm for incompressible flows (with known issues)
  - Finite Difference Method (FDM)
  - Finite Volume Method (FVM)
  - Lattice Boltzmann Method (LBM)
  - PISO algorithm
  - Vorticity-Stream function solver
- **3D Solvers**: Basic implementations of FEM, spectral methods, IBM, Level Set, and VOF
- **Linear Solvers**: Conjugate Gradient, GMRES (needs debugging), BiCGSTAB
- **Mesh Support**: Basic structured grids and some unstructured mesh capabilities

### Known Issues and Limitations

#### Critical Bugs
- ❌ **SIMPLE Solver**: Contains hardcoded grid spacing that was recently fixed but needs more testing
- ❌ **GMRES Solver**: Unstable implementation with very loose test tolerances (0.2)
- ❌ **Validation Suite**: Uses overly loose tolerances (5-20%) that don't properly validate correctness

#### Design Issues
- **Code Duplication**: Multiple inconsistent implementations of the same functionality
- **Incomplete Boundary Conditions**: Many boundary condition types are defined but not implemented
- **Missing Preconditioners**: Only identity preconditioner is implemented, limiting solver performance
- **Misleading Performance Claims**: Code contains "zero-copy" comments despite frequent use of `.clone()`

#### Numerical Methods
- **QUICK Scheme**: Recently corrected from centered to upwinded implementation
- **Momentum Solver**: Uses explicit iteration which limits stability and time step size
- **Convergence Checking**: Now includes momentum residuals but needs more testing

### What's NOT Production-Ready

- No comprehensive validation against experimental data
- Limited test coverage with loose tolerances
- Performance has not been optimized
- Many algorithms are basic implementations without advanced features
- Documentation is incomplete and sometimes misleading

## Building and Running

### Prerequisites

- Rust nightly toolchain (required for some dependencies)
- Basic development tools (git, cargo)

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/cfd-suite.git
cd cfd-suite

# Build the project
cargo build --release

# Run tests (note: some tests may fail or use loose tolerances)
cargo test

# Run examples
cargo run --example simple_pipe_flow
```

## Architecture

The project is organized into several crates:

- `cfd-core`: Core abstractions and interfaces
- `cfd-math`: Mathematical utilities and linear solvers
- `cfd-mesh`: Mesh handling and quality metrics
- `cfd-1d`: 1D solvers for pipe networks
- `cfd-2d`: 2D grid-based solvers
- `cfd-3d`: 3D mesh-based solvers
- `cfd-io`: Input/output operations
- `cfd-validation`: Validation benchmarks (needs improvement)

## Contributing

This is an open-source project and contributions are welcome. However, please be aware that:

1. The codebase has significant technical debt
2. Many features are partially implemented
3. Test coverage needs improvement
4. Performance optimization is needed

Priority areas for contribution:
- Fixing the GMRES implementation
- Implementing proper preconditioners (ILU, SOR, etc.)
- Improving validation suite with tighter tolerances
- Refactoring to eliminate code duplication
- Implementing implicit momentum solvers
- Adding comprehensive boundary condition support

## References

The implementations are based on various CFD textbooks and papers:

- Patankar (1980) - Numerical Heat Transfer and Fluid Flow
- Versteeg & Malalasekera (2007) - An Introduction to Computational Fluid Dynamics
- Anderson (1995) - Computational Fluid Dynamics: The Basics with Applications
- Various research papers for specific algorithms

## License

MIT License - See LICENSE file for details

## Disclaimer

This software is provided as-is for educational and research purposes. It is not suitable for production use or critical applications. The authors make no guarantees about the accuracy, stability, or performance of the implementations.