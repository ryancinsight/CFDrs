# CFD Simulation Suite

A modular Computational Fluid Dynamics (CFD) simulation framework in Rust, emphasizing clean architecture, performance optimization, and extensibility.

## Architecture

The suite is organized into 8 specialized crates:

- **cfd-core**: Core abstractions, fluid properties, boundary conditions, compute dispatch
- **cfd-math**: Numerical methods, linear solvers, SIMD operations (restructured)
- **cfd-mesh**: Mesh generation, topology, quality metrics
- **cfd-io**: File I/O (VTK, CSV, binary), checkpointing, optional HDF5
- **cfd-1d**: 1D pipe networks, microfluidics simulation
- **cfd-2d**: 2D solvers, SIMPLE/PISO algorithms, LBM foundations
- **cfd-3d**: 3D FEM, spectral methods, multiphase foundations
- **cfd-validation**: Convergence studies, error metrics, benchmarks

## Current State: ALPHA (Honest Assessment)

### ✅ Successfully Implemented
- **SIMD Architecture**: Properly restructured with architecture-conditional dispatch (AVX2/SSE/NEON/SWAR)
- **GPU Infrastructure**: WGPU integration with proper shaders (advection, diffusion, pressure, velocity)
- **Modular Design**: Clean separation of concerns, modules under 500 lines
- **Build System**: Fixed HDF5 optional dependency, clean builds without warnings
- **Core Solvers**: SIMPLE/PISO pressure-velocity coupling functional
- **Discretization**: Central, Upwind, Power Law, QUICK schemes
- **Linear Solvers**: CG, BiCGSTAB implementations

### ⚠️ Partially Implemented
- **GPU Kernels**: WGSL shaders present but dispatch not fully connected
- **Turbulence Models**: k-ε, k-ω SST structures in place, validation needed
- **Multiphase**: VOF/Level Set foundations, missing reconstruction
- **LBM**: Collision operators defined, streaming incomplete
- **Examples**: Several examples have API mismatches needing fixes

### ❌ Known Issues
- **Documentation**: Missing docs for many public APIs (670 clippy warnings remaining)
- **Validation**: Literature benchmarks not fully validated
- **Performance**: Some unnecessary clones remain in critical paths
- **GPU Integration**: WGSL shaders present but dispatch not fully connected

## Performance Status

### SIMD Optimization
- **x86_64**: AVX2 (256-bit) and SSE4.2 (128-bit) paths implemented
- **ARM**: NEON (128-bit) support for AArch64
- **Fallback**: SWAR (Software SIMD) for unsupported architectures
- **Zero-copy**: Work in progress, some unnecessary clones remain

### GPU Acceleration
- **Backend**: WGPU for cross-platform support (Vulkan/Metal/DX12)
- **Kernels**: 4 compute shaders implemented (advection, diffusion, pressure, velocity)
- **Status**: Infrastructure ready, dispatch integration incomplete

## Building

### Requirements
- Rust 1.82+ (2021 edition currently, not 2025)
- Optional: HDF5 libraries for HDF5 support (properly feature-gated)

### Build Commands
```bash
# Basic build (no GPU, no HDF5)
cargo build --release --no-default-features

# With GPU support (default)
cargo build --release

# With all features (requires HDF5 system libraries)
cargo build --release --all-features
```

## Design Principles Applied

### Successfully Enforced
- **SSOT**: Single implementation per operation
- **Modular Structure**: simd/operations.rs split into ops/{mod,traits,x86,arm,fallback}.rs
- **Clean Naming**: No adjective-based names (Enhanced*, Optimized* removed)
- **Feature Gates**: Proper conditional compilation for optional dependencies

### In Progress
- **Zero-Copy**: Still have clones in critical paths (e.g., phi_new in solvers)
- **SLAP**: Some functions mix abstraction levels
- **Complete Testing**: Many tests disabled or incomplete

## Quick Start (Working Example)

```rust
use cfd_core::prelude::*;
use cfd_core::error::Result;

fn main() -> Result<()> {
    // Create fluid properties
    let fluid = ConstantPropertyFluid::<f64>::water();
    
    // Set up 2D grid
    let grid = StructuredGrid2D::<f64>::new(
        100, 100,  // nx, ny
        0.0, 1.0,  // x bounds  
        0.0, 1.0   // y bounds
    );
    
    // Note: Full solver integration still needs work
    // See examples directory for current capabilities
    
    Ok(())
}
```

## Development Roadmap

### Immediate Priorities
1. Fix Field2D API (add missing `set` method)
2. Complete GPU kernel dispatch integration
3. Fix example compilation errors
4. Implement missing SWAR operations

### Medium Term
1. Complete turbulence model validation
2. Finish LBM streaming implementation
3. Add unstructured mesh support
4. Comprehensive benchmarking suite

### Long Term
1. Full multiphase flow capability
2. Adaptive mesh refinement
3. Parallel domain decomposition
4. Production-ready API stability

## Contributing

Contributions welcome! Please ensure:
- Code follows Rust idioms and safety guidelines
- Modules stay under 500 lines (enforced)
- No redundant implementations (SSOT principle)
- Tests pass before submitting PRs
- Document public APIs

## Testing

```bash
# Run all tests (currently some fail)
cargo test --workspace --no-default-features

# Run benchmarks
cargo bench --no-default-features
```

## License

MIT OR Apache-2.0

## References

- Patankar, S.V. (1980). Numerical Heat Transfer and Fluid Flow
- Versteeg, H.K. & Malalasekera, W. (2007). An Introduction to Computational Fluid Dynamics
- Leonard, B.P. (1979). A stable and accurate convective modelling procedure

## Acknowledgments

This codebase underwent significant refactoring to enforce clean architecture principles. While functional, it requires additional work to reach production readiness. The honest assessment above reflects the current state as of this refactoring cycle.