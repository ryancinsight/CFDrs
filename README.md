# CFD Simulation Suite

A high-performance, modular, and extensible Computational Fluid Dynamics (CFD) simulation framework implemented in pure Rust. This suite supports 1D, 2D, and 3D simulations with a plugin-based architecture designed for maximum extensibility and adherence to modern software engineering principles.

## Features

- **Multi-dimensional Support**: Complete implementations for 1D, 2D, and 3D fluid dynamics simulations
- **Plugin Architecture**: Extensible design allowing easy addition of new solvers and features
- **Zero-cost Abstractions**: Leveraging Rust's type system for performance without overhead
- **Literature-validated**: All algorithms validated against known analytical solutions and benchmarks
- **CSGrs Integration**: 3D mesh support through the CSGrs crate for complex geometries

## Quick Start

### Prerequisites

- Rust 1.75 or later
- Cargo package manager

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/cfd-suite.git
cd cfd-suite

# Build the project
cargo build --release

# Run tests
cargo test

# Run a simple example
cargo run --example simple_pipe_flow
```

### Basic Usage

```rust
use cfd_suite::prelude::*;

// Create a simple 1D pipe flow simulation
let mut sim = Simulation1D::builder()
    .add_pipe(0.0, 1.0, 0.01)  // length=1m, diameter=1cm
    .set_fluid(Fluid::water())
    .set_boundary_conditions(
        BoundaryCondition::pressure_inlet(101325.0),
        BoundaryCondition::pressure_outlet(100000.0),
    )
    .build()?;

// Run simulation
let results = sim.run()?;

// Export results
results.export_csv("pipe_flow_results.csv")?;
```

## Architecture

The CFD suite is organized as a Rust workspace with the following crates:

- **`cfd-core`**: Core abstractions, plugin system, and common types
- **`cfd-math`**: Mathematical utilities and numerical methods
- **`cfd-io`**: File I/O operations (VTK, CSV, JSON, HDF5)
- **`cfd-mesh`**: Mesh handling and geometry operations
- **`cfd-1d`**: 1D solvers for pipe networks and microfluidic simulations
- **`cfd-2d`**: 2D solvers (FDM, FVM, LBM)
- **`cfd-3d`**: 3D solvers with CSGrs integration (FEM, spectral methods)
- **`cfd-validation`**: Validation framework and benchmark problems

### External Integrations

- **CSGrs**: Used in `cfd-3d` for 3D mesh handling and constructive solid geometry operations
- **scheme**: Used in `cfd-1d` for 2D schematic visualization of microfluidic networks (similar to electronic circuit design tools)
  - **Note**: The scheme integration currently requires nightly Rust due to unstable features. See `crates/cfd-1d/README.md` for details.

### Plugin System

The plugin architecture allows for easy extension:

```rust
use cfd_suite::plugin::{SimulationPlugin, PluginRegistry};

// Define a custom solver plugin
struct MyCustomSolver;

impl SimulationPlugin for MyCustomSolver {
    type Config = MyConfig;
    type State = MyState;
    type Output = MyOutput;
    
    fn initialize(&self, config: Self::Config) -> Result<Self::State> {
        // Initialize solver state
    }
    
    fn step(&self, state: &mut Self::State, dt: f64) -> Result<()> {
        // Perform one time step
    }
    
    fn output(&self, state: &Self::State) -> Self::Output {
        // Generate output
    }
}

// Register the plugin
PluginRegistry::register("my_solver", MyCustomSolver);
```

## Supported Simulations

### 1D Simulations
- **Pipe Networks**: Hagen-Poiseuille flow in complex networks
- **Microfluidics**: Channel-based microfluidic devices (MMFT-compatible)
- **Electrical Analogy**: Fast network solvers using circuit analogies

### 2D Simulations
- **Finite Difference Method (FDM)**: Structured grid simulations
- **Finite Volume Method (FVM)**: Conservative schemes for complex flows
- **Lattice Boltzmann Method (LBM)**: Mesoscopic approach for complex physics

### 3D Simulations
- **Finite Element Method (FEM)**: Unstructured mesh support
- **Spectral Methods**: High-accuracy simulations
- **CSGrs Integration**: Complex geometry handling via Constructive Solid Geometry

## Examples

### 1D Microfluidic Network
```rust
// Create a microfluidic T-junction
let network = NetworkBuilder::new()
    .add_channel("inlet", 100e-6, 50e-6, 1000e-6)  // width, height, length
    .add_channel("outlet1", 100e-6, 50e-6, 500e-6)
    .add_channel("outlet2", 100e-6, 50e-6, 500e-6)
    .connect("inlet", "junction")
    .connect("junction", "outlet1")
    .connect("junction", "outlet2")
    .build()?;
```

### 2D Lid-Driven Cavity
```rust
// Classic benchmark problem
let cavity = Simulation2D::lid_driven_cavity()
    .set_reynolds(1000.0)
    .set_grid_size(128, 128)
    .set_lid_velocity(1.0)
    .build()?;
```

### 3D Flow Around Obstacle
```rust
// Using CSGrs for geometry
use csgrs::prelude::*;

let obstacle = Mesh::sphere(0.1, 32, 16, None);
let domain = Mesh::cuboid(2.0, 1.0, 1.0, None)
    .difference(&obstacle.translate(1.0, 0.5, 0.5));

let simulation = Simulation3D::from_csg(domain)
    .set_inlet_velocity(vec3(1.0, 0.0, 0.0))
    .set_fluid_properties(Fluid::air())
    .build()?;
```

### 1D Microfluidic Network Example

```rust
use cfd_1d::prelude::*;
use cfd_core::prelude::*;

// Create a simple microfluidic network
let mut network = NetworkBuilder::new()
    .add_channel("ch1", 1e-3, 50e-6) // 1mm long, 50μm diameter
    .add_pump("pump1", PumpType::Pressure(1000.0)) // 1 kPa
    .add_junction("j1", JunctionType::TMixer)
    .connect("pump1", "ch1")
    .connect("ch1", "j1")
    .build()?;

// Set fluid properties
network.set_fluid(Fluid::water());

// Solve for steady-state flow
let solution = ElectricalAnalogySolver::new()
    .solve(&network)?;

// Export to 2D schematic (when scheme integration is available)
// let schematic = network.to_scheme()?;
// schematic.save("network.scheme")?;
```

## Validation

All implemented algorithms are validated against:

- **Analytical Solutions**: Poiseuille flow, Couette flow, Stokes flow
- **Benchmark Problems**: Lid-driven cavity, flow over cylinder, backward-facing step
- **Literature References**: Extensive comparison with published results

See the [validation report](docs/validation.md) for detailed comparisons.

## Performance

The suite is designed with performance in mind:

- **Zero-copy abstractions**: Minimal memory allocations
- **Iterator-based algorithms**: Leveraging Rust's iterator optimizations
- **Parallel execution**: Multi-threaded solvers where applicable
- **SIMD optimizations**: Vectorized operations for supported architectures

## Documentation

- [API Documentation](https://docs.rs/cfd-suite)
- [User Guide](docs/user_guide.md)
- [Developer Guide](docs/developer_guide.md)
- [Validation Report](docs/validation.md)

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Setup

```bash
# Install development dependencies
cargo install cargo-watch cargo-expand

# Run tests in watch mode
cargo watch -x test

# Check code quality
cargo clippy -- -W clippy::pedantic
cargo fmt --check
```

## Design Principles

This project adheres to the following principles:

- **SOLID**: Single Responsibility, Open/Closed, Liskov Substitution, Interface Segregation, Dependency Inversion
- **CUPID**: Composable, Unix philosophy, Predictable, Idiomatic, Domain-based
- **Clean Architecture**: Clear separation of concerns and dependencies
- **Zero-cost Abstractions**: Performance without compromise

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Inspired by the [MMFT Simulator](https://github.com/cda-tum/mmft-simulator) for 1D microfluidics
- 3D mesh support through [CSGrs](https://lib.rs/crates/csgrs)
- Mathematical operations powered by [nalgebra](https://nalgebra.org/)

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

## Roadmap

- [ ] GPU acceleration support
- [ ] Adaptive mesh refinement
- [ ] Turbulence modeling
- [ ] Multiphysics coupling
- [ ] Real-time visualization

See the [full roadmap](ROADMAP.md) for more details.