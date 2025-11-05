# MPI Setup Guide for CFD Suite

This guide provides comprehensive instructions for enabling MPI (Message Passing Interface) support in the CFD simulation suite for distributed-memory parallel computing.

## Overview

The CFD suite includes MPI-based domain decomposition capabilities that enable:

- **Scalability**: Solve problems with millions of cells across multiple compute nodes
- **Load Balancing**: Automatic domain decomposition with optimal process distribution
- **Ghost Cell Communication**: Efficient boundary data exchange between processes
- **Parallel I/O**: Distributed data input/output operations
- **Collective Operations**: Global reductions and synchronizations

## System Requirements

### Linux (Ubuntu/Debian)

```bash
# Install OpenMPI development libraries
sudo apt-get update
sudo apt-get install libopenmpi-dev openmpi-bin

# Verify installation
mpirun --version
```

### macOS (Homebrew)

```bash
# Install OpenMPI
brew install open-mpi

# Verify installation
mpirun --version
```

### Windows (Microsoft MPI)

1. Download Microsoft MPI v10.1.2 from:
   https://www.microsoft.com/en-us/download/details.aspx?id=100593

2. Install both `msmpisdk.msi` and `MSMpiSetup.exe`

3. Add to PATH:
   ```
   C:\Program Files\Microsoft MPI\Bin\
   ```

4. Verify installation:
   ```cmd
   mpiexec.exe -help
   ```

## Compilation

### Enable MPI Feature

Compile with MPI support enabled:

```bash
# Build with MPI support
cargo build --features mpi

# Run tests with MPI
cargo test --features mpi

# Run specific MPI tests
cargo test --features mpi mpi_integration
```

### Feature Configuration

The MPI feature is defined in `Cargo.toml`:

```toml
[features]
mpi = ["cfd-core/mpi", "dep:mpi"]
```

## Architecture Overview

### Domain Decomposition Strategies

The MPI implementation supports multiple decomposition strategies:

1. **Simple 1D**: Strip decomposition along x-direction
2. **Cartesian 2D**: Block decomposition in x-y plane
3. **Recursive Bisection**: Advanced load balancing (future)
4. **METIS**: Graph partitioning (future)

### Communication Patterns

- **Ghost Cells**: Halo exchange for boundary conditions
- **Collective Operations**: Reductions, broadcasts, barriers
- **Point-to-Point**: Direct process communication

## Usage Examples

### Basic MPI Initialization

```rust
use cfd_core::compute::mpi::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Initialize MPI
    let universe = MpiUniverse::new()?;
    let world = universe.world();

    println!("Running on {} processes", world.size());

    // Your CFD computation here
    run_parallel_cfd(&world)?;

    Ok(())
}
```

### Domain Decomposition

```rust
use cfd_core::compute::mpi::*;

// Create global domain
let global_extents = GlobalExtents::new_2d(nx, ny, bounds);

// Setup domain decomposition
let decomp = DomainDecomposition::new(
    global_extents,
    &communicator,
    DecompositionStrategy::Cartesian2D,
)?;

// Create distributed grid
let dist_grid = DistributedGrid::new(
    global_extents,
    &communicator,
    DecompositionStrategy::Cartesian2D,
)?;
```

### Parallel CFD Solver

```rust
// Example parallel CFD workflow
let solver = ParallelCfdSolver::new(dist_grid, config)?;

// Time stepping loop
for step in 0..max_steps {
    // Update ghost cells
    dist_grid.update_ghost_cells()?;

    // Solve momentum equation
    solver.solve_momentum()?;

    // Solve pressure correction
    solver.solve_pressure()?;

    // Correct velocities
    solver.correct_velocity()?;

    // Global convergence check
    let residual = dist_grid.global_reduce_max(local_residual);
    if residual < tolerance {
        break;
    }
}
```

## Performance Tuning

### Process Distribution

Optimal process counts depend on problem size:

- **Small problems** (< 10⁴ cells): 1-4 processes
- **Medium problems** (10⁴-10⁶ cells): 4-16 processes
- **Large problems** (>10⁶ cells): 16-256+ processes

### Load Balancing

The decomposition algorithms attempt to:
- Minimize surface-to-volume ratio
- Balance computational load
- Optimize communication patterns

### Communication Optimization

- Use non-blocking communication for overlapping compute/communicate
- Minimize ghost cell layers (typically 1-2 layers)
- Optimize message sizes for network hardware

## Troubleshooting

### Common Issues

#### MPI Library Not Found

**Error**: `Could not find MPI library`

**Solution**: Ensure MPI development libraries are installed and in PATH

#### Compilation Errors

**Error**: `undefined reference to MPI_*`

**Solution**: Check that MPI feature is enabled and libraries are linked

#### Runtime Errors

**Error**: `MPI_COMM_WORLD not initialized`

**Solution**: Ensure `MpiUniverse::new()` is called before MPI operations

### Debugging MPI Programs

Run with debugging output:

```bash
# Linux/macOS
mpirun -np 4 --oversubscribe cargo run --features mpi --bin your_program

# Windows
mpiexec -n 4 cargo run --features mpi --bin your_program
```

### Performance Profiling

Use MPI profiling tools:

```bash
# Intel VTune (with MPI)
mpirun -np 4 vtune -collect hotspots ./your_program

# HPCToolkit
mpirun -np 4 hpcrun ./your_program
```

## Integration with CFD Solvers

The MPI implementation integrates with existing CFD solvers:

### Momentum Solver

```rust
let momentum_solver = ParallelMomentumSolver::new(&dist_grid, &config)?;
momentum_solver.solve_parallel(&mut fields)?;
```

### Energy Solver

```rust
let energy_solver = ParallelEnergySolver::new(&dist_grid)?;
energy_solver.solve_parallel(&mut temperature)?;
```

### Turbulence Models

```rust
let turbulence_solver = ParallelTurbulenceSolver::new(&dist_grid, model_type)?;
turbulence_solver.solve_parallel(&mut fields)?;
```

## Testing

### Unit Tests

Run MPI-specific tests:

```bash
cargo test --features mpi mpi_unit_tests
```

### Integration Tests

```bash
# Run with 4 processes
mpirun -np 4 cargo test --features mpi mpi_integration
```

### Validation Tests

```bash
# Compare serial vs parallel results
mpirun -np 4 cargo test --features mpi validation_parallel
```

## Future Enhancements

### Planned Features

1. **Advanced Load Balancing**
   - Dynamic repartitioning during simulation
   - Zoltan integration for complex geometries

2. **Parallel I/O**
   - HDF5 parallel I/O
   - ADIOS2 support
   - Checkpoint/restart capabilities

3. **Adaptive Meshes**
   - Parallel mesh refinement
   - Dynamic load balancing for AMR

4. **GPU Acceleration**
   - MPI + CUDA/OpenCL integration
   - Multi-GPU per node support

5. **Fault Tolerance**
   - Checkpoint/restart mechanisms
   - Process failure recovery

## References

- **MPI Standard**: https://www.mpi-forum.org/
- **OpenMPI**: https://www.open-mpi.org/
- **Microsoft MPI**: https://docs.microsoft.com/en-us/message-passing-interface/
- **Domain Decomposition Literature**:
  - Smith, Bjørstad, Gropp: "Domain Decomposition"
  - Quarteroni, Valli: "Domain Decomposition Methods for PDEs"

## Support

For MPI-related issues:

1. Check this documentation
2. Verify MPI installation
3. Test with simple MPI examples
4. Check CFD suite issue tracker

The MPI implementation follows best practices from high-performance computing and is designed for production CFD applications requiring distributed computing capabilities.
