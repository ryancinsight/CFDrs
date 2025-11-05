# CFD Suite MPI Deployment Guide

This guide provides comprehensive instructions for deploying the CFD suite with MPI parallelization for production CFD simulations.

## Table of Contents

1. [System Requirements](#system-requirements)
2. [MPI Installation](#mpi-installation)
3. [CFD Suite Installation](#cfd-suite-installation)
4. [Configuration](#configuration)
5. [Running Simulations](#running-simulations)
6. [Performance Optimization](#performance-optimization)
7. [Monitoring and Troubleshooting](#monitoring-and-troubleshooting)
8. [Scaling Guidelines](#scaling-guidelines)

## System Requirements

### Hardware Requirements

#### Minimum Requirements
- **CPU**: 4 cores (Intel/AMD with AVX2 support)
- **RAM**: 8 GB per node
- **Storage**: 10 GB available space
- **Network**: 1 GbE interconnect

#### Recommended Production Setup
- **CPU**: 16-32 cores per node (Intel Xeon/AMD EPYC)
- **RAM**: 64-128 GB per node
- **Storage**: 100 GB SSD per node + shared storage
- **Network**: 10 GbE or higher (Infiniband preferred for >64 cores)

### Software Requirements

#### Operating Systems
- **Linux**: Ubuntu 20.04+, CentOS 7+, RHEL 8+ (recommended)
- **macOS**: 11.0+ (development only)
- **Windows**: WSL2 with Ubuntu (development only)

#### Required Software
- **Rust**: 1.70+ (`rustc --version`)
- **MPI Implementation**: OpenMPI 4.1+ or MPICH 3.4+
- **C Compiler**: GCC 9+ or Clang 11+
- **CMake**: 3.16+ (for some dependencies)
- **Git**: 2.25+

## MPI Installation

### Linux (Ubuntu/Debian)

```bash
# Update package list
sudo apt update

# Install OpenMPI
sudo apt install -y openmpi-bin openmpi-doc libopenmpi-dev

# Verify installation
mpirun --version
mpicc --version
```

### Linux (CentOS/RHEL)

```bash
# Install OpenMPI
sudo yum install -y openmpi openmpi-devel

# Or for newer versions
sudo dnf install -y openmpi openmpi-devel

# Verify installation
mpirun --version
```

### Linux (From Source - Recommended for Performance)

```bash
# Download and extract OpenMPI
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.4.tar.gz
tar -xzf openmpi-4.1.4.tar.gz
cd openmpi-4.1.4

# Configure and build
./configure --prefix=/opt/openmpi-4.1.4 \
            --enable-mpi-cxx \
            --enable-mpi-fortran=no \
            --with-cma \
            --with-verbs
make -j$(nproc)
sudo make install

# Add to PATH
echo 'export PATH=/opt/openmpi-4.1.4/bin:$PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=/opt/openmpi-4.1.4/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc
```

### macOS (Using Homebrew)

```bash
# Install OpenMPI
brew install open-mpi

# Verify installation
mpirun --version
```

### Windows (WSL2)

```bash
# Install Ubuntu on WSL2, then follow Linux instructions
# Note: Windows native MPI support is limited
```

### Alternative: MPICH

```bash
# MPICH installation (alternative to OpenMPI)
sudo apt install -y mpich libmpich-dev

# Verify installation
mpirun --version
mpicc --version
```

## CFD Suite Installation

### Clone Repository

```bash
# Clone the CFD suite repository
git clone https://github.com/your-org/cfd-suite.git
cd cfd-suite
```

### Build with MPI Support

```bash
# Build with MPI feature enabled
cargo build --release --features mpi

# Verify MPI support is enabled
cargo check --features mpi

# Run tests with MPI support
cargo test --features mpi
```

### Build Documentation

```bash
# Generate documentation
cargo doc --features mpi --open
```

### Install System-wide (Optional)

```bash
# Install binaries system-wide
cargo install --path crates/cfd-suite --features mpi

# Verify installation
cfd-suite --help
```

## Configuration

### Environment Variables

```bash
# Set MPI implementation (if multiple installed)
export MPI_ROOT=/opt/openmpi-4.1.4

# Set thread affinity for better performance
export OMP_NUM_THREADS=1  # Disable OpenMP threading in MPI processes
export MPI_THREAD_MULTIPLE=0

# Set memory allocation (if needed)
export MALLOC_ARENA_MAX=1

# Enable core dumps for debugging
ulimit -c unlimited
```

### MPI Hostfile Configuration

Create `/etc/hostsfile` or `hostfile` in working directory:

```
# Hostfile format: hostname slots=max_processes_per_node
node01 slots=16
node02 slots=16
node03 slots=16
node04 slots=16
```

### Batch System Configuration (SLURM)

Example SLURM job script (`submit_job.slurm`):

```bash
#!/bin/bash
#SBATCH --job-name=cfd_simulation
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --account=my_project

# Load modules
module load openmpi/4.1.4
module load gcc/9.3.0

# Set environment
export OMP_NUM_THREADS=1

# Run simulation
mpirun -np $SLURM_NTASKS \
       --hostfile $SLURM_NODELIST \
       ./target/release/cfd-suite \
       --config simulation_config.toml \
       --output results/
```

## Running Simulations

### Basic MPI Execution

```bash
# Run on local machine with 4 processes
mpirun -np 4 ./target/release/cfd-suite --config config.toml

# Run across multiple nodes
mpirun -np 64 --hostfile hostfile ./target/release/cfd-suite --config config.toml

# Run with specific process placement
mpirun -np 16 --map-by node ./target/release/cfd-suite --config config.toml
```

### Configuration File Example

Create `simulation_config.toml`:

```toml
[domain]
nx = 1000
ny = 1000
x_min = 0.0
x_max = 1.0
y_min = 0.0
y_max = 1.0

[physics]
reynolds_number = 1000.0
time_step = 0.001
max_time = 10.0

[solver]
tolerance = 1e-6
max_iterations = 1000
preconditioner = "block_jacobi"

[mpi]
decomposition_strategy = "cartesian_2d"
load_balance_threshold = 1.2
enable_adaptive_mesh = true

[output]
vtk_output = true
output_frequency = 100
checkpoint_frequency = 1000
output_directory = "results/"
```

### Advanced MPI Options

```bash
# Bind processes to specific cores
mpirun -np 16 --bind-to core ./target/release/cfd-suite --config config.toml

# Use specific network interface
mpirun -np 64 --mca btl_tcp_if_include ib0 ./target/release/cfd-suite --config config.toml

# Enable debugging output
mpirun -np 4 --debug-daemons ./target/release/cfd-suite --config config.toml

# Profile MPI communication
mpirun -np 16 --mca coll_base_verbose 1 ./target/release/cfd-suite --config config.toml
```

## Performance Optimization

### Process Placement

```bash
# Optimal placement for multi-socket systems
mpirun -np 32 --map-by socket:PE=8 ./target/release/cfd-suite --config config.toml

# NUMA-aware placement
mpirun -np 16 --bind-to numa ./target/release/cfd-suite --config config.toml
```

### Memory Optimization

```bash
# Use huge pages for better memory performance
echo 1 > /proc/sys/vm/nr_hugepages  # Enable huge pages

# Run with huge page support
mpirun -np 16 --mca mpi_yield_when_idle 1 \
       ./target/release/cfd-suite --config config.toml
```

### Network Optimization

```bash
# Use high-performance interconnect
mpirun -np 64 --mca btl openib,self ./target/release/cfd-suite --config config.toml

# Optimize collective operations
mpirun -np 128 --mca coll_tuned_use_dynamic_rules 1 \
       ./target/release/cfd-suite --config config.toml
```

### Load Balancing

```bash
# Enable dynamic load balancing
export CFD_LOAD_BALANCE_THRESHOLD=1.2
export CFD_MIN_CELLS_PER_PROCESS=1000

mpirun -np 64 ./target/release/cfd-suite \
       --config config.toml \
       --enable-load-balancing
```

## Monitoring and Troubleshooting

### Performance Monitoring

```bash
# Monitor MPI processes
mpirun -np 16 --mca mpi_show_mpit_debug 1 ./target/release/cfd-suite --config config.toml

# Use performance validator
mpirun -np 16 ./target/release/performance_benchmark \
       --cores 1,2,4,8,16 \
       --output scaling_results.json
```

### Common Issues and Solutions

#### 1. MPI Initialization Errors

**Error**: `MPI_Init failed`
**Solution**:
```bash
# Check MPI installation
mpirun --version

# Check network connectivity
ping other_nodes

# Verify SSH key setup for passwordless login
ssh-keygen -t rsa
ssh-copy-id other_nodes
```

#### 2. Memory Issues

**Error**: `Out of memory`
**Solution**:
```bash
# Increase memory limits
ulimit -m unlimited
ulimit -v unlimited

# Reduce problem size or increase nodes
# Check memory usage per process
```

#### 3. Communication Deadlocks

**Error**: Simulation hangs during communication
**Solution**:
```bash
# Enable timeout
mpirun -np 16 --mca mpi_abort_delay 10 ./target/release/cfd-suite --config config.toml

# Check for communication pattern issues
# Verify ghost cell exchange implementation
```

#### 4. Load Imbalance

**Symptoms**: Poor parallel efficiency, high communication overhead
**Solution**:
```bash
# Enable load balancing
export CFD_LOAD_BALANCE_ENABLED=1

# Adjust decomposition strategy
# Use adaptive mesh refinement
```

### Debugging Tools

```bash
# MPI debugging
mpirun -np 4 --debug ./target/release/cfd-suite --config config.toml

# Valgrind for memory debugging (single process)
valgrind --tool=memcheck ./target/release/cfd-suite --config config.toml

# Performance profiling
mpirun -np 16 --mca mpi_paffinity_alone 1 \
       perf record ./target/release/cfd-suite --config config.toml
```

## Scaling Guidelines

### Optimal Core Counts

| Problem Size | Optimal Cores | Efficiency Target |
|-------------|---------------|-------------------|
| 100×100     | 1-4          | >95%             |
| 1K×1K       | 4-16         | >90%             |
| 10K×10K     | 16-64        | >85%             |
| 100K×100K   | 64-256       | >80%             |

### Network Requirements

| Core Count | Network Type | Bandwidth/Core |
|------------|--------------|----------------|
| 1-16      | 1 GbE       | 100 Mb/s      |
| 16-64     | 10 GbE      | 1 Gb/s        |
| 64-256    | 10 GbE/IB   | 10 Gb/s       |
| 256+      | Infiniband  | 100 Gb/s      |

### Memory Scaling

- **Per Core**: 1-4 GB for typical CFD problems
- **Total**: Problem size × 8 bytes × safety factor (2-4)
- **Distribution**: Evenly distribute across nodes

### Checkpoint Recommendations

| Simulation Size | Checkpoint Frequency | Storage Requirements |
|----------------|---------------------|---------------------|
| Small (<1GB)   | Every 100 steps    | 10-100 GB          |
| Medium (1-10GB)| Every 500 steps    | 100-1TB            |
| Large (>10GB)  | Every 1000 steps   | 1-10TB             |

## Production Deployment Checklist

- [ ] MPI installation verified
- [ ] Network connectivity tested
- [ ] Hostfile configured
- [ ] CFD suite built with MPI support
- [ ] Configuration files validated
- [ ] Test run on single node successful
- [ ] Scaling test completed
- [ ] Monitoring tools configured
- [ ] Backup and recovery procedures in place
- [ ] Documentation accessible to users

## Support and Resources

### Documentation
- [CFD Suite User Guide](USER_GUIDE.md)
- [API Documentation](https://docs.rs/cfd-suite/)
- [Performance Validation Guide](PERFORMANCE_VALIDATION.md)

### Community Resources
- GitHub Issues: Report bugs and request features
- Discussions: Ask questions and share experiences
- Wiki: Community-contributed guides and tutorials

### Performance Tuning
- [HPC Best Practices](https://www.mpi-forum.org/)
- [MPI Performance Analysis](https://www.open-mpi.org/faq/)
- [CFD Scaling Studies](https://www.cfd-online.com/)

---

*This deployment guide is maintained as part of the CFD suite documentation. Please report any issues or suggest improvements via GitHub issues.*
