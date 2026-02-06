# pycfdrs

Python bindings for CFD-rs: High-performance computational fluid dynamics in Rust.

## Installation

```bash
cd crates/pycfdrs
maturin develop  # Development mode
# or
maturin build --release  # Production build
pip install target/wheels/pycfdrs-*.whl
```

## Quick Start

```python
import pycfdrs

# Create bifurcation solver
solver = pycfdrs.BifurcationSolver(
    d_parent=100e-6,
    d_daughter1=80e-6,
    d_daughter2=80e-6
)

# Solve with blood
result = solver.solve(
    flow_rate=3e-8,
    blood_type="casson"
)

print(f"Pressure drop: {result.dp_1:.2f} Pa")
print(f"Wall shear stress: {result.wss_1:.2f} Pa")
```

## Features

- 1D bifurcation/trifurcation network solvers
- 2D Poiseuille and Venturi flow solvers
- 3D FEM bifurcation solver
- Non-Newtonian blood rheology models (Casson, Carreau-Yasuda)
- Validated against analytical solutions and FEniCS

## Documentation

See `validation/` directory for comprehensive validation examples.
