# cfd-python

Python bindings for CFD-rs: High-performance computational fluid dynamics in Rust, exposed via [PyO3](https://pyo3.rs/) and [Maturin](https://www.maturin.rs/).

## Installation

```bash
cd crates/cfd-python
pip install maturin
maturin develop              # Development mode (editable install)
# or
maturin build --release      # Production wheel
pip install target/wheels/cfd_python-*.whl
```

Optional dependencies:

```bash
pip install "cfd-python[dev]"         # pytest, numpy, matplotlib, pandas
pip install "cfd-python[validation]"  # adds scipy and jupyter
```

## Quick Start

```python
import cfd_python

# 1D bifurcation solver
bifurc = cfd_python.BifurcationSolver(
    d_parent=100e-6,
    d_daughter1=80e-6,
    d_daughter2=80e-6,
)
blood = cfd_python.CassonBlood()
result = bifurc.solve(flow_rate=3e-8, pressure=40.0, blood=blood)
print(f"Pressure drop 1: {result.dp_1:.2f} Pa")
print(f"Wall shear rate 1: {result.gamma_1:.2f} s⁻¹")
```

## Available Classes

### 1D Solvers
| Class | Description |
|---|---|
| `BifurcationSolver` / `BifurcationResult` | Single bifurcation network |
| `TrifurcationSolver` / `TrifurcationResult` | Trifurcation network |
| `SerpentineSolver1D` / `SerpentineResult1D` | Serpentine channel resistance |
| `VenturiSolver1D` / `VenturiResult1D` | Venturi constriction resistance |

### 2D Solvers
| Class | Description |
|---|---|
| `PoiseuilleSolver` / `PoiseuilleResult` / `PoiseuilleConfig` | Poiseuille pipe flow |
| `Poiseuille2DSolver` / `Poiseuille2DResult` | 2D Poiseuille flow |
| `VenturiSolver2D` / `VenturiResult2D` | 2D Venturi constriction |
| `TrifurcationSolver2D` / `TrifurcationResult2D` | 2D trifurcation |
| `BifurcationSolver2D` / `BifurcationResult2D` | 2D bifurcation |
| `CavitySolver2D` / `CavityResult2D` | Lid-driven cavity (Ghia benchmark) |

### 3D Solvers
| Class | Description |
|---|---|
| `Bifurcation3DSolver` / `Bifurcation3DResult` | 3D FEM bifurcation |
| `Trifurcation3DSolver` / `Trifurcation3DResult` | 3D FEM trifurcation |
| `Poiseuille3DSolver` / `Poiseuille3DResult` | 3D Poiseuille flow |
| `Venturi3DSolver` / `Venturi3DResult` | 3D Venturi constriction |
| `Serpentine3DSolver` / `Serpentine3DResult` | 3D serpentine channel |

### Blood Rheology Models
| Class | Description |
|---|---|
| `CassonBlood` | Casson non-Newtonian model |
| `CarreauYasudaBlood` | Carreau-Yasuda shear-thinning model |
| `CrossBlood` | Cross model |
| `FahraeuasLindqvist` | Fåhraeus-Lindqvist effect |

### Pulsatile Flow
| Class | Description |
|---|---|
| `WomersleyNumber` | Womersley number calculation |
| `WomersleyProfile` | Velocity profile for pulsatile flow |
| `WomersleyFlow` | Full pulsatile flow solution |

## Features

- 1D/2D/3D CFD solvers unified under a single Python module
- Non-Newtonian blood rheology (Casson, Carreau-Yasuda, Cross, Fåhraeus-Lindqvist)
- Pulsatile Womersley flow
- Lid-driven cavity benchmark (Ghia et al.)
- Validated against analytical solutions

## Documentation

See the `crates/cfd-validation/` directory in the monorepo for comprehensive validation examples.
