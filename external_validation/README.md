# CFD-rs External Validation Framework

This directory contains comprehensive validation tools that compare CFD-rs simulation results against established Python CFD packages and analytical solutions.

## Purpose

The external validation framework ensures CFD-rs produces correct results by:

1. **Comparing with reference implementations**: Python CFD packages with known correctness
2. **Validating against analytical solutions**: Mathematical exact solutions for simple cases
3. **Benchmarking against literature**: Published benchmark data (e.g., Ghia et al. 1982)
4. **Quantifying errors**: L2 norms, relative errors, convergence rates

## Reference Packages

### 1. [Python_CFD](https://github.com/DrZGan/Python_CFD)
Educational CFD implementations in Python:
- 1D finite difference methods
- 2D incompressible Navier-Stokes
- Pressure-velocity coupling (SIMPLE, SIMPLER)
- Channel and cavity flows

### 2. [cfd-comparison-python](https://github.com/pmocz/cfd-comparison-python)
Philip Mocz's comparison of numerical methods:
- **Finite Volume Method**: Godunov-type with minmod limiter
- **Spectral Method**: Fourier-based for periodic domains
- **Lattice Boltzmann Method**: D2Q9 lattice
- **Smoothed Particle Hydrodynamics**: Lagrangian particle method

### 3. [Fluidsim](https://fluidsim.readthedocs.io/)
High-performance spectral CFD framework:
- Pseudo-spectral Navier-Stokes solvers
- Turbulence simulations
- Fourier-based spatial discretization

## Validation Test Cases

### 1D Tests

#### Poiseuille Flow (Pipe Flow)
- **Analytical**: Hagen-Poiseuille equation
- **CFD-rs**: `PoiseuilleSolver1D` via pycfdrs
- **Reference**: Python_CFD channel flow
- **Validation**: Flow rate, velocity profile, pressure drop
- **Tolerance**: < 2% relative error

```
Q = π * R⁴ * ΔP / (8 * μ * L)
```

#### Bifurcation Flow
- **Analytical**: Mass conservation (Q_in = Q_out1 + Q_out2)
- **CFD-rs**: `BifurcationSolver` via pycfdrs
- **Physics**: Murray's law for optimal branching
- **Validation**: Flow split ratio, pressure continuity
- **Tolerance**: < 1e-10 mass conservation error

### 2D Tests

#### Poiseuille Flow (Channel)
- **Analytical**: Parabolic velocity profile
- **CFD-rs**: `Poiseuille2DSolver` via pycfdrs
- **Reference**: cfd-comparison-python FVM
- **Validation**: Velocity profile, mean velocity
- **Tolerance**: < 5% relative error

```
u(y) = (ΔP/2μL) * y * (H - y)
```

#### Lid-Driven Cavity
- **Reference**: Ghia et al. (1982) benchmark data
- **CFD-rs**: `CavitySolver2D` via pycfdrs
- **Cases**: Re = 100, 400, 1000, 3200, 5000, 7500, 10000
- **Validation**: Centerline velocity profiles (U and V)
- **Tolerance**: < 5% L2 error vs Ghia data

#### Venturi Flow
- **Analytical**: Bernoulli equation
- **CFD-rs**: `VenturiSolver2D` via pycfdrs
- **Reference**: Bernoulli pressure coefficient
- **Validation**: Throat pressure, pressure recovery
- **Tolerance**: < 10% (accounts for viscous losses)

```
Cp = (P_throat - P_inlet) / (0.5 * ρ * U_inlet²) = 1 - (A_inlet/A_throat)²
```

#### Serpentine Flow
- **Physics**: Dean vortices in curved channels
- **CFD-rs**: `Serpentine2DSolver` via pycfdrs
- **Validation**: Secondary flow patterns, mixing efficiency
- **Reference**: Dean number correlations

### 3D Tests

#### Poiseuille Flow (Circular Pipe)
- **Analytical**: Hagen-Poiseuille
- **CFD-rs**: `Poiseuille3DSolver` via pycfdrs
- **Validation**: Flow rate, velocity profile
- **Tolerance**: < 5% relative error

#### Bifurcation Flow
- **Physics**: Murray's law, 3D flow patterns
- **CFD-rs**: `Bifurcation3DSolver` via pycfdrs
- **Validation**: Flow split, wall shear stress distribution
- **Blood models**: Casson, Carreau-Yasuda

## Directory Structure

```
external_validation/
├── README.md                          # This file
├── python_cfd/                        # Python_CFD reference clone
├── cfd_comparison/                    # cfd-comparison-python reference
│   └── finitevolume_ref.py           # FVM implementation
├── fluidsim/                          # Fluidsim validation scripts
├── scripts/
│   ├── validation_runner.py          # Main validation orchestrator
│   └── compare_with_reference.py     # Direct comparison tool
└── results/                           # Generated validation reports
    ├── validation_report.json
    └── validation_results.png
```

## Usage

### Prerequisites

```bash
# Build pycfdrs Python bindings
cd crates/pycfdrs
maturin develop --release

# Install Python dependencies
pip install numpy matplotlib

# Optional: Install reference packages
pip install fluidsim
```

### Run All Validations

```bash
python external_validation/scripts/validation_runner.py --test all
```

### Run Specific Test

```bash
# 1D validations only
python external_validation/scripts/validation_runner.py --test 1d

# 2D validations only
python external_validation/scripts/validation_runner.py --test 2d

# Individual tests
python external_validation/scripts/validation_runner.py --test poiseuille
python external_validation/scripts/validation_runner.py --test bifurcation
python external_validation/scripts/validation_runner.py --test cavity
python external_validation/scripts/validation_runner.py --test venturi
```

### Compare with Reference Implementations

```bash
python external_validation/scripts/compare_with_reference.py --test all
```

## Interpreting Results

### Validation Report Format

```json
{
  "test_name": "Poiseuille 1D Flow Rate",
  "test_type": "1d",
  "cfrs_value": 1.23e-9,
  "reference_value": 1.25e-9,
  "absolute_error": 2.0e-11,
  "relative_error": 0.016,
  "tolerance": 0.02,
  "passed": true,
  "units": "m³/s",
  "description": "Flow rate in circular pipe vs analytical Hagen-Poiseuille"
}
```

### Success Criteria

| Test Category | Tolerance | Metric |
|--------------|-----------|--------|
| 1D Poiseuille | 2% | Relative error in flow rate |
| 1D Bifurcation | 1e-10 | Mass conservation error |
| 2D Poiseuille | 5% | Relative error in mean velocity |
| 2D Cavity | 5% | L2 error vs Ghia benchmark |
| 2D Venturi | 10% | Pressure coefficient error |
| 3D Poiseuille | 5% | Relative error in flow rate |

## Blood Flow Validation

### Non-Newtonian Models

#### Casson Model
```python
blood = pycfdrs.CassonBlood()
# Yield stress: τ_y = 0.0056 Pa
# Infinite-shear viscosity: μ_∞ = 0.00345 Pa·s
```

#### Carreau-Yasuda Model
```python
blood = pycfdrs.CarreauYasudaBlood()
# Zero-shear: μ_0 = 0.056 Pa·s
# Infinite-shear: μ_∞ = 0.00345 Pa·s
# Relaxation time: λ = 3.313 s
# Power-law index: n = 0.3568
```

### Physiological Validation

- **Shear rates**: 10-1000 s⁻¹ (physiological range)
- **Wall shear stress**: 1-5 Pa in capillaries, 0.5-3 Pa in arteries
- **Fåhræus-Lindqvist effect**: Viscosity reduction in vessels < 300 μm

## Literature References

### Analytical Solutions
1. **Hagen-Poiseuille**: Viscous flow in circular pipes
2. **Bernoulli**: Inviscid flow, pressure-velocity relationship
3. **Dean**: Secondary flow in curved pipes

### Benchmark Data
1. **Ghia et al. (1982)**: Lid-driven cavity benchmark
   - Journal of Computational Physics, 48(3):387-411
2. **Murray (1926)**: Optimal branching (Murray's law
   - Proc. Natl. Acad. Sci. USA, 12(3):207-214

### Blood Rheology
1. **Merrill et al. (1969)**: Casson model parameters
2. **Cho & Kensey (1991)**: Carreau-Yasuda for blood
3. **Pries et al. (1992)**: Fåhræus-Lindqvist effect

## Troubleshooting

### pycfdrs Not Found
```bash
cd crates/pycfdrs
maturin develop --release
```

### Matplotlib Not Available
Validation runs without plotting. Install for visualizations:
```bash
pip install matplotlib
```

### Numerical Instability
- Reduce time step
- Increase grid resolution
- Check Reynolds number (avoid stiffness)

## Contributing

To add new validation tests:

1. Implement test in `validation_runner.py`
2. Add corresponding method to `ReferenceComparator`
3. Document expected results and tolerances
4. Update this README with test description

## License

The validation framework and CFD-rs are MIT/Apache-2.0 licensed.
Reference implementations may have their own licenses:
- cfd-comparison-python: GPL-3.0
- Python_CFD: Check repository
- Fluidsim: Check repository
