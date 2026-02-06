# CFD-RS Validation Suite

This directory contains comprehensive validation scripts and examples for the cfd-rs CFD library, comparing results against analytical solutions, literature data, and established Python CFD packages.

## Validation Approach

Our validation strategy follows ASME V&V 20-2009 guidelines:

1. **Code Verification**: Compare against analytical solutions (method of manufactured solutions)
2. **Calculation Verification**: Grid convergence studies using Richardson extrapolation
3. **Validation**: Compare against experimental data and literature benchmarks

## Validation Categories

### 1D Hemodynamics

| Test Case | Validation Type | Literature Reference | Status |
|-----------|----------------|---------------------|--------|
| Poiseuille Flow (Casson) | Analytical | Merrill et al. (1969) | ✓ PASS |
| Murray's Law | Physical Law | Murray (1926) | ✓ PASS |
| Fåhræus-Lindqvist Effect | Empirical | Pries et al. (1992) | ✓ PASS |
| Womersley Flow | Analytical | Womersley (1955) | ✓ PASS |

### 2D Flows

| Test Case | Validation Type | Literature Reference | Status |
|-----------|----------------|---------------------|--------|
| Poiseuille Flow | Analytical | Hagen-Poiseuille | ✓ PASS |
| Lid-Driven Cavity | Benchmark | Ghia et al. (1982) | ✓ PASS |
| Venturi Flow | Physical Law | ISO 5167-1:2003 | ✓ PASS |
| Bifurcation Flow | Physical Law | Murray (1926) | ✓ PASS |

### Blood Rheology Models

| Model | Validation Type | Literature Reference | Status |
|-------|----------------|---------------------|--------|
| Casson | Literature Data | Merrill et al. (1969) | ✓ PASS |
| Carreau-Yasuda | Literature Data | Cho & Kensey (1991) | ✓ PASS |
| Cross | Literature Data | Cross (1965) | ✓ PASS |
| Fåhræus-Lindqvist | Empirical | Pries et al. (1992) | ✓ PASS |

## Python Comparison Scripts

### `python_cfd_comparison.py`

Comprehensive validation comparing cfd-rs (via pycfdrs) against analytical solutions and literature data.

**Usage:**
```bash
# Install pycfdrs first
cd crates/pycfdrs
maturin develop

# Run validation
python validation/python_cfd_comparison.py
```

**Features:**
- 1D Poiseuille flow with Casson blood model
- Murray's law validation for bifurcations
- 2D Poiseuille analytical comparison
- Venturi flow (Bernoulli equation)
- Blood rheology model validation
- Automated report generation (JSON)

### `cross_package_comparison.py`

Cross-validation against established Python CFD packages:
- [DrZGan/Python_CFD](https://github.com/DrZGan/Python_CFD)
- [pmocz/cfd-comparison-python](https://github.com/pmocz/cfd-comparison-python)
- [fluidsim](https://fluidsim.readthedocs.io/)

**Usage:**
```bash
pip install fluidsim numpy matplotlib scipy
python validation/cross_package_comparison.py
```

**Features:**
- Lid-driven cavity (Ghia et al. 1982 benchmark)
- Poiseuille flow analytical validation
- Velocity profile comparisons
- L2 error computation

## Rust Examples

### 1D Blood Flow Validation

```bash
cargo run --example blood_flow_1d_validation --no-default-features
```

Validates:
- Poiseuille flow with Casson model (Merrill 1969)
- Murray's law for symmetric bifurcations
- Asymmetric bifurcation flow split (Caro 1978)
- Fåhræus-Lindqvist effect (Pries 1992)

### 2D Bifurcation Validation

```bash
cargo run --example bifurcation_2d_blood_validation --no-default-features
```

Validates:
- Symmetric bifurcation (Murray's law)
- Asymmetric bifurcation (carotid geometry)
- Microvascular flow with Fåhræus-Lindqvist
- Wall shear stress calculations

### Venturi Flow Validation

```bash
cargo run --example venturi_blood_flow_validation --no-default-features
```

Validates:
- Bernoulli equation (inviscid limit)
- Pressure recovery coefficient (ISO 5167)
- Blood shear-thinning effects
- Carreau-Yasuda model

## Literature References

### Blood Rheology
1. **Merrill, E.W. et al. (1969)**. "Pressure-flow relations of human blood in hollow fibers at low flow rates". *J. Appl. Physiol.* 27(1):93-98.

2. **Cho, Y.I. & Kensey, K.R. (1991)**. "Effects of the non-Newtonian viscosity of blood on flows in a diseased arterial vessel". *Biorheology* 28(3-4):241-262.

3. **Pries, A.R. et al. (1992)**. "Blood viscosity in tube flow: dependence on diameter and hematocrit". *Am. J. Physiol.* 263(6):H1770-H1778.

### Hemodynamics
4. **Murray, C.D. (1926)**. "The Physiological Principle of Minimum Work". *Proc. Natl. Acad. Sci.* 12(3):207-214.

5. **Caro, C.G. et al. (1978)**. "The Mechanics of the Circulation". Oxford University Press.

6. **Womersley, J.R. (1955)**. "Method for the calculation of velocity, rate of flow and viscous drag in arteries when the pressure gradient is known". *J. Physiol.* 127(3):553-563.

### CFD Benchmarks
7. **Ghia, U.K.N.G. et al. (1982)**. "High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method". *J. Comput. Phys.* 48(3):387-411.

8. **ISO 5167-1:2003**. "Measurement of fluid flow by means of pressure differential devices".

## Validation Criteria

### Tolerances

| Test Type | Tolerance | Rationale |
|-----------|-----------|-----------|
| Analytical solutions | 1% | Round-off error, discretization |
| Physical laws | 5% | Measurement uncertainty |
| Literature data | 10-50% | Model variation, experimental error |
| Cross-package | 5% | Different numerical methods |

### Acceptance Criteria

A validation case passes if:
1. Error < specified tolerance
2. Physical trends are correct (e.g., viscosity decreases with shear rate)
3. Limiting behavior is correct (e.g., μ→μ_∞ as γ̇→∞)

## Report Generation

Validation reports are generated in JSON format:

```json
{
  "timestamp": "2024-01-15T10:30:00",
  "total_tests": 8,
  "passed_tests": 8,
  "cases": [
    {
      "name": "Poiseuille Flow (Casson)",
      "dimension": "1D",
      "test_type": "Analytical",
      "passed": true,
      "error_metric": 0.005,
      "tolerance": 0.01,
      "rust_value": 123.45,
      "reference_value": 123.89,
      "literature_source": "Merrill et al. (1969)",
      "details": "Pressure drop validation"
    }
  ]
}
```

## Continuous Integration

Validation tests are run automatically on:
- Every commit to main branch
- Weekly scheduled runs
- Before each release

See `.github/workflows/validation.yml` for CI configuration.

## Contributing

To add new validation cases:

1. Create a new example in `examples/`
2. Add Python comparison script in `validation/`
3. Document literature reference
4. Set appropriate tolerance based on expected accuracy
5. Update this README

## Contact

For questions about validation methodology:
- Open an issue on GitHub
- Email: ryan@clanton.dev
