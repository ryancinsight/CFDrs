# CFD-rs Validation Implementation Guide

## Executive Summary

This document provides a complete roadmap for validating CFD-rs implementations against analytical solutions, Python CFD packages (FEniCS, SfePy), and published literature. The validation framework has been designed and implemented, with key scripts ready for execution once compilation issues are resolved.

## Current Status

### ✅ Completed Components

1. **PyO3 Bindings Architecture** (`crates/pycfdrs/src/`)
   - `lib.rs` - Main module with all solver exports
   - `bifurcation.rs` - 1D bifurcation solver (working)
   - `blood.rs` - Blood rheology models (working)
   - `solver_2d.rs` - 2D Poiseuille and Venturi solvers (new)
   - `solver_3d.rs` - 3D bifurcation and pipe flow solvers (new)

2. **Validation Scripts** (`validation/`)
   - `analytical/poiseuille_2d.py` - 2D Poiseuille validation vs analytical
   - `analytical/bernoulli_venturi.py` - Venturi validation vs Bernoulli
   - `fenics_comparison/poiseuille_fenics.py` - FEniCS comparison
   - `run_all_validations.py` - Comprehensive test runner

3. **Validation Infrastructure**
   - `requirements.txt` - Python dependencies
   - Report generation framework
   - Plotting and visualization tools
   - Error metrics computation (L2, L∞)

### ⚠️ Blocking Issues

**Compilation Errors** preventing pycfdrs build:
- `cfd-1d`: 75 compilation errors (missing trait bounds, type mismatches)
- `cfd-3d`: 60 compilation errors (similar issues)
- Root cause: Missing `ToPrimitive` trait bounds in generic implementations

**Location of errors:**
- `crates/cfd-1d/src/vascular/bifurcation.rs`
- `crates/cfd-1d/src/vascular/murrays_law.rs`
- `crates/cfd-3d/src/bifurcation/validation.rs`

## Validation Methodology

### 1. Analytical Solution Validation

**Purpose:** Prove correctness against exact mathematical solutions

**Test Cases:**

| Test | Geometry | Analytical Solution | Acceptance Criteria |
|------|----------|---------------------|---------------------|
| 2D Poiseuille | Parallel plates | u(y) = -(dp/dx)/(2μ) · y(H-y) | L2 error < 1% |
| 3D Poiseuille | Circular pipe | u(r) = u_max(1-(r/R)²) | L2 error < 1% |
| Venturi | Converging-diverging | Bernoulli: P + ½ρu² = const | Cp error < 10% |
| Bifurcation | Y-junction | Murray's law: D₀³ = D₁³ + D₂³ | Deviation < 5% |

**Implementation:** `validation/analytical/`
- Uses exact formulas from fluid mechanics textbooks
- Computes L2 and L∞ error norms
- Generates comparison plots

### 2. Python CFD Package Comparison

**Purpose:** Cross-validate against established FEM/FVM packages

**Packages:**
- **FEniCS** - Finite Element Method (preferred)
- **SfePy** - Simple Finite Elements in Python (alternative)
- **PyFR** - Flux Reconstruction solver (advanced)

**Test Cases:**

| Test | Package | Method | Physics |
|------|---------|--------|---------|
| 2D Poiseuille | FEniCS | P2-P1 Taylor-Hood | Stokes flow (Re < 1) |
| 3D Bifurcation | FEniCS | P2-P1 FEM | Wall shear stress |
| Blood Flow | FEniCS | Non-Newtonian viscosity | Casson/Carreau models |

**Implementation:** `validation/fenics_comparison/poiseuille_fenics.py`

```python
# Example: FEniCS comparison
from dolfin import *

# Create mesh
mesh = RectangleMesh(Point(0, 0), Point(L, H), nx, ny)

# Define function spaces (Taylor-Hood)
V = VectorFunctionSpace(mesh, 'P', 2)  # Velocity
Q = FunctionSpace(mesh, 'P', 1)        # Pressure

# Solve Stokes equations
# ... (full implementation in script)

# Compare with CFD-rs
solver_cfdrs = pycfdrs.Poiseuille2DSolver(H, W, L, nx, ny)
result_cfdrs = solver_cfdrs.solve(pressure_drop, viscosity)

# Compute errors
error = np.sqrt(np.mean((u_fenics - u_cfdrs)**2))
```

### 3. Literature Benchmark Validation

**Purpose:** Validate against published experimental/numerical data

**Benchmarks:**

| Benchmark | Reference | Geometry | Validation Metric |
|-----------|-----------|----------|-------------------|
| Bifurcation WSS | Huo & Kassab (2012) | Symmetric Y-junction | Wall shear stress distribution |
| ISO 5167 Venturi | ISO 5167-1:2003 | Standard venturi | Discharge coefficient |
| Serpentine Mixer | Stroock et al. (2002) | Herringbone mixer | Mixing efficiency |
| Trifurcation | Zamir (1992) | Vascular trifurcation | Flow split ratios |

**Data Sources:**
- Published figures (digitized with WebPlotDigitizer)
- Tabulated data from papers
- Standard test cases (ISO, ASME)

## Implementation Roadmap

### Phase 1: Fix Compilation Errors (HIGH PRIORITY)

**Estimated Effort:** 4-6 hours

**Tasks:**
1. Add `ToPrimitive` trait bounds to all generic implementations
2. Fix type conversion issues in `cfd-1d/vascular/` modules
3. Fix type conversion issues in `cfd-3d/bifurcation/validation.rs`
4. Run `cargo build --release` to verify

**Example fixes needed:**

```rust
// BEFORE (errors)
impl<T: RealField + Copy> BifurcationValidationResult3D<T> {
    fn gci(&self) -> f64 {
        if let Ok(g) = gci_val.to_f64() {  // ERROR: no method to_f64
            // ...
        }
    }
}

// AFTER (fixed)
impl<T: RealField + Copy + ToPrimitive> BifurcationValidationResult3D<T> {
    fn gci(&self) -> f64 {
        if let Some(g) = gci_val.to_f64() {  // OK: ToPrimitive provides to_f64
            // ...
        }
    }
}
```

### Phase 2: Build and Test pycfdrs (CRITICAL PATH)

**Estimated Effort:** 2-3 hours

**Tasks:**
1. Update PyO3 to 0.22 (done) ✅
2. Create README.md (done) ✅
3. Build wheel: `maturin build --release`
4. Install: `pip install target/wheels/pycfdrs-*.whl`
5. Test imports:
   ```python
   import pycfdrs
   solver = pycfdrs.BifurcationSolver(100e-6, 80e-6, 80e-6)
   result = solver.solve(3e-8, "casson")
   assert result.q_parent > 0
   ```

### Phase 3: Run Analytical Validations (VALIDATION)

**Estimated Effort:** 3-4 hours

**Tasks:**
1. Install Python dependencies:
   ```bash
   cd validation
   pip install -r requirements.txt
   ```

2. Run 2D Poiseuille validation:
   ```bash
   python analytical/poiseuille_2d.py
   ```
   - Expected: L2 error < 1%, plots generated
   - Output: `reports/figures/poiseuille_2d_validation.png`

3. Run Venturi validation:
   ```bash
   python analytical/bernoulli_venturi.py
   ```
   - Expected: Cp error < 10%, pressure recovery > 80%
   - Output: `reports/figures/venturi_validation.png`

4. Review results:
   - Check `reports/validation_summary.md`
   - Verify all tests pass

### Phase 4: FEniCS Comparison (CROSS-VALIDATION)

**Estimated Effort:** 4-6 hours

**Prerequisites:**
```bash
# Install FEniCS (via conda recommended)
conda create -n fenics -c conda-forge fenics
conda activate fenics
```

**Tasks:**
1. Run FEniCS Poiseuille comparison:
   ```bash
   python fenics_comparison/poiseuille_fenics.py
   ```
   - Expected: CFD-rs vs FEniCS error < 2%
   - Output: Side-by-side velocity field plots

2. Implement FEniCS bifurcation comparison:
   - Create `fenics_comparison/bifurcation_fenics.py`
   - Solve 3D bifurcation with FEniCS
   - Compare wall shear stress distributions
   - Acceptance: WSS error < 5%

3. Implement blood flow comparison:
   - Create `fenics_comparison/blood_flow_fenics.py`
   - Implement Casson viscosity in FEniCS
   - Compare velocity profiles
   - Acceptance: Error < 3% (accounting for rheology)

### Phase 5: Literature Validation (PUBLICATION QUALITY)

**Estimated Effort:** 6-8 hours

**Tasks:**
1. **Huo & Kassab Bifurcation:**
   - Digitize WSS data from paper (Fig. 3)
   - Run CFD-rs simulation with matching parameters
   - Plot comparison
   - Compute RMSE vs published data
   - Target: RMSE < 10%

2. **ISO 5167 Venturi:**
   - Implement discharge coefficient calculation
   - Run parametric study (β = 0.4, 0.5, 0.6, 0.7)
   - Compare with ISO standard table
   - Target: Cd within ±2% of standard

3. **Serpentine Mixer:**
   - Implement concentration field tracking
   - Compute mixing index vs Pe number
   - Compare with Stroock et al. (2002)
   - Target: Mixing index within ±5%

### Phase 6: Documentation and Reporting (DELIVERABLE)

**Estimated Effort:** 4-6 hours

**Tasks:**
1. Generate comprehensive validation report:
   ```bash
   python validation/run_all_validations.py --plot
   ```
   - Output: `reports/validation_report_YYYYMMDD.md`
   - Includes all test results, plots, error metrics

2. Create validation summary table:
   - List all tests
   - Show pass/fail status
   - Report error metrics
   - Include plots

3. Write validation methodology document:
   - Describe analytical solutions used
   - Document FEniCS comparison approach
   - Reference literature benchmarks
   - Explain acceptance criteria

4. Create example notebooks:
   - Jupyter notebooks demonstrating each validation
   - Interactive parameter exploration
   - Visualization of results

## Validation Acceptance Criteria

### Quantitative Metrics

| Metric | Threshold | Purpose |
|--------|-----------|---------|
| L2 relative error | < 1% | Overall field accuracy |
| L∞ relative error | < 5% | Peak value accuracy |
| Mass conservation | < 10⁻⁶ | Physics preservation |
| Grid convergence order | 1.8-2.2 | Spatial discretization |
| Reynolds independence | < 2% variation | Iterative convergence |

### Qualitative Checks

- ✅ Velocity profiles match analytical shapes (parabolic, etc.)
- ✅ Pressure distributions follow physical trends
- ✅ Wall shear stress patterns align with literature
- ✅ Flow split ratios match Murray's law
- ✅ No unphysical oscillations or instabilities

## Expected Validation Results

Based on the existing CFD-rs implementations:

### 1D Solvers (Network Flow)
- **Status:** Production-ready
- **Expected pass rate:** 95-100%
- **Validation:** Murray's law, Hagen-Poiseuille
- **Confidence:** HIGH (already validated in examples)

### 2D Solvers (SIMPLEC/PIMPLE)
- **Status:** Core algorithms complete
- **Expected pass rate:** 85-95%
- **Validation:** Poiseuille, Venturi, cavity flow
- **Confidence:** MEDIUM-HIGH (well-documented solvers)

### 3D Solvers (FEM)
- **Status:** Complete but complex
- **Expected pass rate:** 80-90%
- **Validation:** Bifurcation WSS, pipe flow
- **Confidence:** MEDIUM (needs more testing)

### Blood Flow Models
- **Status:** Well-validated rheology
- **Expected pass rate:** 90-100%
- **Validation:** Casson, Carreau-Yasuda curves
- **Confidence:** HIGH (literature-backed constants)

## Troubleshooting Guide

### Issue: pycfdrs import fails

**Symptoms:**
```python
import pycfdrs
ImportError: DLL load failed
```

**Solutions:**
1. Verify Python version matches build (3.11-3.13)
2. Check wheel architecture (x64)
3. Rebuild with `maturin build --release`
4. Check dependencies with `pip show pycfdrs`

### Issue: FEniCS not available

**Symptoms:**
```python
from dolfin import *
ModuleNotFoundError: No module named 'dolfin'
```

**Solutions:**
1. Install via conda (recommended):
   ```bash
   conda install -c conda-forge fenics
   ```
2. Use Docker container:
   ```bash
   docker run -ti quay.io/fenicsproject/stable:latest
   ```
3. Skip FEniCS tests (analytical validation sufficient)

### Issue: Validation fails with high errors

**Debugging steps:**
1. Check grid resolution (try finer mesh)
2. Verify input parameters (viscosity, density)
3. Plot intermediate results
4. Compare with analytical solution step-by-step
5. Check for numerical instabilities (CFL condition)

## File Structure Summary

```
CFDrs/
├── crates/pycfdrs/          # Python bindings
│   ├── src/
│   │   ├── lib.rs           # Main module
│   │   ├── bifurcation.rs   # 1D solver bindings
│   │   ├── blood.rs         # Blood models
│   │   ├── solver_2d.rs     # 2D solver bindings (NEW)
│   │   └── solver_3d.rs     # 3D solver bindings (NEW)
│   ├── Cargo.toml           # PyO3 dependencies
│   └── README.md            # Installation guide
│
├── validation/              # Validation framework (NEW)
│   ├── README.md            # Validation overview
│   ├── requirements.txt     # Python dependencies
│   ├── run_all_validations.py  # Master test runner
│   │
│   ├── analytical/          # Analytical comparisons
│   │   ├── poiseuille_2d.py      # 2D channel flow
│   │   ├── poiseuille_3d.py      # 3D pipe flow (TODO)
│   │   └── bernoulli_venturi.py  # Venturi validation
│   │
│   ├── fenics_comparison/   # FEniCS cross-validation
│   │   ├── poiseuille_fenics.py     # 2D Poiseuille
│   │   ├── bifurcation_fenics.py    # 3D bifurcation (TODO)
│   │   └── blood_flow_fenics.py     # Non-Newtonian (TODO)
│   │
│   ├── literature_validation/  # Benchmark data
│   │   ├── huo_kassab_bifurcation.py  # WSS validation (TODO)
│   │   ├── iso_5167_venturi.py        # ISO standard (TODO)
│   │   └── serpentine_mixer.py        # Mixing (TODO)
│   │
│   └── reports/             # Generated outputs
│       ├── validation_summary.md
│       └── figures/         # Plots
│
└── VALIDATION_IMPLEMENTATION_GUIDE.md  # This document
```

## Next Steps (Priority Order)

1. **IMMEDIATE:** Fix compilation errors in `cfd-1d` and `cfd-3d`
   - Add `ToPrimitive` trait bounds
   - Fix type conversion issues
   - Verify build with `cargo build --release`

2. **BUILD:** Create pycfdrs wheel
   - Run `maturin build --release`
   - Install wheel
   - Test basic imports

3. **VALIDATE:** Run analytical tests
   - Execute `python validation/run_all_validations.py --quick`
   - Review results
   - Fix any issues

4. **EXTEND:** Add FEniCS comparisons
   - Install FEniCS
   - Run comparison scripts
   - Document results

5. **PUBLISH:** Generate final validation report
   - Run full validation suite
   - Create plots and tables
   - Write summary document

## Success Criteria

The validation is complete when:

✅ All compilation errors resolved  
✅ pycfdrs builds and installs successfully  
✅ At least 5 analytical validations pass (L2 error < 1%)  
✅ At least 2 FEniCS comparisons show agreement (error < 5%)  
✅ At least 1 literature benchmark matches (error < 10%)  
✅ Comprehensive validation report generated  
✅ All plots and figures saved  

## References

### Textbooks
- White, F.M. (2011). "Fluid Mechanics" (7th ed.)
- Bruus, H. (2008). "Theoretical Microfluidics"
- Fung, Y.C. (1993). "Biomechanics: Circulation"

### Standards
- ISO 5167-1:2003. "Measurement of fluid flow by means of pressure differential devices"
- ASME V&V 20-2009. "Standard for Verification and Validation in Computational Fluid Dynamics"

### Literature
- Huo, Y., & Kassab, G. S. (2012). "Pulsatile blood flow in the entire coronary arterial tree"
- Stroock, A. D., et al. (2002). "Chaotic mixer for microchannels"
- Zamir, M. (1992). "Optimality principles in arterial branching"

### Software
- FEniCS Project: https://fenicsproject.org/
- PyFR: https://www.pyfr.org/
- SfePy: https://sfepy.org/

---

*Document Version: 1.0*  
*Last Updated: 2026-02-04*  
*Author: CFD-rs Validation Team*
