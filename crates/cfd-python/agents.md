# cfd-python — Agent Reference

> **Role**: Python bindings for CFDrs via PyO3. Exposes Rust solvers as a `pycfdrs` Python wheel.  
> **Build tool**: `maturin develop` (`.venv`) / `maturin build --release`  
> **Depends on**: all solver crates (`cfd-1d`, `cfd-2d`, `cfd-3d`, `cfd-core`, `cfd-validation`)

---

## Purpose

`cfd-python` makes CFDrs accessible from Python for:

- Rapid prototyping and design exploration
- Integration with scientific Python stack (NumPy, SciPy, Matplotlib)
- External validation scripts and Jupyter notebooks
- Platform for SDT device design Python workflows

All computation happens in Rust; Python sees zero-copy results via `pyo3`.

---

## Module Structure

```
src/
  lib.rs                   #[pymodule] pycfdrs: registers all submodules
  bifurcation.rs           PyBifurcationSolver, PyBifurcationSolver2D, PyBifurcation3DSolver
  blood.rs                 PyCassonBlood, PyCarreauYasudaBlood, PyCrossBlood, PyFahraeusLindqvist
  poiseuille_2d.rs         PyPoiseuilleSolver (2D), PyPoiseuille2DSolver
  result_types.rs          PyFlowResult, PyPressureResult, PyValidationResult (shared result types)
  solver_2d.rs             PyVenturiSolver2D, PyTrifurcationSolver2D, PyCavitySolver2D, PyBifurcationSolver2D
  solver_3d.rs             PyBifurcation3DSolver, PyTrifurcation3DSolver, PyPoiseuille3DSolver,
                           PyVenturi3DSolver, PySerpentine3DSolver
  womersley.rs             PyWomersleyNumber, PyWomersleyProfile, PyWomersleyFlow
```

Additional bindings in `lib.rs`:
- `PySerpentineSolver1D` — 1D serpentine channel
- `PyVenturiSolver1D` — 1D venturi
- `PyTrifurcationSolver` — 1D trifurcation
- Validation submodule re-exports

---

## Exposed Python Classes

### 1D Solvers

| Python class | Rust solver | Module |
|-------------|------------|--------|
| `BifurcationSolver` | `cfd_1d::bifurcation` | `bifurcation.rs` |
| `TrifurcationSolver` | `cfd_1d::trifurcation` | `lib.rs` |
| `SerpentineSolver1D` | `cfd_1d::serpentine` | `lib.rs` |
| `VenturiSolver1D` | `cfd_1d::venturi` | `lib.rs` |

### 2D Solvers

| Python class | Rust solver | Module |
|-------------|------------|--------|
| `PoiseuilleSolver` | `cfd_2d::poiseuille` | `poiseuille_2d.rs` |
| `VenturiSolver2D` | `cfd_2d::venturi` | `solver_2d.rs` |
| `TrifurcationSolver2D` | `cfd_2d::trifurcation` | `solver_2d.rs` |
| `CavitySolver2D` | `cfd_2d::cavity` (Ghia) | `solver_2d.rs` |
| `BifurcationSolver2D` | `cfd_2d::bifurcation` | `solver_2d.rs` |
| `Poiseuille2DSolver` | `cfd_2d::poiseuille` (FVM variant) | `poiseuille_2d.rs` |

### 3D Solvers

| Python class | Rust solver | Module |
|-------------|------------|--------|
| `Bifurcation3DSolver` | `cfd_3d::bifurcation` | `solver_3d.rs` |
| `Trifurcation3DSolver` | `cfd_3d::trifurcation` | `solver_3d.rs` |
| `Poiseuille3DSolver` | `cfd_3d::poiseuille` | `solver_3d.rs` |
| `Venturi3DSolver` | `cfd_3d::venturi` | `solver_3d.rs` |
| `Serpentine3DSolver` | `cfd_3d::serpentine` | `solver_3d.rs` |

### Blood Models

| Python class | Model | Module |
|-------------|-------|--------|
| `CassonBlood` | Casson: √τ = √τ_y + √(μ_∞·γ̇) | `blood.rs` |
| `CarreauYasudaBlood` | Carreau-Yasuda | `blood.rs` |
| `CrossBlood` | Cross power-law | `blood.rs` |
| `FahraeusLindqvist` | Fahraeus-Lindqvist (diameter-dependent viscosity) | `blood.rs` |

### Womersley Flow

| Python class | Purpose | Module |
|-------------|---------|--------|
| `WomersleyNumber` | α = R√(ω/ν) | `womersley.rs` |
| `WomersleyProfile` | u(r,t) via Bessel functions | `womersley.rs` |
| `WomersleyFlow` | Complete pulsatile flow cycle | `womersley.rs` |

---

## Example Usage

```python
import pycfdrs

# 1D bifurcation
solver = pycfdrs.BifurcationSolver(
    inlet_pressure_pa=200.0,
    outlet_pressure_pa=0.0,
    parent_diameter_m=200e-6,
    daughter_diameter_m=140e-6,
)
result = solver.solve()
print(f"Flow rate: {result.flow_rate_m3s * 1e9 * 60:.3f} µL/min")
print(f"Pressure drop: {result.pressure_drop_pa:.1f} Pa")

# Casson blood model
blood = pycfdrs.CassonBlood(yield_stress_pa=0.0015, viscosity_pa_s=0.0035)
viscosity = blood.effective_viscosity(shear_rate_s1=100.0)

# Womersley pulsatile profile
womersley = pycfdrs.WomersleyFlow(radius_m=100e-6, frequency_hz=1.2, nu=3.5e-6)
profile = womersley.velocity_profile(r_points=64, t=0.25)
```

---

## Build Instructions

```sh
# Development install (editable wheel in .venv)
cd crates/cfd-python
maturin develop

# Release wheel (for distribution)
maturin build --release

# Run Python validation suite
python -c "import pycfdrs; print(pycfdrs.__version__)"
python tests/validate_pycfdrs.py
```

---

## Result Types (`result_types.rs`)

All solver `solve()` methods return `PyFlowResult`:

```python
result.flow_rate_m3s        # float — volume flow rate (m³/s)
result.pressure_drop_pa     # float — inlet-outlet ΔP (Pa)
result.reynolds_number      # float — characteristic Re
result.wall_shear_pa        # float — max wall shear stress (Pa)
result.converged            # bool — solver convergence flag
result.iterations           # int — iteration count
```

Validation methods return `PyValidationResult`:

```python
result.l2_error             # float — L₂ norm vs. reference
result.linf_error           # float — L∞ norm
result.observed_order       # float — Richardson extrapolation order
result.pass                 # bool — passes tolerance check
```

---

## Binding Conventions

| Convention | Rationale |
|------------|-----------|
| `#[pyclass(name = "NameWithout Py")]` | Python sees `BifurcationSolver`, not `PyBifurcationSolver` |
| All methods return `PyResult<T>` | Python exceptions map to Rust `OptimError` / `SolverError` |
| No Python sate mutation via `&mut self` — return new objects | Immutability; thread-safe from Python |
| NumPy arrays via `pyo3-numpy` | Zero-copy field arrays for profile outputs |
| SI units everywhere | No unit duality at the binding layer |

---

## Prohibited Patterns

- No solver logic in `cfd-python` source — all physics must live in Rust crates
- Never expose internal Rust types (`SlotMap` keys, mesh handles) to Python — create wrapper value types
- Python class names must match their conceptual Rust solver without the `Py` prefix
- Do not add Python-only convenience methods — keep the binding thin; add logic to the Rust crate
