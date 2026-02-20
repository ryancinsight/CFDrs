# cfd-validation — Agent Reference

> **Role**: Verification & Validation (V&V) framework per ASME V&V 20-2009.  
> **Depends on**: `cfd-core`, `cfd-math`, `cfd-1d`, `cfd-2d`, `cfd-3d`

---

## Purpose

`cfd-validation` provides systematic evidence that CFDrs solvers are:
- **Verified** — solving the equations correctly (MMS, Richardson extrapolation, GCI)
- **Validated** — solving the correct equations (analytical solutions, benchmarks, literature)

The project priority hierarchy is defined in the workspace `agents.md`. No production
solver should be shipped without a corresponding entry in this crate.

---

## Module Structure

```
src/
  lib.rs                    pub mods; V&V philosophy comment block

  analytical/
    mod.rs
    poiseuille.rs           Exact: u(y) = U_max(1 − (2y/H)²)
    couette.rs              Exact: u(y) = U_wall · y/H
    womersley.rs            Pulsatile: u(r,t) via Womersley functions J₀/J₁
    taylor_green.rs         Decaying vortex: u(x,y,t) = cos(x)sin(y)e^{-2νt}
    blasius.rs              Laminar boundary layer (Blasius self-similar)
    stokes.rs               Stokes flow: u = (∇p/2μ)(R²−r²) and oscillatory Stokes
    utils.rs                L2/Linf error helpers; grid convergence utilities

  manufactured/
    mod.rs
    mms_advection.rs        MMS for advection equation (exact polynomial source)
    mms_diffusion.rs        MMS for diffusion equation
    mms_burgers.rs          MMS for viscous Burgers equation
    mms_ns.rs               MMS for steady incompressible N-S
    mms_turbulent.rs        MMS for k-ε / k-ω equations
    mms_multi_physics.rs    Coupled MMS (fluid-thermal)
    richardson/
      mod.rs
      core.rs               Richardson extrapolation: p = log(e₁/e₂)/log(r)
      analysis.rs           Effective p determination; grid refinement study
      types.rs              GridLevel, RichardsonResult
      validation.rs         GCI computation and safety factor

  benchmarks/
    mod.rs
    ghia_cavity.rs          Ghia et al. (1982) lid-driven cavity: Re = 100,400,1000
    cylinder.rs             2D cylinder Re = 40 (drag / Strouhal) + Re = 100 (vortex shedding)
    step.rs                 Backward-facing step: Armaly et al. (1983) reattachment length
    bifurcation.rs          2D Y-bifurcation pressure drop vs. 1D H-P
    serpentine.rs           2D serpentine Dean vortex + pressure drop
    venturi.rs              2D venturi Bernoulli + cavitation inception
    trifurcation.rs         2D 1-to-3 uniformity
    bifurcation_3d.rs       3D bifurcation vs. 1D and analytic
    trifurcation_3d.rs      3D trifurcation uniformity
    serpentine_3d.rs        3D serpentine Dean secondary flow
    venturi_3d.rs           3D venturi pressure distribution
    mod.rs

  convergence/
    mod.rs
    criteria.rs             ConvergenceCriteria: relative L2 residual < tol
    study.rs                GridConvergenceStudy: multi-level h-refinement pipeline
    richardson.rs           Richardson utility functions
    gci.rs                  Grid Convergence Index: GCI = F_s|ε| / (r^p − 1)
    error.rs                ConvergenceError

  conservation/
    mod.rs
    traits.rs               ConservationCheck trait
    mass.rs                 ∫∫ (u·n) dS = 0 (closed domain) / = source (open)
    momentum.rs             d/dt ∫ ρu dV = Σ forces (integral momentum balance)
    energy.rs               d/dt ∫ E dV = Σ fluxes + work (energy balance)
    angular_momentum.rs     Torque-angular impulse conservation
    vorticity.rs            Vorticity transport: Dω/Dt = (ω·∇)u + ν∇²ω
    history.rs              ConservationHistory: timestep-by-timestep integral records
    report.rs               ConservationReport: summary table + pass/fail
    tolerance.rs            ConservationTolerance: tols per quantity
    mod.rs

  error_metrics/
    mod.rs
    norms.rs                L1, L2, L∞ norms between numerical and reference fields
    statistics.rs           Mean, RMS, max relative error
    normalized.rs           Volume-normalised norms; point-wise error fields
    analysis.rs             ErrorAnalysis: spatial distribution, localisation

  geometry/
    mod.rs
    geometry_2d.rs          Standard validation domains: channel, cavity, axi-venturi
    geometry_3d.rs          3D validation domains

  literature/
    mod.rs
    blood_flow_1d.rs        Womersley pulsatile blood flow (Womersley 1955)
    chapman_enskog.rs       LBM: Chapman-Enskog recovery of N-S (Lallemand-Luo)
    patankar_1980.rs        SIMPLE algorithm reference results (Patankar 1980)

  reporting/
    mod.rs
    html.rs                 Self-contained HTML V&V report
    json.rs                 Machine-readable JSON validation record
    markdown.rs             GitHub-compatible Markdown V&V summary

  time_integration/
    mod.rs
    integrators.rs          Euler, RK2, RK4 reference implementations for comparison
    stability_analysis.rs   von Neumann stable region for each scheme
    validation.rs           Order-of-accuracy check for time integrators

  benchmarking/
    mod.rs
    performance.rs          Solver throughput benchmarks (cells/second)
    scaling.rs              Strong / weak scaling (single-node multi-thread)
    memory_profiling.rs     Peak RSS, allocation counts, cache miss estimation
```

---

## Verification Methods

### Method of Manufactured Solutions (MMS)

**Procedure** for any new solver:

1. Choose an exact solution `u_exact(x,t)` (smooth polynomial or trigonometric)
2. Compute the source forcing `f = Lu_exact` analytically (symbolic or by hand)
3. Solve the modified PDE `Lu = f` numerically
4. Confirm `‖u_numerical − u_exact‖ = O(h^p)` where `p` is the expected order

```rust
use cfd_validation::manufactured::mms_ns::NsMmsCase;

let case = NsMmsCase::new(Re: 100.0, order_expected: 2);
let result = case.run_convergence_study(&solver, &[h/4, h/2, h, 2*h])?;
assert!(result.observed_order >= 1.9);  // ≥ 2nd order confirmed
```

### Richardson Extrapolation (`manufactured/richardson/`)

**Theorem (Richardson 1911)**:
```
u_exact ≈ u_h + e·h^p + O(h^{p+1})

p = log(‖e₁‖/‖e₂‖) / log(r)   where r = h₂/h₁ is the refinement ratio
```

**Grid Convergence Index (GCI)**:
```
GCI = F_s · |ε| / (r^p − 1)    where F_s = 1.25 (ASME safety factor)
```

GCI < 1% implies grid-independent solution. Report GCI for all production results.

---

## Benchmark References

| Benchmark | Expected result | Reference |
|-----------|----------------|-----------|
| Ghia cavity Re=100 | u_centreline matches Ghia Table 1 | Ghia et al. 1982, JCP 48:387 |
| Ghia cavity Re=1000 | u_centreline matches Ghia Table 1 | Ghia et al. 1982, JCP 48:387 |
| Backward-facing step Re=800 | x_reattach / h = 11.9 ± 0.5 | Armaly et al. 1983, JFM 127:473 |
| 2D cylinder Re=40 | Cd = 1.52 ± 0.03 | Tritton 1959 |
| 2D cylinder Re=100 | St = 0.165 ± 0.002 | Williamson 1989 |
| Poiseuille 2D | max error < 0.1% | Exact parabolic |
| Womersley α=2  | amplitude + phase match | Womersley 1955 |

---

## Conservation Checking

Every time-accurate simulation should run `ConservationMonitor`:

```rust
use cfd_validation::conservation::{MassConservation, MomentumConservation};

let mass_ok = MassConservation::check(&flow_field, tol: 1e-10)?;
let mom_ok  = MomentumConservation::check(&flow_field, &forces, tol: 1e-8)?;
```

Report failures immediately — a non-conservative solver is unphysical regardless of residual convergence.

---

## Key Mathematical Theorems

| Module | Theorem |
|--------|---------|
| `manufactured/richardson/core.rs` | Richardson extrapolation: observed order recovery |
| `convergence/gci.rs` | GCI = 1.25|ε|/(r^p − 1) — ASME V&V 20-2009 |
| `analytical/poiseuille.rs` | Exact Poiseuille: u = U_max(1 − (y/H)²) |
| `analytical/womersley.rs` | Womersley oscillatory flow: u(r,t) via Bessel functions |
| `analytical/taylor_green.rs` | Taylor-Green vortex: u(x,y,t) decays as e^{-2νt} |
| `conservation/mass.rs` | Discrete conservation: Σ_faces (F_f · A_f) = 0 |
| `conservation/vorticity.rs` | Kelvin circulation theorem (inviscid limit) |

---

## Reporting

```rust
use cfd_validation::reporting::{MarkdownReporter, HtmlReporter};

let report = ValidationReport::new("Ghia Re=100")
    .add_benchmark("u_centreline", &numerical, &reference, &norms)
    .add_gci_result(&gci)
    .add_conservation(&conservation);

MarkdownReporter::write(&report, "validation/ghia_re100.md")?;
HtmlReporter::write(&report, "validation/ghia_re100.html")?;
```

---

## Prohibited Patterns

- **No placeholder benchmarks** — every benchmark must compare against closed-form or published tabulated data
- **No tolerance cheating** — do not loosen tolerances to pass a test; fix the solver
- GCI must use `F_s = 1.25` (ASME) — do not lower the safety factor
- MMS source terms must be derived analytically, not numerically differenced
- Conservation checks must use surface integrals, not volume-averaged residuals
