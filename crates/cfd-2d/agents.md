# cfd-2d — Agent Reference

> **Role**: Full 2D incompressible Navier–Stokes solver suite.  
> **Depends on**: `cfd-core`, `cfd-math`, `cfd-schematics`, `cfd-io`

---

## Purpose

`cfd-2d` implements steady and unsteady incompressible Navier–Stokes in 2D on
structured and unstructured grids, with a full turbulence model zoo and multiple
numerical methods.

### Governing Equations

Incompressible Navier–Stokes:
```
∂u/∂t + (u · ∇)u = −∇p/ρ + ν∇²u + f     (momentum)
∇ · u = 0                                   (continuity)
```

Non-dimensionalised form:
```
Re = ρUL/μ
```

---

## Module Structure

```
src/
  lib.rs                    pub mods; solver table; governing equation docs

  grid/
    mod.rs
    structured.rs           StructuredGrid2D: uniform / stretched Cartesian
    unstructured.rs         UnstructuredGrid2D: triangles / quads
    boundaries.rs           GridBoundary: patch labelling
    refinement.rs           Local refinement: h-refinement of structured quadrilaterals

  solvers/
    mod.rs
    fdm/                    Finite Difference Method
      mod.rs
      advection_diffusion.rs  Central / upwind FD for scalar transport
      poisson.rs              FD Poisson solver (pressure)
      mms_tests.rs            MMS verification (order of accuracy checks)
    fvm/                    Finite Volume Method
      mod.rs
      flux.rs                 Face-flux interpolation (linear, TVD)
      geometry.rs             Cell-centred FVM geometry (face normals, areas)
    lbm/                    Lattice Boltzmann Method (D2Q9)
      mod.rs
      d2q9.rs                 D2Q9 weight/lattice vectors
      bgk.rs                  BGK collision (single relaxation time)
      mrt.rs                  MRT collision (multiple relaxation time, Lallemand-Luo)
      regularized.rs          Regularized LBM (Latt-Chopard)
      streaming.rs            Streaming step + BC application
      equilibrium.rs          Maxwell-Boltzmann equilibrium distribution
      validation.rs           Poiseuille / lid-driven cavity reference
    simple.rs               SIMPLE algorithm (Patankar 1980)
    simplec_pimple.rs       SIMPLEC (Van Doormal & Raithby) + PIMPLE loop
    cavity.rs               Lid-driven cavity solver (Ghia benchmark support)
    bifurcation.rs          Bifurcation channel solver (Y-junction 2D)
    serpentine.rs           Serpentine channel solver (2D bends + Dean vortex)
    venturi.rs              Venturi solvers (converging-diverging, cavitation inception)
    poiseuille.rs           Poiseuille entry + fully-developed validation
    non_newtonian_poiseuille.rs  Power-law / Casson Poiseuille

  piso_algorithm.rs         PISO algorithm (Issa 1986): predictor / corrector loops

  physics/
    mod.rs
    turbulence/
      mod.rs
      k_epsilon.rs          Standard k−ε (Launder-Sharma)
      k_omega_sst.rs        k−ω SST (Menter 1994): blended near/far wall
      les_smagorinsky.rs    LES: classic Smagorinsky CS=0.17
      les_dynamic.rs        Dynamic LES: Germano–Lilly procedure
      les_wale.rs           WALE model (near-wall Cs → 0 automatically)
      les_vreman.rs         Vreman algebraic SGS
      les_sigma.rs          Sigma model (Nicoud)
      les_miles.rs          MILES: Monotone Integrated LES
      des.rs                DES97 / DDES / IDDES (Spalart 1997)
      reynolds_stress.rs    Full RSM: 6-equation Reynolds-stress transport
      spalart_allmaras.rs   One-equation S-A (SA-neg variant)
      wall_functions.rs     Standard / LRN wall functions; y+ evaluation
      mod.rs
    momentum/
      mod.rs
      coefficients.rs       SIMPLE/PISO face coefficients (aₑ, aₙ, …)
      corrections.rs        QUICK / TVD coefficient corrections
      muscl.rs              MUSCL slope limiters (van Leer, minmod, superbee)
    energy/                 Energy equation (buoyancy-driven) — optional
    vorticity/              Vorticity-streamfunction ψ–ω formulation

  discretization/
    mod.rs
    upwind.rs               1st/2nd-order upwind
    central.rs              Central differencing (2nd order, unbounded)
    tvd.rs                  TVD schemes: van Leer, minmod, superbee, UMIST
    weno.rs                 WENO5 on 2D structured grids
    stencil.rs              Extended stencil utilities (non-orthogonal correction)
    high_order.rs           4th/6th-order compact FD stencils

  pressure_velocity/
    mod.rs
    rhie_chow.rs            Rhie-Chow interpolation (pressure-velocity decoupling fix)
    simple_coefficients.rs  aP, aₑ, aₙ, … coefficient assembly

  schemes/
    time/
      explicit.rs           Forward Euler, Adams-Bashforth
      implicit.rs           Backward Euler, Crank-Nicolson
      multistep.rs          BDF2, BDF3
      adaptive.rs           Adaptive Δt from CFL/diffusive limit
```

---

## Solver Reference

| Solver | Module | Algorithm | Notes |
|--------|--------|-----------|-------|
| SIMPLE | `solvers/simple.rs` | Semi-Implicit linked equations | Patankar 1980; steady-state |
| SIMPLEC | `solvers/simplec_pimple.rs` | SIMPLEC variant | Faster convergence via improved p' |
| PISO | `piso_algorithm.rs` | Predictor + 2 correctors | Transient; Issa 1986 |
| PIMPLE | `solvers/simplec_pimple.rs` | PISO + outer SIMPLE loop | Transient with large Δt |
| FDM | `solvers/fdm/` | Central/upwind FD | Structured grids |
| FVM | `solvers/fvm/` | Cell-centred FVM | Unstructured tri/quad |
| LBM D2Q9 | `solvers/lbm/` | BGK / MRT / Regularized | Explicit; Re < 1000 |

---

## Turbulence Model Selection Guide

| Re range | Geometry type | Recommended model |
|---------|--------------|-------------------|
| < 2000 | Any (millifluidic) | Laminar N-S (no turbulence) |
| 2000–10 000 | Channel / pipe | k-ω SST |
| > 10 000 | Bluff body wake | LES Dynamic or DDES |
| Separated flow | Venturi / step | k-ω SST + LRN wall function |
| DNS cross-check | Validation only | — (no turbulence model) |

> **Note**: Most millifluidic SDT channels operate at Re ≪ 100. Turbulence models
> are present for completeness and for the venturi cavitation case, not routine use.

---

## Key Algorithms

### Rhie-Chow Interpolation (`pressure_velocity/rhie_chow.rs`)

Prevents checkerboard pressure oscillations on co-located grids:
```
uₑ = (ū)ₑ − Δt·(dp/dx)ₑ + Δt·(dp/dx)̄ₑ
```
where `(·)ₑ` denotes face value and overbar denotes cell-averaged interpolation.

### SIMPLE Algorithm (`solvers/simple.rs`)

1. Guess `p*`, solve momentum → `u*`  
2. Solve pressure-correction `p'` from continuity violation  
3. Correct: `p = p* + αₚp'`, `u = u* + correction(∇p')`  
4. Repeat until `‖residual‖ < ε`

Relaxation: `αᵤ = 0.7`, `αₚ = 0.3` defaults.

### LBM BGK (`solvers/lbm/bgk.rs`)

**Theorem (Chapman-Enskog expansion)**: The LBM D2Q9-BGK recovers the
incompressible N-S equations to O(Ma²) when `τ = ν/(c_s²·Δt) + 0.5`:
```
fᵢ* = fᵢ − (fᵢ − fᵢ^eq) / τ    (collision)
fᵢ(x + cᵢΔt, t + Δt) = fᵢ*(x, t)  (streaming)
```
MRT variant uses transformation matrix M for improved stability at `Re > 100`.

---

## Key Mathematical Theorems

| Module | Theorem |
|--------|---------|
| `solvers/lbm/bgk.rs` | Chapman-Enskog: BGK → N-S at O(Ma²) |
| `piso_algorithm.rs` | PISO: second-order time accuracy with prediction-correction |
| `physics/turbulence/k_omega_sst.rs` | Menter SST blending: k-ε far, k-ω near wall |
| `discretization/tvd.rs` | TVD condition: ‖Δu^{n+1}‖ ≤ ‖Δu^n‖ (no new extrema) |
| `pressure_velocity/rhie_chow.rs` | Rhie-Chow eliminates pressure-velocity decoupling |

---

## Validation References

| Benchmark | Module | Reference |
|-----------|--------|-----------|
| Ghia lid-driven cavity | `solvers/cavity.rs` | Ghia et al. 1982 |
| Poiseuille (parabolic profile) | `solvers/poiseuille.rs` | Exact: u = Q·y(H-y)/(2νH³/3) |
| LBM Poiseuille | `solvers/lbm/validation.rs` | Chapman-Enskog theory |
| Venturi pressure drop | `solvers/venturi.rs` | Bernoulli + minor losses |

---

## Relationship to Other Crates

| Crate | Relationship |
|-------|-------------|
| `cfd-schematics` | Geometry coordinates seed grid boundaries |
| `cfd-1d` | 1D nodal pressures provide inlet/outlet BCs for 2D channel refinement |
| `cfd-3d` | `cfd-2d` benchmarks validate corresponding 3D slices |
| `cfd-validation` | Ghia cavity, Poiseuille, LBM reference comparisons |
| `cfd-io` | VTK and CSV output for field visualisation |

---

## Prohibited Patterns

- No 3D arrays or 3D index loops — keep 2D
- LBM must not exceed `Ma > 0.1` (compressibility errors) — enforce in inlet BC
- Turbulence models must not be applied at Re < 1000 without explicit justification
- Rhie-Chow is mandatory for co-located SIMPLE/PISO — do not remove it
