# cfd-3d — Agent Reference

> **Role**: High-fidelity 3D CFD — FEM, IBM, Level-Set, VOF multiphase, spectral, and domain-specific solvers.  
> **Depends on**: `cfd-core`, `cfd-math`, `cfd-2d` (turbulence models shared), `cfd-io`

---

## Purpose

`cfd-3d` extends the 2D capability to full three-dimensional simulations:

- **Finite Element Method (FEM)**: Galerkin + SUPG/PSPG stabilised 3D velocity-pressure
- **Immersed Boundary Method (IBM)**: complex solid boundaries on Cartesian grids
- **Level-Set**: sharp interface tracking with Sussman redistancing
- **Volume of Fluid (VOF)**: PLIC-reconstructed free-surface and two-phase flow
- **Spectral Methods**: Fourier (FFT-based) and Chebyshev collocation in 3D
- **Domain solvers**: bifurcation, trifurcation, serpentine channel, venturi (each with geometry + solver + validation)

---

## Module Structure

```
src/
  lib.rs                     pub mods; theorem refs; relationship table

  fem/
    mod.rs
    shape_functions.rs        Trilinear Q8 / quadratic Serendipity shape functions N_i(ξ,η,ζ)
    quadrature.rs             3D Gauss-Legendre quadrature rules (1, 2, 3-point per direction)
    galerkin.rs               Galerkin weak form: ∫ N_i · (u·∇u) dΩ + viscous + pressure terms
    stabilization.rs          SUPG / PSPG stabilization parameters (Brooks-Hughes 1982)
    operators.rs              FEM stiffness K, mass M, gradient G, divergence D assembly
    projection_solver.rs      Chorin projection: predictor u* → pressure Poisson → correct u

  ibm/
    mod.rs
    direct_forcing.rs         Exact immersed boundary force: f_ib = (u_desired − u*) / Δt
    feedback_forcing.rs       PID feedback IBM (Goldstein–Handler–Sirovich 1993)
    lagrangian.rs             Lagrangian marker set (solid surface points)
    interpolation.rs          Discrete delta functions: Roma (2-pt), Peskin (4-pt)
    tests.rs

  level_set/
    mod.rs
    sdf.rs                    Signed Distance Function φ(x): φ > 0 fluid, φ < 0 solid
    sussman.rs                Sussman redistancing: ∂φ/∂τ + S(φ)(|∇φ| − 1) = 0
    narrow_band.rs            Narrow-band acceleration: only update ±4Δx band
    curvature.rs              Interface curvature κ = ∇·(∇φ/|∇φ|)

  vof/
    mod.rs
    reconstruction.rs         PLIC: plane nₖ · x = αₖ in each cell
    advection.rs              Geometric VOF advection (unsplit Lagrangian / operator split)
    plic.rs                   Youngs gradient estimate + PLIC plane parameters
    surface_tension.rs        CSF (Brackbill–Kothe–Zemach): f_σ = σ·κ·δ(φ)·n̂
    cavitation_solver.rs      Rayleigh-Plesset coupled VOF cavitation

  spectral/
    mod.rs
    fourier.rs                3D Fourier spectral method (FFT-based)
    chebyshev.rs              3D Chebyshev collocation
    chebyshev_tests.rs        Exponential convergence verification
    poisson_solver.rs         Spectral Poisson: P̂_k = −f̂_k / (k_x² + k_y² + k_z²)
    collocation.rs            Pseudospectral collocation points

  physics/
    turbulence.rs             Re-exports `cfd-2d` LES/DES/RANS models; 3D field bookkeeping

  bifurcation/
    mod.rs
    geometry.rs               3D bifurcation channel geometry (circular tubes, swept)
    solver.rs                 3D N-S solver for Y-bifurcation (SIMPLE 3D)
    validation.rs             Pressure drop vs. 1D H-P; Dean secondary flow

  trifurcation/
    mod.rs
    geometry.rs               3D trifurcation geometry
    solver.rs                 3D trifurcation solver
    validation.rs             Flow uniformity across three outlet branches

  serpentine/
    mod.rs
    geometry.rs               3D serpentine geometry (rectangular channel, N-bends)
    solver.rs                 3D N-S with Dean vortex validation
    validation.rs             Secondary flow profiles vs. Dean number correlation

  venturi/
    mod.rs
    geometry.rs               3D converging-diverging venturi (axisymmetric → 3D extrusion)
    solver.rs                 3D N-S + cavitation detection
    validation.rs             Bernoulli pressure drop + cavitation number
```

---

## Key Algorithms

### FEM: SUPG/PSPG Stabilisation (`fem/stabilization.rs`)

**Theorem (Lax-Milgram / Galerkin stability)**: The Galerkin weak form
```
a(u, v) = l(v)   ∀v ∈ V
```
has a unique solution when `a(·,·)` is coercive and continuous (Lax-Milgram lemma).

Without stabilisation, equal-order velocity-pressure elements violate the
**inf-sup (LBB) condition**. SUPG/PSPG adds residual-based terms:
```
τ_SUPG = h_e / (2‖u‖) · (4ν/h_e² + 2‖u‖/h_e)^{-1}
τ_PSPG = τ_SUPG
```
Stabilisation restores well-posedness for Q1/P1 and Q2/P1 elements.

### IBM: Direct Forcing (`ibm/direct_forcing.rs`)

Body force that enforces no-slip on Lagrangian surface Γ:
```
f_ib(x) = Σ_l (u_desired(X_l) − ũ(X_l)) · δ_h(x − X_l) · ΔV_l
```
where `δ_h` is the discrete delta function and `ũ` is the predicted velocity.

### Level-Set: Sussman Redistancing (`level_set/sussman.rs`)

**Theorem**: The reinitialization PDE
```
∂φ/∂τ + S(φ)(|∇φ| − 1) = 0,   S(φ) = φ / √(φ² + ε²)
```
converges to the signed distance function in pseudo-time τ.

### VOF: PLIC Reconstruction (`vof/plic.rs`)

**Theorem (Volume Conservation)**: PLIC satisfies
```
F^{n+1} = F^n − Δt · Σ_faces F_flux
```
Volume fraction `F` is preserved to machine precision (geometric flux calculation).

Surface tension via CSF: `f_σ = σ · κ · ∇F / ρ̄` — requires smooth curvature estimate.

### Spectral: Fourier Poisson (`spectral/poisson_solver.rs`)

**Theorem (Spectral Poisson)**: In Fourier space,
```
∇²p = f  ⟺  P̂_k = −f̂_k / (|k|²)
```
Solved in O(N log N) via FFT. Valid for periodic domains.

---

## Key Mathematical Theorems

| Module | Theorem |
|--------|---------|
| `fem/galerkin.rs` | Lax-Milgram: unique solution for coercive bilinear form |
| `fem/stabilization.rs` | SUPG/PSPG restores LBB stability for equal-order elements |
| `ibm/direct_forcing.rs` | Direct forcing enforces no-slip to O(Δt) accuracy |
| `level_set/sussman.rs` | Sussman redistancing converges to exact SDF |
| `vof/advection.rs` | PLIC geometric advection conserves volume to machine precision |
| `spectral/poisson_solver.rs` | Spectral Poisson: exact diagonalisation in Fourier space |
| `spectral/chebyshev.rs` | Chebyshev exponential convergence O(e^{-cN}) for analytic fields |

---

## Domain Solver Validation

| Solver | Validation Target | Module |
|--------|-------------------|--------|
| `bifurcation` | 1D H-P pressure drop; Dean secondary flow Re > 100 | `bifurcation/validation.rs` |
| `trifurcation` | Flow rate uniformity < 1% deviation | `trifurcation/validation.rs` |
| `serpentine` | Dean number secondary flow; pressure drop vs. Ito correlation | `serpentine/validation.rs` |
| `venturi` | Bernoulli ΔP; cavitation inception number | `venturi/validation.rs` |

---

## Relationship to Other Crates

| Crate | Relationship |
|-------|-------------|
| `cfd-2d` | Turbulence models shared (LES/DES/RANS); 2D benchmarks cross-validate 3D mid-planes |
| `cfd-mesh` | `HalfEdgeMesh` / `IndexedMesh` provides unstructured mesh input for FEM |
| `cfd-core` | Fluid properties, BCs, `Problem` abstraction |
| `cfd-math` | Sparse linear solvers, preconditioners, spectral operators |
| `cfd-validation` | 3D benchmarks (Ghia 3D, `Bifurcation3DValidation`, `Venturi3DValidation`) |
| `cfd-io` | VTK XML output (`.vts` / `.vtu`) for ParaView visualisation |

---

## Prohibited Patterns

- FEM assembly must not store the full dense stiffness matrix for large problems (use CSR)
- IBM Lagrangian markers must be placed at least 1.5Δx from grid boundaries
- Level-Set redistancing must run every 5–10 time steps, not every step (cost)
- VOF surface tension via CSF requires smoothed curvature; raw `∇F` is too noisy — always smooth
- Spectral methods are only valid on periodic or Chebyshev collocation domains; do not apply to non-periodic BCs without Galerkin basis transformation
