# cfd-math — Agent Reference

> **Role**: Numerical-methods library. All PDE solvers depend on this crate.  
> **Depends on**: `cfd-core` only.

---

## Purpose

`cfd-math` provides the complete numerical backbone:

- **Sparse linear algebra** — CSR/CSC assembly + SpMV
- **Iterative linear solvers** — GMRES, BiCGSTAB, Conjugate Gradient
- **Preconditioners** — ILU(0/k), AMG, Jacobi, SOR/SSOR, Schwarz, block
- **High-order methods** — Discontinuous Galerkin (DG), WENO5, spectral (Fourier + Chebyshev)
- **Time integration** — Euler, Runge-Kutta, IMEX, RK-Chebyshev (explicit stability), exponential integrators, adaptive stepping
- **Differentiation** — finite-difference stencils, spectral derivatives, gradient operators
- **Quadrature** — 1D/3D composite rules, tensor-product, Gauss-Legendre, variable quadrature
- **Interpolation** — linear, Lagrange, cubic spline, ENO reconstruction
- **Iterator utilities** — norm iterators, parallel reductions, stencil windows, running statistics
- **SIMD** — AVX2/SSE/NEON/SWAR dispatch (see note below)
- **Performance monitoring** — operation counts, throughput telemetry

---

## Module Structure

```
src/
  lib.rs
  error.rs                    MathError variants
  performance_monitor.rs      OperationCounter, ThroughputMonitor

  sparse/
    mod.rs
    assembly.rs               SparseMatrixAssembler (COO → CSR)
    builder.rs                SparseMatrixBuilder (incremental)
    operations.rs             SpMV, AXPY, dot, norm on CSR matrices
    patterns.rs               Sparsity pattern analysis
    tests.rs

  linear_solver/
    mod.rs
    config.rs                 LinearSolverConfig (tol, max_iter, restart)
    traits.rs                 IterativeLinearSolver, Preconditioner
    gmres/
      mod.rs
      arnoldi.rs              Modified Gram-Schmidt Arnoldi process
      givens.rs               Givens QR for upper Hessenberg system
      solver.rs               GMRES(m) restart loop
    bicgstab/mod.rs           BiCGStab(ℓ)
    conjugate_gradient/mod.rs CG + PCG
    direct_solver.rs          LU factorisation wrapper (small systems)
    block_preconditioner.rs   Block-diagonal preconditioner
    matrix_free/
      mod.rs                  Matrix-free Krylov interface
      parallel_solvers.rs     Rayon-parallel SpMV
    operators/
      gpu.rs                  GPU operator wrappers
      momentum.rs             Discretised momentum operator
      poisson.rs              Poisson operator (pressure)
    preconditioners/
      mod.rs
      basic.rs                Identity, Jacobi, SOR, SSOR
      cholesky.rs             Incomplete Cholesky
      deflation.rs            Deflation preconditioner
      ilu/                    ILU(0), ILU(k), triangular solve, parallel
      multigrid/
        mod.rs
        amg.rs                Ruge-Stüben AMG + aggregation coarsening
        coarsening.rs         Strength-of-connection, independent sets
        cycles.rs             V-cycle, W-cycle, F-cycle
        gmg.rs                Geometric multigrid
        interpolation.rs      Prolongation / restriction operators
        restriction.rs        Full-weighting, injection
        smoothers.rs          Gauss-Seidel, Jacobi, Chebyshev smoother
      schwarz.rs              Additive Schwarz decomposition

  high_order/
    mod.rs
    dg/
      mod.rs
      basis.rs                Legendre polynomial basis
      flux.rs                 Upwind, Lax-Friedrichs, Roe flux
      limiter.rs              slope limiter (minmod, superbee)
      operators.rs            DG mass matrix, stiffness, lift
      solver.rs               DG time-stepping loop
      tests.rs
      validation.rs
    spectral/
      mod.rs
      assembly.rs             Spectral stiffness/mass assembly
      element.rs              Spectral element basics
      operators.rs            D, D², Chebyshev differentiation matrices
    weno.rs                   WENO3/WENO5 reconstruction

  time_stepping/
    mod.rs
    traits.rs                 TimeStepper, AdaptiveStepper
    runge_kutta.rs            Classical RK4, SSP-RK3, Dormand-Prince
    imex.rs                   IMEX-RK (stiff diffusion + explicit convection)
    adaptive.rs               Error-based step control (PI controller)
    exponential.rs            Exponential integrators (Krylov φ-functions)
    rk_chebyshev.rs           Runge-Kutta-Chebyshev (extended stability region)
    stability.rs              CFL / von Neumann stability analysis

  differentiation/
    mod.rs
    finite_difference.rs      FD stencils (1st–4th order)
    gradient.rs               Gradient operators on structured grids
    schemes.rs                Scheme comparison utilities
    tests.rs

  integration/
    mod.rs
    traits.rs                 Integrable, QuadratureRule
    quadrature.rs             Gauss-Legendre 1D (up to P=20)
    quadrature_3d.rs          3D tensor-product Gauss rules
    composite.rs              Composite trapezoid / Simpson
    tensor.rs                 Tensor-product Gauss rules
    variable.rs               Adaptive Gauss-Kronrod
    utils.rs

  interpolation/
    mod.rs
    traits.rs                 Interpolator<T>
    linear.rs                 1D/2D/3D trilinear
    lagrange.rs               Lagrange polynomial on arbitrary nodes
    cubic_spline.rs           Natural cubic spline (tridiagonal solve)

  iterators/
    mod.rs
    norms.rs                  NormIteratorExt: L1, L2, L∞ in one pass
    parallel.rs               Rayon parallel reduction helpers
    statistics.rs             Mean, variance, running stats
    stencils.rs               Sliding stencil window iterator
    windows.rs                Overlapping window slices

  simd/
    mod.rs
    arch_detect.rs            CPUFeatures runtime detection
    vector.rs                 SimdVector<T, N> abstraction
    vectorization.rs          Vectorized loop tiling
    cfd.rs                    CFD-specific SIMD kernels
    ops/
      mod.rs
      traits.rs               SimdOps
      x86.rs                  AVX2 / SSE4.1 f32/f64 kernels
      arm.rs                  NEON kernels
    swar/                     Software-SIMD fallback (f32/f64/integer)
    tests.rs
```

---

## Linear Solver Stack

### Solver Selection Guide

| Problem type | Recommended solver | Preconditioner |
|-------------|-------------------|----------------|
| SPD (Poisson, Laplace) | CG | AMG or ILU(0) |
| Non-symmetric NS pressure-correction | GMRES(30) | ILU(k) or AMG |
| Non-symmetric momentum | BiCGStab | ILU(0) |
| Symmetric indefinite (saddle point) | BiCGStab | Block diagonal |
| Small dense (< 200×200) | Direct (LU) | — |

### AMG (`preconditioners/multigrid/amg.rs`)

**Theorem (Ruge-Stüben Convergence)**: For M-matrices, classical AMG with
Ruge-Stüben coarsening achieves O(1) convergence per V-cycle independent of
problem size.

Coarsening strategies:
- `ClassicalCoarsening` — strength-of-connection + independent sets (Ruge-Stüben)
- `AggregationCoarsening` — aggregation-based (smoother transition for non-M matrices)

Smoothers: Gauss-Seidel (default), Jacobi, SOR, Chebyshev polynomial.

Cycles: V-cycle, W-cycle, F-cycle.

### GMRES (`linear_solver/gmres/`)

**Theorem (GMRES Optimality)**: GMRES minimises `‖b − Axₘ‖₂` over the Krylov
subspace `Kₘ(A, r₀)`. Each restart resets the subspace; restart = 30 is the
workspace default.

The Arnoldi process uses **modified Gram-Schmidt** (`arnoldi.rs`); QR is
updated online via **Givens rotations** (`givens.rs`).

---

## High-Order Methods

### WENO5 (`high_order/weno.rs`)

**Theorem (WENO Accuracy)**: WENO5 achieves 5th-order accuracy in smooth regions
and 3rd-order near discontinuities via adaptive smoothness indicators.

Smoothness indicators (Jiang-Shu):
```
βₖ = Σⱼ ∫ (d^j pₖ/dx^j)² dx    (sum over derivative orders j=1..r-1)
```

### Spectral Methods (`high_order/spectral/`)

**Theorem (Spectral Convergence)**: For `u ∈ Cᵏ`, the truncated Chebyshev
expansion converges at rate `O(N⁻ᵏ)`; for analytic `u`, convergence is
exponential `O(e^{-cN})`.

Chebyshev differentiation matrix `D` satisfies: `Dᵢⱼ = cᵢ/cⱼ · (−1)^{i+j} / (xᵢ − xⱼ)`.

### Discontinuous Galerkin (`high_order/dg/`)

**Theorem (DG Stability)**: With appropriate upwind numerical fluxes and
slope limiting, DG methods are stable for nonlinear conservation laws.

Energy estimate: `d/dt ‖u‖² ≤ C · ‖f‖²` (discrete entropy inequality).

---

## Time Integration

| Method | Module | Order | Stability |
|--------|--------|-------|-----------|
| Explicit Euler | `runge_kutta` | 1 | Conditionally stable (CFL) |
| SSP-RK3 | `runge_kutta` | 3 | Strong stability preserving |
| Classic RK4 | `runge_kutta` | 4 | A-stable for small CFL |
| Dormand-Prince | `adaptive` | 4(5) | Adaptive step (PI controller) |
| IMEX-RK | `imex` | 2–4 | L-stable diffusion, ERK convection |
| RK-Chebyshev | `rk_chebyshev` | 2 | Extended real-axis stability (s² stages) |
| Exponential | `exponential` | ≥2 | Exact for linear problems |

---

## SIMD Status

**SIMD SpMV is 27–32% slower than scalar** for CSR matrices due to
irregular gather access `x[col_indices[j]]` (Sprint 1.43.0 / 1.55.0 benchmarks).

Recommendation: **Use Rayon parallelism** (`matrix_free/parallel_solvers.rs`)
for 5–20× speedup on multi-core hardware. The `simd` feature flag is disabled by
default; do NOT enable it expecting performance gains on CSR problems.

---

## Key Mathematical Theorems

| Module | Theorem |
|--------|---------|
| `linear_solver/gmres/` | GMRES minimises residual over Kₘ subspace |
| `linear_solver/conjugate_gradient/` | CG converges in ≤n steps for SPD A |
| `linear_solver/preconditioners/ilu/` | ILU stability under appropriate fill-in |
| `linear_solver/preconditioners/multigrid/amg.rs` | AMG O(1) convergence for M-matrices |
| `high_order/weno.rs` | WENO5 5th-order smooth / 3rd-order discontinuous |
| `high_order/spectral/` | Exponential spectral convergence for analytic u |
| `high_order/dg/` | DG discrete entropy stability |
| `time_stepping/stability.rs` | von Neumann / CFL conditions per scheme |

---

## Prohibited Patterns

- No direct physics (fluid properties, BCs) — import from `cfd-core`
- No file I/O — use `cfd-io`
- Avoid introducing new SIMD SpMV paths without benchmark justification
- AMG must remain M-matrix-safe; do not relax coarsening thresholds without proof
