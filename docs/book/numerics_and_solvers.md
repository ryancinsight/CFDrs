# Chapter 5 — Spectral, FEM, MUSCL, and Matrix-Free Methods

## Overview

CFDrs's solver layer decomposes into **discretization** (stencil/FEM/spectral),
**linear algebra** (matrix-free, sparse direct, Krylov), and **time
integration** (explicit, implicit, adaptive). This chapter presents the
mathematical foundations for each discretization family, derives the weak forms
and flux reconstruction schemes, explains the SIMPLEC/PIMPLE pressure-velocity
coupling loop, formalizes matrix-free Jacobian action, states the CFL and
Richardson analysis, and maps each concern to its `cfd-math` / `cfd-2d` /
`cfd-3d` API.

Three discretization surfaces dominate CFDrs:

| Surface | Crate | Backend trait | Formal order |
|---|---|---|---|
| Stencil (finite-volume) | `cfd-2d::stencil`, `cfd-3d::stencil` | `StencilOp<F>` | 2nd central; TVD with MUSCL |
| FEM (P1, P2 tetrahedral) | `cfd-3d::fem` | `FEMOp<F>` | P1 2nd, P2 3rd in $L_2$ |
| Spectral (FFT-based) | `cfd-3d::spectral`, `cfd-2d::spectral` | `SpectralOp<F>` | Exponential on smooth fields |

All three consume the boundary-and-state triple from
[Foundations](foundations.md). Commonalities — trait frontiers, error budgets,
and Atlas backends — live in [`cfd-math`] so they are not redeclared per
solver.

---

## 1 — Finite-Element Weak Form

### Strong to Weak

Consider a scalar advection-diffusion prototype $-\nabla\cdot(k\nabla u - \mathbf{b} u) = f$ in
$\Omega$ with $u = g$ on $\Gamma_D$ and $k\partial u/\partial n = h$ on
$\Gamma_N$. Multiply by a test function $v \in V_0 = \{v \in H^1(\Omega):
v|_{\Gamma_D} = 0\}$ and integrate by parts:

$$\int_\Omega k\nabla u \cdot \nabla v\, d\Omega - \int_\Omega (\mathbf{b}\cdot\nabla u) v\, d\Omega + \int_{\Gamma_N} h v\, d\Gamma = \int_\Omega f v\, d\Omega \tag{1}$$

i.e. seek $u \in V_g = \{u \in H^1(\Omega): u|_{\Gamma_D}=g\}$ such that

$$a(u,v) = \ell(v) \quad \forall v \in V_0 \tag{2}$$

with bilinear form

$$a(u,v) = \int_\Omega k\nabla u \cdot \nabla v\, d\Omega - \int_\Omega (\mathbf{b}\cdot\nabla u) v\, d\Omega \tag{3}$$

and linear functional $\ell(v) = \int_\Omega fv\,d\Omega - \int_{\Gamma_N} hv\,d\Gamma$.

For Stokes flow the mixed weak form uses Taylor-Hood P2-P1 elements:
velocity $ \mathbf{u}_h \in [P_2]^d$, pressure $p_h \in P_1$, satisfying the
inf-sup (LBB) condition:

$$\text{Find } (\mathbf{u}_h, p_h) \in V_h \times Q_h \text{ s.t. } \quad a(\mathbf{u}_h, \mathbf{v}) + b(\mathbf{v}, p_h) = \langle \mathbf{f}, \mathbf{v}\rangle \quad \forall \mathbf{v} \in V_h$$

$$b(\mathbf{u}_h, q) = 0 \quad \forall q \in Q_h \tag{4}$$

with $a(\mathbf{u},\mathbf{v}) = \int_\Omega 2\mu\,\mathbf{S}(\mathbf{u}):\mathbf{S}(\mathbf{v})$ and
$b(\mathbf{v},q) = -\int_\Omega q\,\nabla\cdot\mathbf{v}$. The equal-order P1-P1 pair
violates LBB and requires stabilization (Brezzi-Pitkaranta PSPG term in
`cfd-3d::fem::stabilization`).

### Discrete System

Expanding $u_h = \sum_j U_j \phi_j$ over basis functions $\phi_j$ yields
$K_{ij} = a(\phi_j,\phi_i)$, $b_i = \ell(\phi_i)$ and linear system
$K\mathbf{U} = \mathbf{b}$. CFDrs assembles $K$ via Gauss quadrature over
tetrahedral (Wandzurat) or hexahedral reference elements, with SIMD-vectorized
quadrature loops via `hermes-simd` lanes.

$$K_{ij} = \sum_e \sum_q w_q\, \left[k \nabla\phi_j(\boldsymbol{\xi}_q)\cdot\nabla\phi_i(\boldsymbol{\xi}_q) - (\mathbf{b}\cdot\nabla\phi_j(\boldsymbol{\xi}_q))\phi_i(\boldsymbol{\xi}_q)\right] |J_e(\boldsymbol{\xi}_q)| \tag{5}$$

where $|J_e|$ is the Jacobian determinant mapping reference to physical element.

### Conforming Assembly and Penalty

Strong Dirichlet BCs are enforced via row-zeroing:
$K_{i_G,*}=0$, $K_{i_G,i_G}=1$, $b_{i_G}=g(x_{i_G})$. Weak enforcement via Nitsche
penalty is available for non-conforming cut-cell FEM — see
[geometry chapter](geometry_and_meshing.md).

---

## 2 — Spectral Tau Method

### Basis and Collocation

For periodic $x \in [0, 2\pi]$ or channel with Chebyshev non-periodic
direction $y \in [-1,1]$, CFDrs expands fields as

$$u(x, y, t) = \sum_{k_x=-N_x/2}^{N_x/2-1}\sum_{m=0}^{N_y-1} \hat{u}_{k_x, m}(t)\, T_m(y)\, e^{ik_x x} \tag{6}$$

where $T_m$ is the Chebyshev polynomial of the first kind,
$T_m(y) = \cos(m\arccos y)$. Collocation points are Chebyshev-Gauss-Lobatto
$y_j = \cos(\pi j / (N_y-1))$, clustering at walls for boundary layer
resolution.

### Chebyshev Differentiation Matrix

The first-derivative matrix satisfies $\mathbf{u}' = D \mathbf{u}$ with

$$D_{ij} = \begin{cases} \frac{\bar{c}_i}{\bar{c}_j} \frac{(-1)^{i+j}}{y_i - y_j} & i \ne j \\ -\frac{y_i}{2(1-y_i^2)} & i=j \ne 0, N-1 \\ \frac{2(N-1)^2+1}{6} & i=j=0 \\ -\frac{2(N-1)^2+1}{6} & i=j=N-1 \end{cases} \tag{7}$$

where $\bar{c}_0 = \bar{c}_{N-1} = 2$, $\bar{c}_j =1$ otherwise. Second derivative
$D_2 = D \cdot D$ (with care at boundary rows). CFDrs stores $D$ and $D_2$ as
`leto::NdArray` with atlas backend; they are precomputed once per resolution.

### Tau Equations

In spectral space, applying the Laplacian and projection gives for the 3-D
periodic Poisson $-\nabla^2 \phi = f$:

$$-(k_x^2 + k_z^2)\hat{\phi} + D_2 \hat{\phi} = -\hat{f}, \quad \text{each } (k_x,k_z) \text{ slice} \tag{8}$$

Boundary conditions determine the last two rows: e.g. Dirichlet
$\phi(y=\pm 1)=0$ becomes $\sum_m (\pm 1)^m \hat{\phi}_m = 0$. The tau-correction
discards two highest-mode residual equations and replaces them with boundary
conditions, yielding a pentadiagonal system solved in $O(N_y)$ per mode.

### FFT Plane Decoupling

In $(x, z)$ periodic directions, steps are:

1. Forward 2-D real FFT on each $y$-plane (via `apollo::FftPlan`).
2. For each Fourier mode $(k_x, k_z)$ solve the 1-D tau system along $y$.
3. Inverse FFT the result to physical space.

Cost: $O(N_x N_z \log N_x \log N_z + N_x N_z N_y)$ per solve, as opposed to
$O(N^3)$ for direct 3-D Laplacian. Validated in
[`spectral_3d_poisson`](examples/spectral_3d_poisson.md).

---

## 3 — MUSCL Limiter Derivation

### Finite-Volume Upwind Reconstruction

On a 1-D cell-face $i+1/2$, second-order linear reconstruction from cell $i$
uses slope $s_i$:

$$u_{i+1/2}^L = u_i + \tfrac12 s_i \Delta x, \quad u_{i+1/2}^R = u_{i+1} - \tfrac12 s_{i+1}\Delta x \tag{9}$$

Unlimited slopes $s_i = (u_{i+1} - u_{i-1})/(2\Delta x)$ produce overshoots
near discontinuities (Godunov's theorem: linear second-order schemes are not
monotone). MUSCL (Monotone Upstream-centered Schemes for Conservation Laws,
van Leer 1979) limits slopes:

$$s_i = \mathrm{limiter}\left(\frac{u_i - u_{i-1}}{\Delta x},\ \frac{u_{i+1}-u_i}{\Delta x}\right) \tag{10}$$

where the limiter function $\Phi(r)$ with $r = (\Delta u_L)/(\Delta u_R)$ enforces TVD:

- **minmod**: $\Phi(r) = \max(0, \min(1, r))$ — most diffusive but monotone.
- **van Leer**: $\Phi(r) = (r+|r|)/(1+|r|)$ — smooth.
- **superbee** (Roe): $\Phi(r) = \max(0, \min(2r,1), \min(r,2))$ — compressive.
- **van Albada**: $\Phi(r) = (r+r^2)/(1+r^2)$.
- **Koren**: $\Phi(r) = \max(0, \min(2r, (2+r)/3, 2))$ — third-order upwind TVD.

CFDrs `MusclLimiter` trait:

```rust
pub trait MusclLimiter<F: FloatElement> {
    fn phi(&self, r: F) -> F;
    fn limited_slope(&self, delta_l: F, delta_r: F) -> F {
        if delta_r.abs() < F::EPSILON { return F::zero(); }
        let r = delta_l / delta_r;
        self.phi(r) * delta_r
    }
}
```

Registered implementations: `MinMod`, `VanLeer`, `Superbee`, `VanAlbada`,
`Koren`. The solver consumes `&dyn MusclLimiter<F>` per flux evaluation.

### Flux Evaluation

Riemann flux at the face uses the limited left/right states with a numerical
flux function (Rusanov / Roe / HLLC):

$$F_{i+1/2} = F^{\mathrm{Riemann}}(u_{i+1/2}^L, u_{i+1/2}^R) \tag{11}$$

For incompressible flow the advective flux $F = \mathbf{u} \otimes \mathbf{u}$
is often treated with a central flux for low-Mach plus MUSCL correction for
sharper gradients, balancing kinetic-energy preservation and monotonicity.

---

## 4 — Matrix-Free Operator Application

For high-resolution flows, CFDrs prefers **matrix-free** evaluation: the
linear operator is a trait, and the solver allocates only the action function
$A\cdot x$. Matrix-free applications decouple memory footprint from
discretization size — e.g. a $256^3$ spectral Laplacian would need
$>4$ GB as a sparse CSR but the matrix-free fft-based action uses only
$O(N)$ with $O(N\log N)$ per matvec.

```rust
pub trait LinearOp<F: FloatElement> {
    fn apply(&self, x: &NdArray<F, Ix3>, y: &mut NdArray<F, Ix3>);
    fn diag(&self) -> Option<CowArray<F, Ix1>> { None }
    fn shape(&self) -> (usize, usize);
}
```

A Krylov solver asks only for `apply`, so it composes with any of the
three discretizations.

### Jacobian Action without Jacobian Assembly — JFNK

For nonlinear problems $\mathbf{F}(\mathbf{u})=0$, the Jacobian-vector product
is approximated by directional finite difference:

$$\mathbf{J}\mathbf{v} \equiv \frac{\partial \mathbf{F}}{\partial \mathbf{u}}\mathbf{v} \approx \frac{\mathbf{F}(\mathbf{u} + \varepsilon\mathbf{v}) - \mathbf{F}(\mathbf{u})}{\varepsilon} \tag{12}$$

with

$$\varepsilon = \frac{\sqrt{\epsilon_{\mathrm{mach}}}\, \|\mathbf{u}\|}{ \|\mathbf{v}\| + \sqrt{\epsilon_{\mathrm{mach}}}} \tag{13}$$

($\epsilon_{\mathrm{mach}} = 2.2\times10^{-16}$ for `f64`). This is embedded in
`cfd-math::jfnk::JacobianFreeOp`. Preconditioning uses an approximate Jacobian
(typically inviscid / frozen-viscosity) computed once, while the JFNK
matvec stays exact.

For exact analytical Jacobian action (e.g. in 2-D vorticity-stream), CFDrs
also provides `ExplicitJacobianOp<F>` backed by `leto::CsrMatrix`.

---

## 5 — SIMPLEC and PIMPLE Pressure-Velocity Coupling

### SIMPLEC Loop

The SIMPLE family enforces $\nabla\cdot\mathbf{u}=0$ through a predictor-corrector:

```
// Predictor
1. Solve momentum with current p*:  A_u u* = -∇p* + RHS
2. Decompose diagonal:  A = diag(A) + N  (off-diagonal neighbor contribution)
3. For SIMPLEC consistent approx: d_u = 1 / (A_ii - Σ A_ij) vs SIMPLE d_u = 1/A_ii

// Pressure correction
4. Solve pressure-correction Poisson:  ∇·(d_u ∇p') = ∇·u*
5. Correct velocity:  u = u* - d_u ∇p'
6. Correct pressure:  p = p* + α_p p'   (α_p ≈ 0.3 SIMPLE, ≈ 1.0 SIMPLEC)
7. Repeat until ‖∇·u‖ < ε_continuity
```

The consistent correction in SIMPLEC ($d_u = 1/(A_{ii} - \sum A_{ij})$ rather
than $1/A_{ii}$) allows $\alpha_p \approx 1$ and improves convergence on
non-orthogonal meshes. Formalization in Patankar (1980) with SIMPLEC variant
by van Doormaal & Raithby (1984).

### PIMPLE (PISO + SIMPLE Outer Loops)

For transient flows, PIMPLE nests:

```
for outer_correctors (SIMPLE-like):
    solve momentum predictor
    for inner_correctors (PISO nCorrectors ≈ 2):
        solve pressure-correction
        explicit velocity correction
    solve turbulence / species if present
```

CFDrs `PimpleConfig`:

```rust
pub struct PimpleConfig {
    pub n_outer_correctors: usize,   // typically 1-3
    pub n_inner_correctors: usize,   // PISO inner loops, typically 2
    pub n_non_orthogonal_correctors: usize, // mesh non-orth correction
    pub residual_control: ResidualControl,
}
```

Demonstrated in [`simplec_pimple_demo`](examples/simplec_pimple_demo.md).

---

## 6 — CFL Condition and Adaptive Time-Stepping

### CFL Stability

Explicit time-stepping for advection-diffusion satisfies

$$\mathrm{CFL} = \frac{|u|\Delta t}{\Delta x} \leq \mathrm{CFL}_{\max} \tag{14}$$

$$\mathrm{CFL}_{\mathrm{diff}} = \frac{\nu \Delta t}{\Delta x^2} \leq \frac{1}{2d} \quad (d\text{-dimensional explicit diffusion}) \tag{15}$$

Overall $\Delta t = \min(\mathrm{CFL}_{\max}\Delta x/|u|_{\max}, \,
\mathrm{CFL}_{\max}^{\mathrm{diff}} \Delta x^2 /\nu)$. CFDrs default
$\mathrm{CFL}_{\max} = 0.5$ (central advection with RK4) or $0.8$ (MUSCL upwind).

### Adaptive Stepping

```rust
let stepping = AdaptiveStepping::builder()
    .cfl_max(0.5)
    .atol_field(1e-6)
    .build();

stepping.run(&mut solver, &state)?;
```

The adaptive scheme refines $\Delta t$ until local truncation error falls below
`atol_field` **per cell** via a PI controller:

$$\Delta t_{n+1} = \Delta t_n \left(\frac{\mathrm{tol}}{\mathrm{err}_n}\right)^{k_P} \left(\frac{\mathrm{err}_{n-1}}{\mathrm{err}_n}\right)^{k_I} \tag{16}$$

with $k_P = 0.3$, $k_I = 0.1$ (Gustafsson controller). The `err` is the
embedded RK pair difference $\|\mathbf{u}^{(p)} - \mathbf{u}^{(p-1)}\|_\infty$.

---

## 7 — Iterative Solvers

`cfd-math::krylov` exposes BiCGSTAB, GMRES, and IDR(s) behind a common trait:

```rust
pub trait KrylovSolver<F: FloatElement> {
    fn solve(&self, op: &dyn LinearOp<F>, rhs: &NdArray<F, Ix1>,
             precond: &dyn Preconditioner<F>, tol: F, max_iter: usize)
             -> KrylovResult<F>;
}
```

All three reach the same convergence target at the same iteration count
regardless of which discretization backs them.

```rust
let krylov = BiCgStab::<f64>::new();
let solution = krylov.solve(&op, &rhs, &preconditioner, 1e-6, 5000)?;
```

Preconditioners: Jacobi (`DiagJacobi`), ILU(0) (`Ilu0`), SSOR, and multigrid
`V-cycle` for Poisson-dominated blocks. The pressure-correction system is
typically solved with CG or BiCGSTAB with a multigrid preconditioner; momentum
systems with GMRES + ILU(0).

### Convergence Tolerances

Solver tolerances are relative residual norms:

$$ \frac{\|A\mathbf{x}_k - \mathbf{b}\|}{\|\mathbf{b}\| + \|\mathbf{b}_0\|} < \mathrm{tol} \tag{17}$$

with defaults $\mathrm{tol}_{p} = 10^{-7}$ for pressure Poisson and
$\mathrm{tol}_{u} = 10^{-6}$ for momentum, consistent with the Atlas trait
frontier so callers swap solvers without re-tuning.

---

## 8 — Richardson Extrapolation and Grid Convergence

For three grids with refinement ratio $r = h_2/h_1 = h_3/h_2 < 1$, order

$$p = \frac{\ln((f_2 - f_3)/(f_1 - f_2))}{\ln r}, \quad f_{\mathrm{exact}} \approx f_1 + \frac{f_1 - f_2}{r^p - 1} \tag{18}$$

and Grid Convergence Index (Roache 1994):

$$\mathrm{GCI}_{12} = \frac{F_s |f_1 - f_2| / |f_1|}{r^p - 1}, \quad F_s = 1.25 \text{ (three-grid), } 3.0 \text{ (two-grid)} \tag{19}$$

CFDrs runs (18)-(19) automatically in `cfd-validation::richardson`.

---

## API Mapping

| Concern | Path |
|---|---|
| Stencil ops | `cfd-2d::stencil::{StencilOp, CentralDiff, UpwindMUSCL}` |
| FEM assembly | `cfd-3d::fem::{FemAssembler, P1Tet, P2Tet, TaylorHoodP2P1}` |
| Spectral ops | `cfd-3d::spectral::{SpectralOp, ChebyshevDiffMat, TauSolver}` |
| FFT | `apollo::FftPlan` (Atlas) / `cfd-math::spectral` |
| MUSCL | `cfd-2d::physics::muscl::{MusclLimiter, MinMod, VanLeer, Superbee, Koren}` |
| Matrix-free | `cfd-math::linear_op::{LinearOp, JacobianFreeOp}` |
| SIMPLEC/PIMPLE | `cfd-3d::pressure::{SimplecConfig, PimpleConfig, PimpleLoop}` |
| Krylov | `cfd-math::krylov::{BiCgStab, Gmres, Idrs, Preconditioner}` |
| Adaptive time | `cfd-2d::adaptive::AdaptiveStepping`, `cfd-3d::adaptive` |
| Richardson | `cfd-validation::richardson::{RichardsonExtrapolation, GCI}` |

---

## Examples Referenced by This Chapter

- [Example: spectral_3d_poisson](examples/spectral_3d_poisson.md) — spectral Poisson with FFT-plane decoupling and tau solve.
- [Example: fem_3d_stokes](examples/fem_3d_stokes.md) — 3-D Stokes FEM, P1 tetrahedral and Taylor-Hood P2-P1.
- [Example: matrix_free_demo](examples/matrix_free_demo.md) — matrix-free operator application and JFNK Jacobian-free Newton–Krylov.
- [Example: muscl_schemes_demo](examples/muscl_schemes_demo.md) — MUSCL limiters (minmod through Koren) for the 2-D solver.
- [Example: simplec_pimple_demo](examples/simplec_pimple_demo.md) — SIMPLEC and PIMPLE pressure-correction families; convergence history.
- [Example: 2d_heat_diffusion](examples/2d_heat_diffusion.md) — explicit diffusion on a coarse grid; CFL adaptation demo.

## Further Reading

- Zienkiewicz & Taylor, *The Finite Element Method*, Vol.1, 7th ed. (2013) — FEM weak form.
- Canuto et al., *Spectral Methods: Fundamentals in Single Domains* (2006) — Chebyshev tau method.
- van Leer, JCP 32, 101 (1979) — original MUSCL.
- Patankar, *Numerical Heat Transfer and Fluid Flow* (1980) — SIMPLE loop.
- van Doormaal & Raithby, Numer. Heat Transfer 7(2), 147 (1984) — SIMPLEC.
- Issa, JCP 62, 65 (1986) — PISO/PIMPLE.
- Knoll & Keyes, JCP 193, 357 (2004) — JFNK review.
- `cfd-math` source: `crates/cfd-math/src/{krylov,linear_op,jfnk,spectral}.rs`.
- [Geometry chapter](geometry_and_meshing.md) for boundary generation upstream of solvers.
- [Core flow benchmarks](core_flows.md) for Richardson extrapolation application.
