# Chapter 19 — 3-D Flows

## Overview

Three-dimensional CFD examples live in `cfd-3d`. They exercise the full 3-D
pressure-velocity coupling, turbulence, cavitation, spectral and FEM
discretization pathways, and the Atlas GPU opt-in behind `hephaestus`. Where
`cfd-2d` validates with Ghia-class 2-D cases, `cfd-3d` introduces genuinely
tri-dimensional physics: spanwise variation in bifurcations, Dean-vortex
torsion in serpentine channels, cavitation inception loci that span a
cross-section, and spectral Poisson solves coupling all three spatial
dimensions via FFT-plane decoupling and Chebyshev tau solves.

This chapter documents the four canonical 3-D examples shipped as cargo
examples in `cfd-3d`, their governing equations, validation oracles, tested
invariants, convergence order claims, and Atlas migration status.

---

## Crate Layout — `cfd-3d`

```
crates/cfd-3d/
  src/
    lib.rs               — crate root, feature gates
    grid.rs              — structured + tet mesh I/O
    mesh/                — TetMesh, quality, laplacian smoothing
    stencil/             — 3-D finite-volume stencil operators (Atlas + legacy)
    fem/
      mod.rs
      p1.rs              — P1 tet stiffness + mass assembly
      p2.rs              — P2 tet
      taylor_hood.rs     — Taylor-Hood P2-P1 Stokes pair
      stabilization.rs   — Brezzi-Pitkaranta PSPG for P1-P1
    spectral/
      mod.rs
      diff_mat.rs        — Chebyshev differentiation matrix
      tau.rs             — spectral tau solver (tridiag/pentadiag)
      fft_bridge.rs      — apollo::FftPlan integration
      poisson.rs         — 3-D periodic / channel spectral Poisson
    turbulence/
      kepsilon.rs
      komega_sst.rs
      smagorinsky.rs
      traits.rs          — TurbulenceModel trait
    multiphase/
      vof.rs             — VOF advection with interface compression
      exchange.rs        — drag / momentum exchange
    cavitation/
      rayleigh_plesset.rs
      eulerian_eulerian.rs
      damage.rs
    pressure/
      simplec.rs
      pimple.rs
    time/
      ab2.rs, crank_nicolson.rs, adaptive.rs
    parallel/
      domain_decomp.rs   — 1-D/2-D block decomposition
      mpi_bridge.rs      — optional MPI via moirai::Executor
```

Atlas-owned types in `cfd-3d` are gated as:

- `leto::NdArray<F, Ix3>` — every field (`u, v, w, p, omega`).
- `leto::Point3/Floor / Vector3` — mesh vertex positions (Atlas type after migration).
- `hermes-simd::Avx2Lane<F>` / `NeonLane<F>` — vectorized FEM quadrature and stencil flux.
- `apollo::FftPlan` — spectral FFT calls.
- `moirai::Executor` — domain decomposition and parallel transport assembly.
- `hephaestus::Backend` — GPU data residency opt-in via `gpu` feature.

---

## 1 — `spectral_poisson_3d` — Spectral Poisson on a Periodic/Channel Domain

### Governing Equation

$$-\nabla^2 \phi = f \quad \text{in } \Omega = [0, 2\pi]^2 \times [-1, 1], \qquad \phi|_{y=\pm 1} = 0 \tag{1}$$

Periodic in $x, z$ (Fourier), non-periodic in $y$ (Chebyshev). Right-hand side
$f$ manufactured to produce a known $\phi$ for convergence verification.

### Discretization

For each Fourier mode $(k_x, k_z)$ define

$$-(k_x^2\mathbb{I} + k_z^2\mathbb{I} - D_2)\hat{\phi}_{k_x,k_z} = \hat{f}_{k_x,k_z} \tag{2}$$

where $D_2$ is the Chebyshev second-derivative matrix (see
[numerics chapter](numerics_and_solvers.md)). Boundary conditions modeled as
the last two equations of the tau system ($y=\pm 1$ Dirichlet). Solved as a
pentadiagonal linear system $O(N_y)$ per mode.

The FFT plane decoupling path:

```
physical f  ──[2-D real FFT x,z]──►  f̂(kx,kz,y_j)
                               ──[Chebyshev tau solve y]──►  φ̂(kx,kz,y_j)
                               ──[2-D real iFFT x,z]─────►  physical φ
```

Complexity $O(N_x N_z \log N_x \log N_z + N_x N_z N_y)$ per solve vs
$O(N^3)$ for direct 3-D Laplacian.

### Tested Invariants

- **Analytical Manufactured Solution**: $\phi_{\mathrm{exact}}(x,y,z) = \sin x
  \cos(\pi y/2) \sin z$ recovers $f_{\mathrm{exact}} = -\nabla^2 \phi_{\mathrm{exact}}$
  exactly; pointwise error $\|\phi_h - \phi_{\mathrm{exact}}\|_\infty < 10^{-10}$
  on $N_x=N_z=32, N_y=65$ (spectral super-convergence).

- **Exponential convergence**: error $\epsilon(N) \sim e^{-cN}$ for smooth
  $f$; validated by doubling $N_y$ and plotting $\log \epsilon$ vs $N$.

- **Condition number**: spectral differentiation condition $\kappa(D_2) \sim
  N_y^4$ managed via diagonal preconditioner for large $N_y$.

### Convergence Order

Spectral: exponential (formally infinite-order) on $C^\infty$ right-hand sides.
In practice $p \approx 10$–$16$ on log-error vs $N$ for moderate $f$ smoothness.
Validated in [`spectral_performance` example](examples/spectral_performance.md).

### Atlas Status

- `apollo::FftPlan` — **graduating** (Atlas path equals legacy `rustfft` path
  within $1\%$ throughput; switching is a type-level flag).
- Chebyshev $D$ / $D_2$ stored as `leto::NdArray` — **in-progress** (legacy
  `nalgebra::DMatrix` being converted).

---

## 2 — `bifurcation_3d_blood` — Non-Newtonian Blood in a Symmetric 3-D Bifurcation

### Governing Setup

Domain: Y-shaped pipe bifurcation built via `BifurcationBuilder` with optional
stenosis tag $S$. Equations:

$$\nabla\cdot\mathbf{u}=0, \quad \rho\frac{D\mathbf{u}}{Dt}= -\nabla p + \nabla\cdot[2\mu(\dot{\gamma})\mathbf{S}(\mathbf{u})] \tag{3}$$

with $\mu(\dot{\gamma})$ Casson or Carreau-Yasuda (see
[biomedical chapter](biomedical_flows.md)). Reynolds number
$\mathrm{Re}_{\mathrm{eff}} \approx 100$–$500$ (parent inlet).

Boundary conditions: parabolic inlet profile (`Inlet{ parabolic: ParabolicProfile
{ umax } }`), Murray-law flow split at child exits ($Q_a/Q_b = r_a^3/r_b^3$),
no-slip walls with optional Fahraeus-Lindqvist viscosity correction for
$D < 300\,\mu\mathrm{m}$ channels.

### Tested Invariants

- **Murray law pressure**: in non-stenosed, Newtonian cases, pressure drop per
  edge scales with Poiseuille law and flow split matches (13)-(15) of
  [biomedical chapter](biomedical_flows.md) within $1\%$.

- **Newtonian limit**: Casson with $\tau_y\to0$ or Carreau with $\lambda\to0$
  recovers Newtonian profile; $L_\infty$ difference $<10^{-6}$ verified.

- **Wall shear stress distribution**: peaks at inner walls of child bifurcation
  separation zones; TAWSS/OSI maps for pulsatile runs flag atherosclerotic
  regions (see biomedical chapter).

- **Mass conservation**: $\sum_{\mathrm{exits}} Q_i = Q_{\mathrm{inlet}}$ within
  solver tolerance $10^{-8}$ (nondim).

### Convergence Order

FEM P1: $O(h^2)$ in $L_2$ energy norm (linear elliptic theory); P2 or
Taylor-Hood P2-P1 Stokes: $O(h^3)$ for velocity, $O(h^2)$ for pressure.

### Atlas Status

- `leto::Point3` for mesh vertices — **bulk** (mid-migration, fallback to
  `nalgebra::Point3` via feature).
- `hermes-simd` quadrature — **in-progress**.
- Rheology dispatch — **complete** at the `FloatElement` seam.

---

## 3 — `venturi_3d_cavitation` — Cavitation Inception in a 3-D Venturi

### Governing Setup

Venturi nozzle built as two frustum cones + cylindrical throat via
`VenturiBuilder`. Equations: mixture form

$$\nabla\cdot(\rho_m\mathbf{u})=0, \quad \rho_m\frac{D\mathbf{u}}{Dt}= -\nabla p + \nabla\cdot(\mu_m(\nabla\mathbf{u}+\nabla\mathbf{u}^T)) \tag{4}$$

$$\frac{\partial\alpha}{\partial t} + \mathbf{u}\cdot\nabla\alpha = \Gamma(\alpha,p) \tag{5}$$

where $\alpha$ is vapour volume fraction, $\rho_m = \alpha\rho_v +
(1-\alpha)\rho_l$, and $\Gamma = \dot{m}^+ - \dot{m}^-$ is the Schnerr-Sauer
mass transfer

$$\dot{m}^+ = C_{\mathrm{evap}}\frac{\rho_v\rho_l}{\rho_m}\alpha(1-\alpha)\sqrt{\frac{2}{3}\frac{p_v-p}{\rho_l}}\,\mathbf{1}_{p<p_v} \tag{6}$$

$$\dot{m}^- = C_{\mathrm{cond}}\frac{\rho_v\rho_l}{\rho_m}\alpha(1-\alpha)\sqrt{\frac{2}{3}\frac{p-p_v}{\rho_l}}\,\mathbf{1}_{p>p_v} \tag{7}$$

with $C_{\mathrm{evap}}=50$, $C_{\mathrm{cond}}=0.01$ (defaults), from Zwart
et al. (2004). Cavitation number $\sigma = (p_\infty-p_v)/(0.5\rho U^2)$.

### Tested Invariants

- **Cavitation inception locus**: vapour first appears at throat minimum
  pressure point; inception cavitation number $\sigma_i$ agrees with
  Franc & Michel experiment Figure 4.12 within $5\%$ for given geometry.

- **Pressure recovery**: pressure recovery length $L_r$ beyond throat matches
  venturi benchmark $L_r \approx 5$–$8 D_{\mathrm{throat}}$.

- **Rayleigh-Plesset consistency**: single-bubble collapse time (equation
  (5) of [turbulence chapter](turbulence_multiphase.md)) matches numerically
  resolved bubble trajectory within $10\%$.

- **Damage index**: $\int \dot{D} dt$ converged spatially and temporally with
  resolution refinement; worst-cell damage correlates with cavitation hot-spots.

### Convergence Order

Second-order central + MUSCL for VOF advection: $p\approx1$–$1.5$ for $\alpha$
(sharp interface smearing), $p\approx2$ for pressure and velocity excluding the
interface zone. Improved by interface compression
($\mathbf{u}_c$ term in Weller VOF, Section 3 of turbulence chapter).

### Atlas Status

- VOF field — `leto::NdArray<F, Ix3>`, **in-progress**.
- Cavitation mass transfer in Atlas kernels — **bulk** (mid-migration).

---

## 4 — `serpentine_3d_dean` — Dean-Vortex Mixing in a 3-D Serpentine

### Governing Setup

Serpentine channel built as chain of capsules (straight + half-torus turns) via
`SerpentineBuilder`. Equations: incompressible N-S plus scalar concentration
equation for mixing diagnostics

$$\frac{\partial c}{\partial t} + \mathbf{u}\cdot\nabla c = D_m \nabla^2 c, \quad c=0,1\text{ at two inlets} \tag{8}$$

with molecular diffusivity $D_m$. Dean number

$$\mathrm{De} = \mathrm{Re}\sqrt{D_h/R_c} \tag{9}$$

with hydraulic diameter $D_h = 4A/P$, $R_c$ radius of curvature of the
U-turn. $\mathrm{De}>100$ triggers a pair of counter-rotating Dean vortices
in each bend, enhancing transverse mixing.

### Tested Invariants

- **Dean vortex topology**: two counter-rotating secondary vortices in each
  U-turn for $\mathrm{De} > 40$, quantified by streamfunction of secondary
  flow in cross-section; intensity $\propto \mathrm{De}^2$ for
  $\mathrm{De} < 100$.

- **Mixing efficiency**: $M = 1 - \sigma_c/\max\sigma_c$ where $\sigma_c$ is
  concentration variance at the outlet plane; $M$ grows as $1-e^{-\alpha N_{\mathrm{turns}}}$
  with Dean-enhanced $\alpha$.

- **Pressure drop vs Dean number**: $\Delta p \propto \mathrm{De}^{1/3}$
  scaling at high De validated against 2-D axisymmetric torus reference
  (Deo et al. 1975).

- **Newtonian limit** (microfluidic variant): same governing equations apply
  with $D_h \sim 100$–$250\,\mu\mathrm{m}$, $Re \sim 1$–$10$, negligible
  Dean vortex — Stokes regime dominated by entrance effects.

### Convergence Order

Second-order central, same as 2-D stencil path (same kernel). Periodicity in
straight sections relaxed by mixing progress.

### Atlas Status

- `moirai::Executor::scope` parallel assembly — **in-progress**.
- Concentration equation — `leto::NdArray`, **bulk**.

---

## Common Validation Infrastructure for 3-D Examples

All four examples share:

```bash
cargo run -p cfd-3d --example spectral_poisson_3d
cargo run -p cfd-3d --example bifurcation_3d_blood
cargo run -p cfd-3d --example venturi_3d_cavitation
cargo run -p cfd-3d --example serpentine_3d_dean
```

and report via `cfd-io::VtkWriter` for ParaView visual inspection, plus JSON
summary for CI. Each validates against an oracle defined in `cfd-validation`
(or `cfd-3d::validation_tests` internal module) and participates in the
`comprehensive_validation_suite` VTK regression set.

---

## Further Reading

- Boyd, *Chebyshev and Fourier Spectral Methods*, 2nd ed., Dover (2001).
- Franc & Michel, *Fundamentals of Cavitation*, Kluwer (2005) — Ch.4 venturi.
- Deo, Streamwise vortices (Dean) — Dean, Proc. Roy. Soc. A 121, 402 (1928).
- `cfd-3d` source: `crates/cfd-3d/src/{fem,spectral,turbulence,cavitation,stencil}.rs`.
- [Numerics chapter](numerics_and_solvers.md) — Chebyshev tau and CFL.
- [Turbulence chapter](turbulence_multiphase.md) — closures and Rayleigh-Plesset.
- [Biomedical chapter](biomedical_flows.md) — Casson / Carreau-Yasuda parameters.
- [Geometry chapter](geometry_and_meshing.md) — CSG pipeline for these examples.
- [Validation suite chapter](crate_validation.md) for Richardson and oracle theory.
