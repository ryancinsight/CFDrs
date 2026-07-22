# Chapter 2 — Canonical Incompressible Benchmarks

## Overview

Canonical incompressible flows are the regression anchors of CFDrs. They
provide **analytical or high-fidelity numerical oracles** against which every
discretization, solver backend, and Atlas migration stage is validated. This
chapter covers three benchmarks that span the Reynolds-number range and geometric
complexity that CFDrs must handle: lid-driven cavity (Ghia), Hagen-Poiseuille
pipe flow, and turbulent channel references. Each section states the governing
equations specialized to the benchmark, derives or cites the analytical profile
where it exists, describes vortex / boundary-layer physics, and maps to the
`cfd-core` and `cfd-validation` APIs.

---

## Governing Setup

All three benchmarks solve the incompressible Navier-Stokes system

$$\nabla \cdot \mathbf{u} = 0 \tag{1}$$

$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu \nabla^2\mathbf{u} \tag{2}$$

in nondimensional form with $\mathrm{Re} = UL/\nu$ the governing similarity
parameter. The three cases differ in domain, boundary conditions, and which
terms balance in (2):

| Benchmark | Domain | Re range | Dominant balance |
|---|---|---|---|
| Lid-driven cavity | Unit square, lid moving at $U=1$ | 100–10000 | Advection ↔ diffusion, secondary vortices |
| Hagen-Poiseuille pipe | Circular cross-section, periodic or pressure-driven | 100–2100 (laminar) | Pressure gradient ↔ viscous diffusion |
| Turbulent channel | Plane channel $[-h, h] \times [0, L_x] \times [0, L_z]$ | $Re_\tau = 180$–590 | Mean shear ↔ Reynolds stress |

---

## 1 — Lid-Driven Cavity (Ghia et al. 1982)

### Problem Definition

The cavity is $\Omega = [0,1]^2$ with boundaries

- Top wall ($y=1$): $\mathbf{u} = (U_{\mathrm{lid}}, 0) = (1, 0)$ — the moving lid.
- Bottom, left, right walls: $\mathbf{u} = \mathbf{0}$ — no-slip.
- No body force: $\mathbf{f} = 0$.

The steady-state form of (2) for $\mathrm{Re}$ moderate ($<1000$) has a unique
stable solution; at higher $\mathrm{Re}$ multiple steady states and Hopf
bifurcations appear. CFDrs validates at $\mathrm{Re} = 100, 400, 1000$.

### Vorticity-Streamfunction Formulation

In 2-D the system reduces to

$$\frac{\partial \omega}{\partial t} + u \frac{\partial \omega}{\partial x} + v \frac{\partial \omega}{\partial y} = \frac{1}{\mathrm{Re}}\nabla^2\omega \tag{3}$$

$$\nabla^2\psi = -\omega, \quad u = \frac{\partial\psi}{\partial y}, \; v = -\frac{\partial\psi}{\partial x} \tag{4}$$

Thom's first-order wall vorticity condition at a no-slip wall:

$$\omega_{\mathrm{wall}} = -\frac{2\psi_{m+1}}{\Delta n^2} \tag{5}$$

Briley's second-order condition improves accuracy for high-Re corners. CFDrs
implements both with a trait switch:

```rust
pub enum WallVorticity { Thom, BrileySecondOrder }

let solver_cfg = VorticityStreamConfig {
    wall_vorticity: WallVorticity::BrileySecondOrder,
    re: 400.0,
    grid: (129, 129),
};
```

### Vortex Formation

At $\mathrm{Re} = 100$, a single primary vortex fills the cavity with two
small corner vortices (bottom-left, bottom-right). As $\mathrm{Re}$
increases:

- $\mathrm{Re} = 400$: bottom corner vortices enlarge; a weak tertiary vortex
  appears in the bottom-right corner at very high resolution.
- $\mathrm{Re} = 1000$: secondary vortices are prominent; the centerline
  velocity profile develops inflection points.
- $\mathrm{Re} > 5000$: multiple secondary and tertiary vortices; boundary
  layers thin as $\delta \sim \mathrm{Re}^{-1/2}$; unsteady transition begins.

CFDrs grid requirements scale as $\Delta x \sim \delta$ in the wall-normal
direction, so non-uniform clustering (hyperbolic tangent or geometric) is needed
above $\mathrm{Re} \approx 2000$. The grid utility `StructuredGrid2D::clustered`
provides this.

### Ghia Oracle

`cfd-validation` ships tabulated centerline velocity profiles from Table I–III
of Ghia et al. (1982) at $\mathrm{Re} = 100, 400, 1000$:

```rust
use cfd_validation::ghia::{GhiaOracle, Centerline};

let oracle = GhiaOracle::re_400();
let u_profile = solver.velocity_along(Centerline::VerticalMid);
let err_l2 = oracle.l2_error(&u_profile);
assert!(err_l2 < 5e-3, "Ghia L2 deviation too large: {err_l2}");
```

The oracle provides:

- $u$ along $x = 0.5$ (vertical centerline) and $v$ along $y = 0.5$ (horizontal).
- Primary vortex center $(x_c, y_c)$ and streamfunction $\psi_c$ at that point.
- Published tolerances: $L_2 < 5\times 10^{-3}$ on $129 \times 129$ grid,
  $L_2 < 1\times 10^{-3}$ on $257 \times 257$.

Convergence order is verified through Richardson extrapolation (see
[numerics_and_solvers](numerics_and_solvers.md)).

### API Mapping

| Concern | Crate / path |
|---|---|
| Grid | `cfd-2d::grid::StructuredGrid2D::unit_square` |
| Solver | `cfd-2d::physics::vorticity_stream::VorticityStreamSolver` |
| BC types | `cfd_core::boundary::{BoundaryKind::Lid, BoundaryKind::Wall}` |
| State | `cfd_core::state::StateVec` |
| Oracle | `cfd-validation::ghia::GhiaOracle` |
| Extrapolation | `cfd-validation::richardson::RichardsonExtrapolation` |

---

## 2 — Hagen-Poiseuille Pipe Flow

### Problem Definition

Fully developed laminar flow in a straight circular pipe of radius $R$,
length $L \gg R$, driven by a constant pressure gradient $\partial p / \partial
z = -G$ (or prescribed volumetric flowrate $Q$). Axisymmetric, steady,
incompressible with no swirl. Cylindrical coordinates $(r, \theta, z)$ with $z$
along the pipe axis.

Under these assumptions $u_r = u_\theta = 0$, $\partial/\partial z = 0$ for
fully developed flow (neglecting entrance length), and (2) reduces in the
$z$-direction to

$$0 = -\frac{1}{\rho}\frac{\partial p}{\partial z} + \nu \frac{1}{r}\frac{\partial}{\partial r}\left(r\frac{\partial u_z}{\partial r}\right) \tag{6}$$

i.e.

$$\frac{1}{r}\frac{d}{dr}\left(r\frac{du_z}{dr}\right) = -\frac{G}{\mu} \tag{7}$$

### Parabolic Profile Derivation

Integrate (7):

$$r\frac{du_z}{dr} = -\frac{G}{2\mu} r^2 + C_1$$

Finiteness at $r=0$ requires $C_1 = 0$. Integrate once more:

$$u_z(r) = -\frac{G}{4\mu} r^2 + C_2$$

No-slip at $r=R$ gives $u_z(R)=0$, hence $C_2 = G R^2 / (4\mu)$. Therefore

$$u_z(r) = \frac{G}{4\mu}\left(R^2 - r^2\right) = 2U_{\mathrm{avg}}\left(1 - \frac{r^2}{R^2}\right) \tag{8}$$

where $U_{\mathrm{avg}} = G R^2 / (8\mu)$ is the mean velocity. The centerline
(maximum) velocity $U_{\max} = G R^2 / (4\mu) = 2 U_{\mathrm{avg}}$.

Volumetric flowrate

$$Q = \int_0^R u_z(r) 2\pi r\, dr = \frac{\pi G R^4}{8\mu} = \frac{\pi R^2 U_{\max}}{2} \tag{9}$$

Pressure drop for a pipe of length $L$:

$$\Delta p = G L = \frac{8\mu L Q}{\pi R^4} = \frac{32 \mu L U_{\mathrm{avg}}}{D^2} \tag{10}$$

with $D = 2R$. In friction-factor form using $\mathrm{Re}_D = U_{\mathrm{avg}} D / \nu$,

$$f = \frac{\Delta p \, D}{\frac{1}{2}\rho U_{\mathrm{avg}}^2 L} = \frac{64}{\mathrm{Re}_D} \tag{11}$$

which is the laminar Moody-chart branch that `cfd-validation` checks against.

### Entrance Length

The fully developed profile assumes $L$ exceeds the hydrodynamic entrance length
$L_h$. For laminar pipe flow

$$L_h \approx 0.06 \, \mathrm{Re}_D \, D \tag{12}$$

CFDrs slender-channel examples set the domain length to $2 L_h$ to ensure the
validation slice is taken in the fully developed region.

### API Mapping and Oracle

```rust
use cfd_validation::poiseuille::{PoiseuilleOracle, PipeParams};

let params = PipeParams { radius: 1.0, grad_p: 8.0, mu: 1.0, rho: 1.0 };
let oracle = PoiseuilleOracle::new(params);
let profile = solver.radial_velocity_profile(z_slice = 0.5);
let err = oracle.max_abs_error(&profile);
```

`PoiseuilleOracle` computes (8) analytically at the grid points used by the
solver, so the comparison is pointwise exact up to solver discretization error.

| Concern | Crate / path |
|---|---|
| Grid | `cfd-1d::grid::UniformGrid1D` (radial), `cfd-2d::grid::PipeGrid` |
| Solver | `cfd-2d::physics::poiseuille`, `cfd-1d::poiseuille` |
| BC | `Inlet { parabolic }`, `OutflowPressure`, `Symmetry` at $r=0$ |
| Oracle | `cfd-validation::poiseuille::PoiseuilleOracle` |

---

## 3 — Turbulent Channel Flow References

### Setup

Plane channel: half-height $h$, friction Reynolds number

$$\mathrm{Re}_\tau = \frac{u_\tau h}{\nu}, \quad u_\tau = \sqrt{\tau_w / \rho} \tag{13}$$

where $\tau_w = \mu \, \partial U / \partial y|_{y=\pm h}$ is wall shear
stress. CFDrs validates against canonical DNS databases:

- Moser, Kim & Mansour, Phys. Fluids 11(4):943 (1999) at $\mathrm{Re}_\tau =
  180, 395, 590$.
- Iwamoto et al. (2002) and Hoyas & Jiménez (2006) at higher $\mathrm{Re}_\tau$.

The domain is periodic in streamwise ($x$) and spanwise ($z$) directions, with
no-slip at $y = \pm h$. The mean streamwise velocity $U(y) = \langle u_x \rangle_{x,z,t}$
exhibits:

- Viscous sublayer $y^+ < 5$: $U^+ = y^+$, where $y^+ = y u_\tau / \nu$ and $U^+ = U / u_\tau$.
- Buffer layer $5 < y^+ < 30$: strong curvature.
- Log-law $30 < y^+ < 0.2 \mathrm{Re}_\tau$: $U^+ = \frac{1}{\kappa}\ln y^+ + B$ with $\kappa \approx 0.41$, $B \approx 5.2$.

### Validation Metrics

`cfd-validation::turbulent_channel` checks:

```rust
let dns = DnsChannel::moser_re_tau_180();
let stats = solver.statistics_after_transient(transient_steps = 50000);

// Mean profile
let err_mean = L2Metric::compute(&stats.u_mean, &dns.u_mean_profile());
assert!(err_mean < 0.05, "mean profile deviation");

// Friction coefficient
let cf = stats.friction_coefficient();
assert!((cf - dns.cf_ref()).abs() / dns.cf_ref() < 0.02);
```

Quantities checked include $U^+(y^+)$, $\langle u'_i u'_j \rangle$, TKE
budget, and centerline Reynolds number. Spatial convergence under LES uses the
$\Delta \sim (\Delta_x \Delta_y \Delta_z)^{1/3}$ metric.

### Connection to Turbulence Models

The channel flow case is the primary testbed for
[turbulence model validation](turbulence_multiphase.md): the $k$-$\epsilon$
RANS closure recovers the log law only with wall functions; $k$-$\omega$ SST
improves near-wall accuracy; Smagorinsky LES with Van Driest damping recovers
the buffer layer but needs resolution $y^+_1 < 1$ at the wall for wall-resolved
LES. See [turbulence chapter](turbulence_multiphase.md) for closure details.

---

## Convergence and Richardson Extrapolation

For both cavity and pipe benchmarks, CFDrs establishes spatial order $p$ via
Richardson extrapolation on three grids with refinement ratio $r$:

$$p = \frac{\ln\frac{f_2 - f_3}{f_1 - f_2}}{\ln r}, \quad f_{\mathrm{exact}} \approx f_1 + \frac{f_1 - f_2}{r^p - 1} \tag{14}$$

where $f_1, f_2, f_3$ are values on fine, medium, coarse grids. For the cavity
solver using second-order central differences and first-order Thom BC, global
second order is expected ($p \approx 2$). The
`comprehensive_validation_suite` in `cfd-validation` runs (14) automatically.

---

## Examples Referenced by This Chapter

- [Example: cavity_validation](examples/cavity_validation.md) — lid-driven cavity against Ghia oracle at $\mathrm{Re}=100,400,1000$; vortex center tracking.
- [Example: pipe_flow_validation](examples/pipe_flow_validation.md) — Hagen-Poiseuille against analytical parabolic profile; friction factor check.
- [Example: turbulent_channel_flow](examples/turbulent_channel_flow.md) — plane channel TKE budget vs Moser DNS.

## Further Reading

- Ghia, Ghia & Shin, JCP 48 (1982).
- Moser, Kim & Mansour, Phys. Fluids 11(4):943 (1999).
- Pope, *Turbulent Flows*, Cambridge (2000) — Chapter 7 for channel flow scaling.
- `cfd-validation` source: `crates/cfd-validation/src/{ghia,poiseuille,turbulent_channel}.rs`.
- [Turbulence and Cavitation](turbulence_multiphase.md) for closure models.
- [Numerics and Solvers](numerics_and_solvers.md) for discretization order.
