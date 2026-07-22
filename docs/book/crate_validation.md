# Chapter 18 — Validation Suite

## Overview

The `cfd-validation` crate is CFDrs's **orthogonal oracle harness**: it
provides tabulated reference data (Ghia cavity, Hagen-Poiseuille, DNS channel,
Rayleigh-Plesset collapse, Womersley) and verification utilities (Richardson
extrapolation, GCI) that every solver must pass. The harness is intentionally
**crate-agnostic**: it depends only on `cfd-core` scalar traits and `leto`
arrays, never on `cfd-2d` / `cfd-3d` internals. This ensures the oracle
definition cannot drift to accommodate a solver bug.

This chapter documents the three crate-level validation examples that run
inside `cfd-validation` and the broader set of oracles available as library
API.

---

## Crate Layout

```
crates/cfd-validation/
  src/
    lib.rs                  — re-exports all oracles
    ghia.rs                 — lid-driven cavity (Ghia et al. 1982)
    poiseuille.rs           — Hagen-Poiseuille parabolic profile
    turbulent_channel.rs    — Moser DNS RANS statistics
    rayleigh_plesset.rs     — single-bubble Rayleigh collapse time
    womersley.rs            — Womersley pulsatile profile for biomedical
    cavitation.rs           — Franc & Michel venturi inception data
    richardson.rs           — Richardson extrapolation + GCI
    tolerance.rs            — adaptive tolerance scaling
    report.rs               — consolidated validation report
```

Each module exposes the oracle as an owned struct with a typed constructor
and a comparison method that returns a quantified residual — never a raw
table. The report module aggregates multiple oracle results into a JSON report
written to `validation_reports/`.

---

## 1 — Ghia Lid-Driven Cavity Oracle

### Governing Setup

Lid-driven cavity on $\Omega = [0,1]^2$, lid velocity $U_{\mathrm{lid}} = 1$,
no-slip elsewhere, $\mathrm{Re} = 100, 400, 1000$ validated.

Ghia, Ghia & Shin (1982) solved (2-D streamfunction-vorticity) on $129^2$ or
$257^2$ grids with multigrid-accelerated SOR. Tabulated $u(0.5, y)$ and
$v(x, 0.5)$ are shipped in `ghia::table_re_XXX` as `&[(f64,f64)]` (position,
velocity component) with linear interpolation for intermediate points.

### Tested Invariants

- **Centerline profiles**: $L_2$ error over the tabulated points

  $$\|u_h - u_{\mathrm{Ghia}}\|_{L_2} = \sqrt{\sum_i w_i (u_h(x_i)-u_G(x_i))^2}, \quad w_i = \Delta y_i \tag{1}$$

  Pass criterion: $L_2 < 5\times10^{-3}$ on $129^2$, $< 1\times10^{-3}$ on
  $257^2$.

- **Primary vortex center**: $(x_c, y_c)$ is argmin of streamfunction; reference
  $x_c = 0.6172, y_c = 0.7344$ at $\mathrm{Re}=100$; deviation $< 1\%$.

- **Secondary vortex topology**: number and approximate location of corner
  vortices matches published counts for the given Re.

- **Mass conservation**: $\max|\nabla\cdot\mathbf{u}| < 10^{-8}$ (nondim).

### Convergence Order Claim

Second-order central discretization with Thom wall vorticity: expected $p \approx
2$ spatial. Demonstrated via Richardson extrapolation (Section 6 below).

### API

```rust
use cfd_validation::ghia::{GhiaOracle, Centerline};

let oracle = GhiaOracle::re_400();  // load Re=400 table
let u_prof = solver.centerline_profile(Centerline::VerticalMid);
let err = oracle.l2_error(&u_prof);
assert!(err < 5e-3);
let (xc, yc) = solver.primary_vortex_center();
assert!((xc - 0.5547).abs() < 0.01); // Re=400 reference
```

---

## 2 — Hagen-Poiseuille Pipe Oracle

### Governing Setup

Fully developed laminar flow in a circular pipe radius $R$, uniform pressure
gradient $\partial p/\partial z = -G$, parabolic velocity

$$u_z(r) = \frac{G}{4\mu}(R^2-r^2) = 2U_{\mathrm{avg}}(1-r^2/R^2) \tag{2}$$

with $U_{\mathrm{avg}} = GR^2/(8\mu)$. Friction factor $f = 64/\mathrm{Re}_D$
where $\mathrm{Re}_D = U_{\mathrm{avg}} D / \nu$.

### Tested Invariants

- **Velocity profile**: $\max_r |u_h(r)-u_{\mathrm{exact}}(r)|/U_{\max} <
  10^{-6}$ on a $65^2$ grid with second-order central scheme (smooth profile
  is resolved accurately; error $\sim O(\Delta r^2)$).

- **Pressure drop**: $\Delta p = f (L/D) (\rho U^2/2)$ matches analytical
  $G L$ within $\epsilon = O(\Delta r^2)$.

- **Entrance length**: when the domain length $L \ge 2 L_h$ with
  $L_h = 0.06\,\mathrm{Re}_D D$, the validation slice (at $z = L/2$) velocity
  defect from parabolic is $<0.1\%$.

### Convergence Order Claim

Second-order spatial for 2-D pipe with axis refinement: $p=2$ in $L_\infty$ on
radial coordinate. Demonstrated both in `cfd-validation` example and the 2-D
pipe example in `cfd-2d`.

### API

```rust
use cfd_validation::poiseuille::{PoiseuilleOracle, PipeParams};

let params = PipeParams { radius: 1.0, grad_p: 8.0, mu: 1.0, rho: 1.0 };
let oracle = PoiseuilleOracle::new(params);
let err = oracle.max_abs_error(&solver_profile);
assert!(err < 1e-10);
```

---

## 3 — Blood Poiseuille (Non-Newtonian Extension)

### Governing Setup

Same pipe geometry as (2) but with $\mu = \mu(\dot{\gamma})$ given by the
Carreau model

$$\mu(\dot{\gamma}) = \mu_\infty + (\mu_0-\mu_\infty)[1+(\lambda\dot{\gamma})^2]^{(n-1)/2} \tag{3}$$

with blood parameters $\mu_\infty=0.00345$, $\mu_0=0.056$, $\lambda=3.313$,
$n=0.3568$. Known semi-analytical solution via quadrature for $u(r)$; CFDrs
validates against this closed-form quadrature oracle, not DNS.

### Tested Invariants

- **Newtonian limit recovery**: as $\lambda \to 0$, viscosity $\mu \to \mu_\infty$
  and profile converges to (2) — validated at $\lambda=10^{-6}$.

- **Shear-thinning peak broadening**: flattened core, higher near-wall shear.

- **Fahraeus-Lindqvist adjustment**: apparent viscosity modulation when
  enabled, $O(1)$ effect at $D < 300\,\mu\mathrm{m}$.

---

## 4 — Turbulent Channel Oracle (DNS Reference)

### Governing Setup

Fully developed plane channel, half-height $h$, friction Reynolds
$Re_\tau = u_\tau h/\nu$. Periodic streamwise and spanwise, no-slip at
$y=\pm h$.

Reference data from Moser, Kim & Mansour (1999) at $Re_\tau=180,395,590$:
mean profile $U^+(y^+)$, Reynolds stresses $\langle u'_i u'_j\rangle$, TKE
budget terms.

### Tested Invariants (RANS path)

- **Log-law recovery**: $U^+ = \kappa^{-1}\ln y^+ + B$ with $\kappa\approx
  0.41$, $B\approx 5.2$ in $30<y^+<0.2 Re_\tau$; mean absolute deviation
  $< 3\%$ from Moser tabulated curve.

- **$k$ and $\epsilon$ budgets**: production equals dissipation at peak of
  buffer layer ($y^+ \approx 12$) within $5\%$ for well-resolved LES.

- **Friction coefficient**: $C_f = 2u_\tau^2/U_b^2$ within $2\%$ of DNS.

### Convergence Order Claim

For RANS $k$-$\epsilon$ with high-Re wall functions, first-order wall-adjacent
layer; for wall-integrated $k$-$\omega$ SST, second-order. LES with Smagorinsky
wall-resolved: formal order $2$ but degraded by filter-width commutation error.

---

## 5 — Rayleigh-Plesset Cavitation Oracle

### Governing Setup

Single spherical bubble radius $R(t)$, liquid pressure $p_\infty(t)$, vapour
pressure $p_v$, surface tension $S$:

$$R\ddot{R} + \tfrac32 \dot{R}^2 = \frac{p_v - p_\infty - 2S/R - 4\mu_l \dot{R}/R}{\rho_l} \tag{4}$$

Neglecting $S,\mu_l$ for the growth phase gives Rayleigh collapse time

$$\tau_c = 0.915\, R_0 \sqrt{\rho_l/\Delta p}, \quad \Delta p = p_\infty - p_v \tag{5}$$

### Tested Invariants

- **Collapse time**: computed $\tau_c$ matches (5) within $2\%$ for
  $\Delta p/p_v > 10$.
- **Energy scaling**: collapse pressure peak $\propto R_0 \Delta p^2 / \rho_l$.

---

## 6 — Richardson Extrapolation and GCI

For three grids $h_1 < h_2 < h_3 = r h_2 = r^2 h_1$:

$$p = \frac{\ln\frac{f_2-f_3}{f_1-f_2}}{\ln r}, \quad f_{\mathrm{exact}}\approx f_1 + \frac{f_1-f_2}{r^p-1} \tag{6}$$

Grid convergence index (Roache):

$$\mathrm{GCI}_{12} = \frac{F_s|f_1-f_2|/|f_1|}{r^p-1}, \quad F_s=1.25\text{ (3-grid)} \tag{7}$$

Pass criteria: $p \in [0.9p_{\mathrm{expected}}, 1.1p_{\mathrm{expected}}]$,
$\mathrm{GCI}_{12} + \mathrm{GCI}_{23}$ gives a confidence interval.

### API

```rust
use cfd_validation::richardson::RichardsonExtrapolation;

let extrap = RichardsonExtrapolation::new(vec![f_coarse, f_medium, f_fine], r = 2.0);
println!("observed order p = {}", extrap.observed_order());
println!("extrapolated f_exact = {}", extrap.extrapolated_value());
println!("GCI_12 = {}", extrap.gci());
```

---

## 7 — Comprehensive Validation Orchestration

`comprehensive_validation_suite` in `cfd-validation` runs every oracle above in
sequence and writes a consolidated JSON report:

```json
{
  "ghia_re_100":  { "l2_error": 0.003, "pass": true },
  "ghia_re_400":  { "l2_error": 0.004, "pass": true },
  "poiseuille":   { "max_error": 1e-10, "pass": true },
  "richardson_p": { "observed_order": 1.98, "expected": 2.0, "pass": true },
  "gci":          { "gci_12": 0.02, "pass": true }
}
```

Used both as an example and as a CI gate (`cargo run -p cfd-validation --example
comprehensive_validation_suite` returns exit code $1$ on failure).

---

## Examples

| Example | Validates | Key check |
|---------|-----------|-----------|
| `comprehensive_validation_suite` | All oracles | Exit-code CI gate |
| `richardson_convergence` | Spatial order $p$ via (6) | $p \approx 2$ (2-D stencil) |
| `blood_poiseuille_2d` | 2-D blood Poiseuille vs analytical parabolic / semi-analytical Carreau | $L_\infty$ and $f=64/\mathrm{Re}_D$ |

Run all three with:

```bash
cargo run -p cfd-validation --example comprehensive_validation_suite
cargo run -p cfd-validation --example richardson_convergence
cargo run -p cfd-validation --example blood_poiseuille_2d
```

---

## Further Reading

- Ghia, Ghia & Shin, JCP 48 (1982).
- Roache, J. Fluids Eng. 116(2), 405 (1994) — GCI original paper.
- Moser, Kim & Mansour, Phys. Fluids 11(4), 943 (1999).
- `cfd-validation` source: `crates/cfd-validation/src/`.
- [Core flow benchmarks](core_flows.md) for benchmark physics.
- [Numerics chapter](numerics_and_solvers.md) for Richardson theory.
- [Migration validation](migration_validation.md) for legacy-vs-Atlas parity.
