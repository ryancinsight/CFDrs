# Chapter 3 — Turbulence Models and Cavitation

## Overview

CFDrs models three coupled physical regimes: **turbulence** (RANS / LES
closures), **multiphase flows** (phases with distinct equations of state), and
**cavitation** (liquid–vapour phase transition at low pressure). Turbulence
and cavitation interact strongly — cavitation inception is seeded in
low-pressure turbulent eddies, while cavitation bubble collapse feeds back into
turbulence kinetic energy. This chapter presents the Reynolds-averaged and
large-eddy foundations, derives the key closure models shipped in `cfd-3d`,
introduces the multiphase mixture framework, and gives the Rayleigh-Plesset
cavitation formulation with damage accumulation.

---

## Governing Equations — Full System

The compressible mixture form that underlies all three regimes is

$$\frac{\partial (\alpha_k \rho_k)}{\partial t} + \nabla \cdot (\alpha_k \rho_k \mathbf{u}_k) = \Gamma_k \tag{1}$$

$$\frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) = -\nabla p + \nabla \cdot \boldsymbol{\tau} + \mathbf{f} \tag{2}$$

where $k$ indexes phases, $\alpha_k$ is the phase volume fraction with
$\sum_k \alpha_k = 1$, $\rho = \sum_k \alpha_k \rho_k$ is mixture density,
and $\Gamma_k$ is inter-phase mass transfer (nonzero only under cavitation).
For incompressible single-phase without cavitation, (1) reduces to
$\nabla \cdot \mathbf{u} = 0$.

---

## 1 — Reynolds Averaging and the Closure Problem

### Reynolds Decomposition

Decompose velocity and pressure into mean and fluctuation:

$$\mathbf{u} = \bar{\mathbf{u}} + \mathbf{u}', \qquad p = \bar{p} + p' \tag{3}$$

where the overbar denotes ensemble or time average with
$\overline{\mathbf{u}'} = 0$, and the average commutes with differentiation.

### RANS Equations

Averaging the incompressible N-S equations yields

$$\nabla \cdot \bar{\mathbf{u}} = 0 \tag{4}$$

$$\frac{\partial \bar{\mathbf{u}}}{\partial t} + (\bar{\mathbf{u}} \cdot \nabla)\bar{\mathbf{u}} = -\frac{1}{\rho}\nabla\bar{p} + \nu\nabla^2\bar{\mathbf{u}} - \nabla \cdot \mathbf{R} + \bar{\mathbf{f}} \tag{5}$$

with the Reynolds stress tensor

$$\mathbf{R} = \overline{\mathbf{u}' \otimes \mathbf{u}'} = \left[ \overline{u'_i u'_j} \right] \tag{6}$$

Equations (4)-(5) have $4$ equations and $10$ unknowns (3 mean velocities +
mean pressure + 6 independent Reynolds stress components). This is the
**closure problem**: a model for $\mathbf{R}$ must be supplied.

### Boussinesq Eddy-Viscosity Hypothesis

The dominant closure route approximates $\mathbf{R}$ via a turbulent viscosity
$\nu_t$:

$$\mathbf{R} - \frac{2}{3}k\mathbf{I} = -\nu_t \left(\nabla\bar{\mathbf{u}} + (\nabla\bar{\mathbf{u}})^T\right) \tag{7}$$

where $k = \tfrac{1}{2}\mathrm{tr}(\mathbf{R})$ is turbulence kinetic energy
and $k$ enters the trace. Modellers solve transport equations for $k$ and a
second scale variable ($\epsilon$ or $\omega$) to determine $\nu_t$.

---

## 2 — Closure Models Shipped in CFDrs

### k-$\epsilon$ (Standard and RNG)

The standard $k$-$\epsilon$ closure solves

$$\frac{\partial k}{\partial t} + \bar{\mathbf{u}}\cdot\nabla k = P_k - \epsilon + \nabla\cdot\left[\left(\nu + \frac{\nu_t}{\sigma_k}\right)\nabla k\right] \tag{8}$$

$$\frac{\partial \epsilon}{\partial t} + \bar{\mathbf{u}}\cdot\nabla\epsilon = C_{\epsilon 1}\frac{\epsilon}{k}P_k - C_{\epsilon 2}\frac{\epsilon^2}{k} + \nabla\cdot\left[\left(\nu + \frac{\nu_t}{\sigma_\epsilon}\right)\nabla\epsilon\right] \tag{9}$$

with production $P_k = \nu_t S_{ij}S_{ij}$ where
$S_{ij} = \tfrac{1}{2}(\partial_i \bar{u}_j + \partial_j \bar{u}_i)$ is the mean
strain rate, and eddy viscosity

$$\nu_t = C_\mu \frac{k^2}{\epsilon} \tag{10}$$

Constants (standard model, Launder & Spalding 1974):
$C_\mu = 0.09$, $C_{\epsilon 1} = 1.44$, $C_{\epsilon 2} = 1.92$,
$\sigma_k = 1.0$, $\sigma_\epsilon = 1.3$.

CFDrs implements both the standard and RNG variants (Yakhot & Orszag 1986) with
a modified destruction coefficient. The closure lives at
`cfd-3d::turbulence::KEpsilon { variant: Standard | RNG }`.

Wall treatment: standard $k$-$\epsilon$ is a high-Re model — a wall function
bridges the viscous sublayer:

$$U^+ = \frac{1}{\kappa}\ln(E y^+), \quad k = \frac{u_\tau^2}{\sqrt{C_\mu}}, \quad \epsilon = \frac{u_\tau^3}{\kappa y} \tag{11}$$

with $\kappa = 0.41$, $E = 9.8$ (smooth wall).

### k-$\omega$ SST (Menter 1994)

$k$-$\omega$ transport:

$$\frac{\partial k}{\partial t} + \bar{\mathbf{u}}\cdot\nabla k = P_k - \beta^* k\omega + \nabla\cdot[(\nu + \sigma_k \nu_t)\nabla k] \tag{12}$$

$$\frac{\partial \omega}{\partial t} + \bar{\mathbf{u}}\cdot\nabla\omega = \alpha \frac{\omega}{k}P_k - \beta\omega^2 + \nabla\cdot[(\nu + \sigma_\omega \nu_t)\nabla\omega] + 2(1-F_1)\sigma_{\omega 2}\frac{\nabla k \cdot \nabla\omega}{\omega} \tag{13}$$

Eddy viscosity with SST limiter:

$$\nu_t = \frac{a_1 k}{\max(a_1\omega, S F_2)} \tag{14}$$

where $F_1$ blends between near-wall $k$-$\omega$ and far-wall $k$-$\epsilon$
formulation (transformed), and $F_2$ is the shear-stress transport blending
function. The SST modification guarantees $\nu_t$ does not exceed
$k / \omega$ times the strain-rate ratio, preventing overshoot in adverse
pressure-gradient boundary layers — critical for venturi pressure recovery
validation.

CFDrs constants follow Menter's 2003 update at
`cfd-3d::turbulence::KOmegaSST`.

### Smagorinsky LES

Filtering the N-S equations with a spatial filter of width $\Delta$ yields the
subgrid-scale (SGS) stress $\tau_{ij}^{\mathrm{SGS}} = \overline{u_i u_j} -
\bar{u}_i \bar{u}_j$. The Smagorinsky (1963) model closes it as

$$\tau_{ij}^{\mathrm{SGS}} - \tfrac{1}{3}\tau_{kk}^{\mathrm{SGS}}\delta_{ij} = -2 \nu_{\mathrm{sgs}} \bar{S}_{ij} \tag{15}$$

$$\nu_{\mathrm{sgs}} = (C_s \Delta)^2 |\bar{S}|, \quad |\bar{S}| = \sqrt{2\bar{S}_{ij}\bar{S}_{ij}} \tag{16}$$

with $C_s = 0.16$ (isotropic turbulence, Lilly 1967) or dynamically computed via
Germano identity. The filter width $\Delta = (\Delta_x \Delta_y \Delta_z)^{1/3}$
on structured grids.

Van Driest damping near walls:

$$\nu_{\mathrm{sgs}}^{\mathrm{damped}} = \nu_{\mathrm{sgs}}\left[1 - \exp(-y^+ / 25)\right]^2 \tag{17}$$

CFDrs sends LES through `cfd-3d::turbulence::Smagorinsky { cs: f64, van_driest:
bool }`. The discretization trait returns $\nu_t = \nu_{\mathrm{sgs}}$ per cell
that the momentum solver adds to molecular viscosity.

```rust
pub trait TurbulenceModel {
    type Field;
    fn step<F: FloatElement>(&mut self, u: &Array3<F>, v: &Array3<F>, w: &Array3<F>,
                             nu_t: &mut Array3<F>, dt: F);
    fn label(&self) -> &'static str;
}
```

---

## 3 — Multiphase Coupling

For immiscible liquid-gas (or liquid-solid) flows CFDrs uses a
**volume-of-fluid** (VOF) style mixture: each cell carries $\alpha \in [0,1]$,
liquid fraction in the liquid-gas case, with mixture properties

$$\rho_m = \alpha \rho_l + (1-\alpha)\rho_g, \quad \mu_m = \alpha \mu_l + (1-\alpha)\mu_g \tag{18}$$

and phase advection

$$\frac{\partial \alpha}{\partial t} + \mathbf{u}\cdot\nabla\alpha + \nabla\cdot\mathbf{u}_c\alpha(1-\alpha) = 0 \tag{19}$$

where the compressive term $\mathbf{u}_c \propto |\mathbf{u}| \hat{\mathbf{n}}$
sharpens the interface (Weller-VOF method).

The `MomentumCoupling` enum selects co-located vs staggered velocity/pressure
coupling. SSE-style inter-phase momentum exchange:

$$\mathbf{F}_D = C_D \rho_c |\mathbf{u}_c - \mathbf{u}_d| (\mathbf{u}_c - \mathbf{u}_d) \tag{20}$$

with $C_D$ from Schiller-Naumann drag correlation.

---

## 4 — Cavitation Physics — Rayleigh-Plesset and Damage

### Cavitation Inception

Cavitation occurs when local static pressure falls below vapour pressure $p_v$:

$$p(\mathbf{x}, t) < p_v(T) \tag{21}$$

The nondimensional cavitation number

$$\sigma = \frac{p_\infty - p_v}{\tfrac{1}{2}\rho U_\infty^2} \tag{22}$$

characterizes cavitation susceptibility. Low $\sigma$ means higher cavitation
risk. In venturi geometries $\sigma$ varies along the axis — CFDrs reports the
minimum $\sigma$ as the cavitation risk indicator.

### Rayleigh-Plesset Equation

Bubble radius dynamics follow Rayleigh-Plesset:

$$R\ddot{R} + \frac{3}{2}\dot{R}^2 = \frac{1}{\rho_l}\left(p_v - p_\infty - \frac{2S}{R} - \frac{4\mu_l \dot{R}}{R}\right) \tag{23}$$

where $S$ is surface tension, $\mu_l$ liquid viscosity. Neglecting surface
tension and viscous damping for the growth phase:

$$R\ddot{R} + \tfrac{3}{2}\dot{R}^2 = \frac{p_v - p_\infty}{\rho_l} \tag{24}$$

Rayleigh collapse time for a void collapsing from $R_0$ under constant
overpressure $\Delta p = p_\infty - p_v$ is

$$\tau_c = 0.915\, R_0 \sqrt{\frac{\rho_l}{\Delta p}} \tag{25}$$

CFDrs uses (23)-(25) for the single-bubble growth-collapse cycle rate.

### Damage Index

Material damage from cavitation follows a Wöhler-like accumulation. CFDrs
defines per-cell damage rate

$$\dot{D} = \mathbf{1}_{p < p_v} \cdot \frac{|\dot{R}|_{\mathrm{collapse}}}{R_{\mathrm{char}}} \cdot f_{\mathrm{impact}}(p) \tag{26}$$

integrated as $D(t) = \int_0^t \dot{D}\, dt$ with units kg/m²/s after
dimensional recovery. The damage hook `cfd-3d::posterior::damage` reports it.

```rust
pub trait Cavitation {
    fn inception_pressure<F: FloatElement>(&self) -> F;
    fn collapse_rate<F: FloatElement>(&self, p: F) -> F;
}
```

### Integration Contract

All three regimes compose through the same `Integrate` trait:

```rust
pub trait Integrate {
    fn step<F: FloatElement>(&mut self) -> Result<(), CfdError>;
}
```

The implementation multiplexes the closure choice at the call site. Regime
switching — RANS to LES or cavitation on/off — is a configuration change, not
a code rewrite.

---

## Validation Oracles

| Oracle | Source | Quantity checked |
|---|---|---|
| Turbulent channel RANS | Moser et al. $Re_\tau=180,395,590$ DNS | $U^+(y^+)$, wall shear, $\langle u_i'u_j'\rangle$ |
| Venturi cavitation inception | Franc & Michel (2005) experiments, Figure 4.12 | Inception $\sigma_i$ vs geometry |
| Rayleigh-Plesset single bubble | Rayleigh (1917) collapse time (25) | $\tau_c$ vs $\Delta p$, $R_0$ |
| VOF interface sharpness | Hirt & Nichols (1981) slotted-disk advection | $\|\nabla\alpha\|$ conservation |

---

## API Mapping

| Concern | Crate / path |
|---|---|
| Turbulence closures | `cfd-3d::turbulence::{KEpsilon, KOmegaSST, Smagorinsky}` |
| VOF / mixture | `cfd-3d::multiphase::{VofField, MomentumCoupling}` |
| Cavitation model | `cfd-3d::cavitation::{RayleighPlesset, EulerianEulerian, Cavitation}` |
| Damage posterior | `cfd-3d::posterior::damage` |
| Integration contract | `cfd-3d::integrate::Integrate` |

---

## Examples Referenced by This Chapter

- [Example: turbulence_models_demo](examples/turbulence_models_demo.md) — three closures compared on a turbulent channel; $k$ and $\nu_t$ profiles plotted.
- [Example: turbulence_validation_demo](examples/turbulence_validation_demo.md) — closure against Moser DNS at $Re_\tau=180$.
- [Example: turbulence_momentum_integration_demo](examples/turbulence_momentum_integration_demo.md) — closure-integrated momentum solve showing log-law recovery.
- [Example: venturi_cavitation](examples/venturi_cavitation.md) — venturi cavitation inception and collapse; $\sigma$ along axis.
- [Example: simple_cavitation](examples/simple_cavitation.md) — single-bubble Rayleigh-Plesset in a constrained channel.
- [Example: cavitation_damage_simulation](examples/cavitation_damage_simulation.md) — long-time damage integration under turbulent inflow.
- [Example: 1d_venturi_cavitation](examples/1d_venturi_cavitation.md) — 1-D reduced cavitation model with $\sigma$-based regime classification.

## Further Reading

- Pope, *Turbulent Flows* (2000) — Chapters 4–6 for Reynolds averaging; Chapter 10 for $k$-$\epsilon$.
- Menter, AIAA J. 32(8), 1598 (1994) for $k$-$\omega$ SST.
- Smagorinsky, MWR 91(3), 99 (1963); Lilly, Proc. IBM Sci. Symp. (1967).
- Rayleigh, Phil. Mag. 34, 94 (1917) — original collapse solution.
- Franc & Michel, *Fundamentals of Cavitation*, Kluwer (2005).
- `cfd-3d` source: `crates/cfd-3d/src/{turbulence,multiphase,cavitation}.rs`.
- [Core flow benchmarks](core_flows.md) for DNS wall-scaling validation context.
- [Numerics and Solvers](numerics_and_solvers.md) for $\nu_t$-implicit treatment.
