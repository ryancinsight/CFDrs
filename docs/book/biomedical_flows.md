# Chapter 4 — Blood Flow and Rheology Workflows

## Overview

CFDrs supports a focused subset of physiologically relevant geometries —
**vascular bifurcations**, **venturi / stenosis regions**, **microfluidic
screens**, and **TPMS (triply-periodic minimal-surface) scaffolds**. Each
geometry couples a `cfd-schematics` CSG primitive to a low-dimensional CFD
solve (1-D, 2-D, or 3-D), with **non-Newtonian rheology** for blood. The
chapter covers the fluid-mechanics scaling (Womersley number, Murray's law),
blood rheology hierarchy (Newtonian through Carreau-Yasuda and Casson),
microvascular effects (Fahraeus-Lindqvist), wall shear stress, hemolysis
assessment, and the `cfd-1d` screening pipelines that walk whole-device chips
and flag regulatory violations.

---

## 1 — Governing Equations and Physiological Scaling

Blood is incompressible, so the foundation is still

$$\nabla \cdot \mathbf{u} = 0, \qquad \rho\left(\frac{\partial\mathbf{u}}{\partial t} + \mathbf{u}\cdot\nabla\mathbf{u}\right) = -\nabla p + \nabla\cdot\boldsymbol{\tau} + \mathbf{f} \tag{1}$$

with non-Newtonian extra stress $\boldsymbol{\tau} = 2\mu(\dot{\gamma})\mathbf{S}$
where $\dot{\gamma} = \sqrt{2\,\mathbf{S}:\mathbf{S}}$ is the scalar shear
rate, $\mathbf{S} = \tfrac12(\nabla\mathbf{u} + \nabla\mathbf{u}^T)$, and
$\mu(\dot{\gamma})$ is a viscosity function introduced in Section 2.

### Dimensionless Groups in Biomedical Flows

#### Womersley Number

For pulsatile flow driven at angular frequency $\omega = 2\pi f$, the
Womersley number is

$$\alpha = R\sqrt{\frac{\omega}{\nu}} = R\sqrt{\frac{\omega \rho}{\mu}} \tag{2}$$

where $R$ is vessel radius, $\nu = \mu/\rho$ kinematic viscosity, $f$ heart
rate (or pump pulsatility). Interpretation:

- $\alpha \ll 1$: viscous-dominated — parabolic profile in phase with
  instantaneous pressure gradient (quasi-steady Hagen-Poiseuille at each
  instant).
- $\alpha \sim 1$–$10$: intermediate — phase lag between pressure and flow,
  annular peak shifts toward the wall.
- $\alpha \gg 10$: inertia-dominated — velocity profile nearly flat in the
  core with thin Stokes layers $\delta = \sqrt{2\nu/\omega}$ at walls.

In systemic circulation $\alpha \approx 20$ in aorta, $\approx 2$–$5$ in
medium arteries, $<1$ in arterioles and microfluidics. CFDrs microfluidic chip
examples operate at $\alpha \approx 0.1$–$1$, so quasi-steady viscous models
are valid without Womersley correction. Venturi stenosis models at
$\alpha \approx 3$–$6$ include the full analytical Womersley solution for the
inlet BC.

#### Analytical Womersley Solution for a Straight Vessel

For fully developed axisymmetric pulsatile flow with pressure gradient
$\partial p/\partial z = \mathrm{Re}[\hat{G} e^{i\omega t}]$, the velocity is

$$u_z(r, t) = \mathrm{Re}\left[\frac{\hat{G}}{i\omega\rho}\left(1 - \frac{J_0(i^{3/2}\alpha r/R)}{J_0(i^{3/2}\alpha)}\right) e^{i\omega t}\right] \tag{3}$$

where $J_0$ is the Bessel function of the first kind, order 0, and
$i^{3/2}\alpha = (1+i)\alpha/\sqrt{2}$. Equation (3) is used as:

- An **inlet boundary condition** for pulsatile venturi/bifurcation solves.
- A **validation oracle** for `cfd-1d::WomersleyInlet` correctness.
- A **posterior formula** for wall shear stress in rigid straight vessels.

CFDrs implements (3) in `cfd-1d::womersley::WomersleyProfile` with truncation of
the Bessel series to machine epsilon; accuracy is verified against the exact
expansion at $\alpha \in [0.1, 20]$.

#### Reynolds Number in Biomedical Context

Including non-Newtonian effects, CFDrs defines an *effective* Reynolds number

$$\mathrm{Re}_{\mathrm{eff}} = \frac{\rho U D}{\mu_{\mathrm{eff}}}, \quad \mu_{\mathrm{eff}} = \mu(\dot{\gamma}=U/D) \tag{4}$$

reported by every solve so that Newtonian thresholds (e.g. $\mathrm{Re} < 2100$
for laminar pipe) can be re-interpreted. Typical values: $Re \approx 1$ in
microfluidics, $\approx 100$–$2000$ in medium arteries and venturis.

---

## 2 — Blood Rheology Hierarchy

Blood shear-thins: red blood cells aggregate at low shear and disaggregate at
high shear. CFDrs provides the following hierarchy, ordered by complexity:

### Newtonian (constant $\mu$)

$$\mu = \mu_{\infty} \approx 3.45\times 10^{-3}\,\mathrm{Pa\cdot s} \tag{5}$$

For CFDrs tests that need a linear baseline or where shear rates exceed
$500\,\mathrm{s}^{-1}$ uniformly.

### Power-Law

$$\mu(\dot{\gamma}) = K \dot{\gamma}^{n-1} \tag{6}$$

$K = 0.035$ Pa·s$^n$, $n = 0.6$ (Chien et al.). Cheap, widely used, but
predicts unbounded viscosity at $\dot{\gamma} \to 0$.

### Carreau and Carreau-Yasuda

$$\mu(\dot{\gamma}) = \mu_\infty + (\mu_0 - \mu_\infty)\left[1 + (\lambda\dot{\gamma})^2\right]^{\frac{n-1}{2}} \tag{7}$$

Carreau constants (Cho & Kensey): $\mu_\infty = 0.00345$ Pa·s, $\mu_0 = 0.056$
Pa·s, $\lambda = 3.313$ s, $n = 0.3568$. Carreau-Yasuda adds a second exponent

$$\mu(\dot{\gamma}) = \mu_\infty + (\mu_0 - \mu_\infty)\left[1 + (\lambda\dot{\gamma})^a\right]^{\frac{n-1}{a}} \tag{8}$$

with $a = 2$ often, allowing independent control of low-shear plateau broadness.

### Casson

$$\sqrt{\tau} = \sqrt{\tau_y} + \sqrt{\mu_c \dot{\gamma}} \tag{9}$$

where $\tau$ is shear stress magnitude, $\tau_y$ yield stress ($\approx 0.004$
Pa for blood), $\mu_c$ Casson viscosity. In apparent-viscosity form

$$\mu(\dot{\gamma}) = \frac{(\sqrt{\tau_y} + \sqrt{\mu_c\dot{\gamma}})^2}{\dot{\gamma}} \tag{10}$$

Casson introduces a true yield stress — relevant for near-stagnation zones
in microfluidic reservoirs where $\dot{\gamma} \to 0$ (e.g. expansion chambers).

### API

```rust
use cfd_1d::blood_rheology::Rheology;

let rheo = Rheology::carreau {          // μ_0, μ_∞, λ, n
    mu_inf: 0.00345,
    mu_0:   0.056,
    lambda: 3.313,
    n:      0.3568,
};

let rheo_cy = Rheology::carreau_yasuda {
    mu_inf: 0.00345, mu_0: 0.056, lambda: 3.313, n: 0.3568, a: 2.0,
};

let rheo_casson = Rheology::casson {
    tau_y: 0.004,  // Pa
    mu_c: 0.00345,
};
```

The solver consumes `&dyn RheologyModel` per timestep — switching models is a
single field change; no template re-instantiation needed since the bound is
dynamic at the screening layer and Atlas-monomorphized in the kernel layer.

---

## 3 — Microvascular Effects

### Fahraeus-Lindqvist Effect

Apparent blood viscosity drops as tube diameter decreases below about $300\,\mu\mathrm{m}$,
because red blood cells migrate toward the axis, leaving a cell-depleted layer
at the wall. CFDrs encodes the Pries-Secomb empirical correlation:

$$\mu_{\mathrm{rel}}(D) = 1 + (\mu_{0.45} - 1)\frac{(1-H_D)^C - 1}{(1-0.45)^C - 1} \tag{11}$$

$$\mu_{0.45} = 220 e^{-1.3 D} + 3.2 - 2.44 e^{-0.06 D^{0.645}} \tag{12}$$

where $D$ is tube diameter in $\mu\mathrm{m}$, $H_D$ discharge haematocrit,
$\mu_{0.45}$ the relative apparent viscosity at $H_D = 0.45$, and
$C = 0.07 e^{-0.06 D} + 0.8 -0.017(e^{-0.15 D}) + 6.0 e^{-0.05 D}$ (Pries et al.
1992). Implementation lives in `cfd-1d::microvascular::FahraeusCorrection`.

For microfluidic chip designs with channel widths $100$–$250\,\mu\mathrm{m}$, FL
correction can shift apparent viscosity by factor $1.5$–$2$, directly impacting
pressure drop and hemolysis predictions. CFDrs chip screening runs both with
and without FL correction to quantify sensitivity.

### Plasma Skimming and Phase Separation at Bifurcations

At bifurcations, uneven haematocrit splitting (Zweifach-Fung effect) is modelled
via Pries' bifurcation law: the higher-flow child receives disproportionately
high $H_D$. This feeds into Murray's law adaptation of effective hydraulic
conductance per branch.

---

## 4 — Vessel Sizing — Murray's Law

Murray's law gives the optimal radius ratio at a vascular bifurcation that
minimizes work (viscous dissipation + metabolic cost of blood volume):

$$r_p^3 = r_a^3 + r_b^3 \tag{13}$$

where $r_p$ is parent radius, $r_a, r_b$ child radii. More generally for a
bifurcation with $n$ children,

$$r_p^g = \sum_i r_i^g \tag{14}$$

with exponent $g = 3$ classically (laminar flow), $g = 2.33$ for turbulent
flow, and $g \approx 2.7$ empirically for systemic arteries. CFDrs bifurcation
builder uses $g = 3$ by default with $g$ configurable:

```rust
use cfd_schematics::bifurcation::BifurcationBuilder;
use leto::Point3;

let bif = BifurcationBuilder::new()
    .parent(Point3::new(0.0, 0.0, 0.0), 2.5)            // 5 mm diameter, mm units
    .child_a(Point3::new(0.0, 4.0, 0.0), 1.8, 0.35)    // 3.6 mm, 35° angle
    .child_b(Point3::new(0.0, 4.0, 0.0), 1.8, -0.35)   // mirror child
    .stenosed_a(0.55)                                    // 55% area stenosis
    .murray_exponent(3.0)
    .build();
```

Flow division follows from continuity and the principle of equal wall shear
stress in an optimal tree:

$$Q_a / Q_b = r_a^3 / r_b^3 \tag{15}$$

CFDrs validates (13)-(15) against the reference phantom Poiseuille pressure drops
— deviation above $1\%$ in non-stenosed optimal bifurcations signals a mesh or
wiring error.

### Stenosis

Stenosed area fraction $S = 1 - A_{\mathrm{stenosis}}/A_{\mathrm{normal}}$.
Clinical grading: mild $S < 50\%$, moderate $50\% \le S < 70\%$, severe $S \ge
70\%$. The pressure drop across an axisymmetric stenosis follows the empirical
Young-Tsai relation $\Delta p = K_v \mu Q/A_0^2 + K_t \rho(Q/A_0)^2$ where $K_v$
and $K_t$ are geometry coefficients computed by CFDrs from the venturi builder.

---

## 5 — Wall Shear Stress (WSS) and Mechanobiology

Wall shear stress magnitude

$$\tau_w = \mu_{\mathrm{eff}} \left|\frac{\partial \mathbf{u}_t}{\partial n}\right|_{\mathrm{wall}} \tag{16}$$

drives endothelial mechanotransduction. CFDrs computes it per wall-adjacent
cell via one-sided differencing using the effective non-Newtonian viscosity at
the wall.

Time-averaged WSS (TAWSS) and oscillatory shear index (OSI) for pulsatile runs:

$$\mathrm{TAWSS} = \frac{1}{T}\int_0^T |\tau_w(t)| dt, \quad \mathrm{OSI} = \frac{1}{2}\left(1 - \frac{|\int_0^T \tau_w(t) dt|}{\int_0^T |\tau_w(t)| dt}\right) \tag{17}$$

OSI ranges $0$ (unidirectional) to $0.5$ (pure oscillation). Low TAWSS with high
OSI correlates with atherosclerotic plaque development — CFDrs 1-D posteriors
flag regions with TAWSS $< 0.5$ Pa and OSI $> 0.15$.

Internal elastic lamina stress and longitudinal WSS gradients (WSSG) are also
tracked for aneurysm growth analysis.

---

## 6 — Hemolysis Assessment

### Shear-Stress and Exposure-Time Models

Hemolysis risk is quantified via a power-law model

$$H = C\, \tau^b\, t^a \tag{18}$$

where $H = \Delta Hb / Hb$ is normalized plasma free haemoglobin, $\tau$ is
scalar shear stress (often von Mises or Tresca equivalent), $t$ is exposure
time along a fluid particle trajectory, and $(C, a, b)$ are empirically fitted
constants. Canonical constants (Giersiepen-Wurzinger 1990): $C = 3.62\times
10^{-5}$ (SI: Pa, s), $a = 0.785$, $b = 2.416$.

Alternative calibration by Zhang et al. (2011): $C = 1.228\times10^{-5}$,
$a = 0.6606$, $b = 1.9918$.

### FDA Shear Thresholds

The FDA critical path for blood-contacting devices specifies screening against
peak scalar shear stress:

- **Design limit**: $\tau_{\max} < 150$ Pa (arterial devices), $< 50$ Pa for
  sensitive regions.
- **FDA nominal hemolysis threshold**: sustained shear $> 150$–$400$ Pa for
  exposure times $> 10$ ms is hemolytic.

CFDrs FDA screening walks the chip graph:

```rust
use cfd_1d::screening::shear_violations;

let report = shear_violations(&state, &geometries)?;
for branch in report.violations {
    eprintln!("branch {}: peak shear {} Pa (limit 9.0)", branch.id, branch.peak);
}
```

Report aggregates per-branch and per-particle-line hemolysis indices with the
worst-case $H$ over all integrated trajectories.

### Exposure-Time and Lagrangian Tracking

Because (18) depends nonlinearly on $t$, CFDrs does Lagrangian integration along
streamlines (1-D network) or particle paths (2-D/3-D) accumulating
$H_{\mathrm{particle}} = \int_{\mathrm{path}} C \dot{\gamma}^b dt^a $. The 1-D
implementation uses per-edge mean shear velocity and per-edge residence time,
avoiding explicit particle pushing; the 2-D/3-D implementation advects tracer
particles with adaptive Runge-Kutta 4.

### API Mapping

| Concern | Crate / path |
|---|---|
| Bifurcation builder | `cfd-schematics::bifurcation::BifurcationBuilder` |
| Serpentine builder | `cfd-schematics::serpentine::SerpentineBuilder` |
| Rheology | `cfd-1d::blood_rheology::{Rheology, RheologyModel}` |
| Womersley inlet | `cfd-1d::womersley::WomersleyProfile` |
| FL correction | `cfd-1d::microvascular::FahraeusCorrection` |
| FDA screening | `cfd-1d::screening::{shear_violations, HemolysisReport}` |
| Hemolysis model | `cfd-1d::hemolysis::{GiersiepenWurzinger, Zhang}` |
| WSS posterior | `cfd-1d::posterior::{wall_shear_stress, tawss_osi}` |

---

## 7 — Microfluidic Mixing Screens and TPMS Scaffolds

### Serpentine Mixers

CFDrs handles a serpentine chip as a chain of bifurcations with Dean number

$$\mathrm{De} = \mathrm{Re}\sqrt{D_h / R_c} \tag{19}$$

where $R_c$ is the radius of curvature of the U-turn. De $> 100$ triggers Dean
vortex pairs that enhance transverse mixing; mixing efficiency is quantified as
$M = 1 - \sigma_c / \max \sigma_c$ where $\sigma_c$ is concentration variance.

### TPMS Scaffolds

Triply-periodic minimal surfaces (Schwarz P, Gyroid, Diamond) are implicit
$SDF < 0$ solids tiled for `cfd-3d` FEM. Porosity $\phi$ and specific surface
area $S_v$ follow from the implicit-function offset. CFDrs Gyroid TPMS for tissue
perfusion is built as

$$F(x,y,z) = \sin x \cos y + \sin y \cos z + \sin z \cos x - t = 0 \tag{20}$$

with $t \in [-1,1]$ controlling porosity.

---

## Examples Referenced by This Chapter

- [Example: blood_flow_1d_validation](examples/blood_flow_1d_validation.md) — 1-D reduced model against reference phantom; Womersley vs steady.
- [Example: bifurcation_2d_blood_validation](examples/bifurcation_2d_blood_validation.md) — 2-D bifurcation with Carreau rheology; Murray validation.
- [Example: blood_rheology_models](examples/blood_rheology_models.md) — three rheologies compared on patient-specific geometry; FL correction sensitivity.
- [Example: venturi_blood_flow_validation](examples/venturi_blood_flow_validation.md) — venturi stenosis Newtonian vs Carreau.
- [Example: serpentine_mixing_comprehensive](examples/serpentine_mixing_comprehensive.md) — serpentine mixer design space 1-D to chip-scale; Dean number sweep.
- [Example: cross_fidelity_branching](examples/cross_fidelity_branching.md) — concurrent 1-D/2-D/3-D bifurcation comparison.

## Further Reading

- Cho & Kensey, *Biorheology* 28, 241 (1991) — Carreau-Yasuda parameters for blood.
- Pries, Neuhaus & Gaehtgens, *Am. J. Physiol.* 264, H1170 (1993) — Fahraeus-Lindqvist correlation.
- Murray, *PNAS* 12, 207 (1926) — original optimal-bifurcation derivation.
- Giersiepen et al., *Int. J. Artif. Organs* 13(5), 300 (1990) — hemolysis power-law.
- Womersley, *Philos. Mag.* 46(8), 199 (1955) — pulsatile flow theory.
- `cfd-1d` and `cfd-schematics` source: `crates/cfd-1d/src/{rheology,womersley,microvascular,screening,hemolysis}.rs`.
