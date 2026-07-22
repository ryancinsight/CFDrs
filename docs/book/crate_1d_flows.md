# Chapter 20 — 1-D Biomedical Flows

## Overview

Reduced-order 1-D haemodynamic examples live in `cfd-1d`. They cover the core
biomedical-flow library: bifurcations, microfluidics, FDA device screening,
TPMS scaffolds, cavitation analysis, and hemolysis assessment. The 1-D approach
reduces per-edge 3-D Navier-Stokes to a 1-D network model (pressure $p_j(t)$ and
flowrate $Q_j(t)$ per edge) with Poiseuille-based hydraulic resistance and
Womersley-based inertance, coupled with non-Newtonian rheology correction and
Fahraeus-Lindqvist microvascular effects. A full medical device chip with
$>100$ branches solves in milliseconds at 1-D, enabling rapid design-space
screening before any 2-D or 3-D refinement.

This chapter documents the eight canonical 1-D examples shipped as cargo
examples in `cfd-1d`, their governing equations, validation oracles, tested
invariants, convergence analysis, and Atlas migration stage.

---

## Governing Equations — 1-D Network Model

For each edge $j$ of the vascular/microfluidic graph $G = (V, E)$,

$$\frac{\partial A_j}{\partial t} + \frac{\partial Q_j}{\partial z} = 0 \tag{1}$$

$$\frac{\partial Q_j}{\partial t} + \frac{\partial}{\partial z}\left(\frac{Q_j^2}{A_j}\right) + \frac{A_j}{\rho}\frac{\partial p_j}{\partial z} = -f_j + g_j \tag{2}$$

where $A_j(z,t)$ is cross-sectional area, $Q_j(z,t)$ volumetric flowrate,
$p_j(z,t)$ pressure, $f_j$ friction, $g_j$ body force. For rigid walls
(microfluidics, TPMS scaffolds) $A_j = const$, so (1) reduces to $\partial_z
Q_j = 0$ per edge and (2) becomes the algebraic $\Delta p_j = R_j Q_j + L_j
\dot{Q}_j$.

Constitutive closure for area:

$$p_j(A_j) = p_{j,0} + \beta_j(\sqrt{A_j} - \sqrt{A_{j,0}}) + \zeta_j(\partial_t A_j) \tag{3}$$

where for rigid walls $\beta_j \to \infty$ (or equivalently $C_j = dA_j/dp_j =
0$) so (3) is replaced by fixed $A_j$. For compliant vessels $\beta_j =
\sqrt{\pi} h_j E_j / (1-\nu_j^2)A_{j,0}$ with wall thickness $h_j$, Young
modulus $E_j$, Poisson ratio $\nu_j$.

Friction via Poiseuille:

$$f_j = K_R \frac{\mu_{\mathrm{eff}}}{\rho}\frac{Q_j}{A_j} \tag{4}$$

with $K_R = 22$ for 1-D consistent with Smith et al. (2002). Effective viscosity:

$$\mu_{\mathrm{eff}} = \mu(\dot{\gamma}_{\mathrm{char}})\times F_{\mathrm{FL}}(D_j, H_D) \tag{5}$$

where $\mu(\cdot)$ is non-Newtonian (Carreau etc., see
[biomedical chapter](biomedical_flows.md)) and $F_{\mathrm{FL}}$ is the
Fahraeus-Lindqvist correction.

Kirchhoff junction coupling: per vertex $v\in V$, mass conservation
$\sum_{j\in\mathrm{in}(v)} Q_j = \sum_{j\in\mathrm{out}(v)} Q_j$ and pressure
continuity $p_{j\in\mathrm{in}(v)} = p_{k\in\mathrm{out}(v)}$ (Bernoulli-based
total-pressure junction with loss coefficient $K_{\mathrm{loss}}$ for
stenosis / sharp-bend vertices).

### Time Integration

For unsteady Womersley inlets, backward Euler ($p=1$ A-stable) or BDF2:

$$\frac{3 Q_j^{n+1} - 4Q_j^n + Q_j^{n-1}}{2\Delta t} + R_j Q_j^{n+1} = \Delta p_j^{n+1}/L_j \tag{6}$$

Steady state solved as a sparse linear graph Laplacian $L_G \mathbf{p} =
\mathbf{b}$ via `cfd-math::krylov::BiCgStab`.

---

## 1 — `blood_bifurcation` — Murray's Law + Casson at Symmetric Bifurcation

### Setup

One parent branch → two symmetric children with Murray-optimal radii
$r_p^3 = r_a^3+r_b^3$, $r_a=r_b$. Casson yield-stress rheology. Inlet mean
velocity $U = 0.2\,\mathrm{m/s}$, parent diameter $5\,\mathrm{mm}$.

Discretized with `cfd-1d::network::NetworkSolver` from `NetworkBlueprint`
converted from `BifurcationBuilder`.

### Tested Invariants

- **Murray flow split**: symmetric bifurcation $Q_a/Q_b = 1$ within machine
  epsilon (symmetry of $R_j$ values).

- **Casson pressure**: total pressure drop scales as $\Delta p = R(\mu_{\mathrm{eff}})\,Q$
  where $R \propto \mu_{\mathrm{eff}}/D^4$; Casson pressures are $\sim 20\%$
  higher than Newtonian at $\dot{\gamma} \sim 50\,\mathrm{s}^{-1}$ due to yield
  stress contribution.

- **WSS flag**: wall shear stress $\tau_w > 0.5$ Pa everywhere (non-stenosed
  healthy bifurcation); flagged by `cfd-1d::posterior::wall_shear_stress`.

### Convergence

Edge count per branch $N_j \approx 20$ sufficient: grid-independent to
$<10^{-6}$ fractional error on $\Delta p$ because the 1-D system is one
spatial degree per edge (parabolic profile integrated exactly via resistance
formula; only junction loss $K_{\mathrm{loss}}$ requires resolution).

---

## 2 — `microfluidic_chip` — Pressure-Driven Multi-Channel Chip

### Setup

Serpentine / parallel network of $N \sim 20$–$100$ channels on a chip, average
width $100$–$250\,\mu\mathrm{m}$. Inlet pressure $30\,\mathrm{kPa}$, outlet
atmosphere; $\mathrm{Re} \sim 1$–$10$ (Stokes regime).

### Tested Invariants

- **Perfusion uniformity**: channel-to-channel flowrange ratio
  $\max Q_j / \min Q_j < 3$ for well-designed manifolds (geometry
  scoring hook in `cfd-schematics::serpentine`).

- **Pressure drop**: analytical $f = 64/\mathrm{Re}_D$ for each straight segment
  plus minor-loss $K$-values for U-turns (Idelchik coefficients).

- **Fahraeus sensitivity**: pressure drop sensitivity $S_{\mathrm{FL}} =
  \Delta p_{\mathrm{with\,FL}}/\Delta p_{\mathrm{without\,FL}}$ reported in
  screen for $D_h=100\,\mu\mathrm{m}$ channels as a measure of blood model
  uncertainty.

---

## 3 — `fda_shear_limit_screening` — FDA Hemolysis Threshold Screening

### Setup

FDA critical path: scalar shear stress screening $\tau_{\mathrm{eq}} < 150$ Pa
throughout device at $Q_{\mathrm{design}}$. Equivalent shear from von Mises:

$$\tau_{\mathrm{eq}} = \sqrt{\tfrac12 (\tau_{ij}\tau_{ij})} \tag{7}$$

with per-edge estimate $\tau = 2\mu_{\mathrm{eff}}(Q/A) / (D_h/2)$ (Poiseuille
wall shear proxy).

### Tested Invariants

- **Threshold flag count**: per threshold scan $T\in\{50, 100, 150\}$ Pa,
  violation fraction reported. Zero violations at $150$ Pa for healthy
  bifurcation chip.

- **Stenosed flag**: recirculation zone downstream of stenosis exhibits
  low shear but high exposure time; both $\tau$ and $H$ flags reported
  separately.

- **CSV/JSON report**: violations machine-readable via
  `cfd-1d::screening::ShearReport { violations, max_shear_per_branch,
  threshold }`.

### Implementation

```rust
use cfd_1d::screening::shear_violations;

let report = shear_violations(&state, &geometries)?;
for branch in report.violations {
    eprintln!("branch {}: peak shear {} Pa (limit {})", branch.id, branch.peak, report.threshold);
}
```

Runs in $O(|E|)$ time.

---

## 4 — `venturi_parallel_analysis` — Parallel Venturi Pressure Drop Mapping

### Setup

Parallel network of venturi constrictions (converge-diverge axisymmetric) each
with its own $\sigma$ (cavitation number). Inlet $Q$ prescribed; pressure drop
$\Delta p_j$ per venturi from semi-empirical Young-Tsai:

$$\Delta p_j = K_v^{(j)}\frac{\mu Q_j}{A_0^2} + K_t^{(j)}\frac{\rho}{2}\left(\frac{Q_j}{A_0}\right)^2 \tag{8}$$

with coefficients $K_v, K_t$ from `cfd-schematics::venturi::VenturiSpec`.

### Tested Invariants

- **Venturi pressure recovery**: pressure recovers to $\sim 80\%$ of inlet
  $p$ downstream; validated against single-venturi 2-D reference.

- **Cavitation hot-spot**: minimum $p$ location at throat; $\sigma =
  (p_{\min}-p_v)/(0.5\rho U^2)$ reported per venturi for SDT device regime
  classification (see next example).

- **Mass conservation**: parent-to-children flow split sums to inlet.

---

## 5 — `tpms_blood_1d` — Gyroid TPMS Network with Carreau-Yasuda Blood

### Setup

TPMS scaffold approximated as a 3-D lattice of curved edges (Gyroid implicit
SDF offset $t$, porosity $\phi$). Each edge type has its own tortuosity
$\tau$ and specific surface area $S_v$ from the SDF.

Flow equations same as Section 1 but with scaffold-corrected hydraulic diameter

$$D_h = \frac{4\phi}{S_v} \tag{9}$$

and tissue-perfusion metric $J_v = L_p S_v (p_c - p_i - \sigma(\Pi_c - \Pi_i))$
for transcapillary flux (Starling equation) reported as a posterior.

### Tested Invariants

- **Porosity scaling**: $\Delta p \propto (1-\phi)/\phi^3$ per Ergun equation
  at high $\phi$.

- **Carreau vs Newtonian**: Carreau blood raises $\Delta p$ by factor
  $f_{\mathrm{CY}} \approx 1.1$–$1.4$ depending on $D_h$.

- **Perfusion uniformity**: coefficient of variation of edge velocities
  $\mathrm{CV} = \sigma_Q / \bar{Q} < 0.3$ for well-connected $\phi > 0.5$.

---

## 6 — `cavitation_venturi_analysis` — SDT Device Cavitation Screening

### Setup

Sonodynamic therapy (SDT) device: array of parallel venturi channels each at
potentially different throat diameters. Cavitation number per channel

$$\sigma_j = \frac{p_{\mathrm{inlet}}-p_{v,\mathrm{eff}}}{0.5\rho U_j^2}, \quad p_{\mathrm{eff}} = p_{\mathrm{throat},j} \tag{10}$$

Regime classification per channel:

| $\sigma$ range | Regime |
|---|---|
| $\sigma > 2$ | No cavitation |
| $1 < \sigma \le 2$ | Incipient (first bubbles) |
| $0.5 < \sigma \le 1$ | Developed cloud cavitation |
| $\sigma \le 0.5$ | Supercavitation (throat fully vapour) |

### Tested Invariants

- **$\sigma$ vs geometry**: for fixed $Q$, $\sigma$ scales as $\sigma \propto
  (R_{\mathrm{throat}}^4)/(R_{\mathrm{inlet}}^4) \times ... $ i.e. small throat
  drops pressure $p_{\mathrm{throat}}$ driving down $\sigma$.

- **Regime classification report**: JSON output with per-channel
  $\sigma_j$ and classified regime name.

- **Vapour volume estimate**: $V_{\mathrm{vapour}} \approx A_{\mathrm{cav},j}
  \times L_{\mathrm{cav},j}$ with cloud length $L_{\mathrm{cav}}$ from
  cavitation closure model (see [turbulence chapter](turbulence_multiphase.md)).

---

## 7 — `medical_millifluidic_screening` — Full Device Screening

### Setup

Combines all screening sub-modules — hemolysis (Section 3), WSS, cavitation
(Section 6) — into one screening pass over a medical device chip (serpentine
+ bifurcation + venturi mix). Intended as a pre-clinical report generator.

### Tested Invariants

- **Sieving orthogonality**: hemolysis violations, shear hot-spots, and
  cavitation risk involve different edges (e.g. stenosis throat = cavitation
  risk + hemolysis risk; long straight segment = low both).

- **Regression snapshot**: golden values for total screening metrics on
  reference chip stored in `tests/data/reference_chip_report.json`; deviations
  above tolerance indicate upstream rheology or network-solver regression.

### Implementation

```rust
use cfd_1d::screening::{full_device_screening, ScreeningConfig};

let cfg = ScreeningConfig { tau_limit: 150.0, sigma_limit: 2.0, wss_lower: 0.5 };
let report = full_device_screening(&network, &cfg)?;
report.write_json("reports/screening.json")?;
```

---

## 8 — `hemolysis_serpentine_analysis` — Giersiepen-Wurzinger Model on Serpentine

### Setup

Giersiepen-Wurzinger hemolysis model (1990)

$$H = C\,\tau^b\, t^a \tag{11}$$

$C=3.62\times10^{-5}$ (SI), $a=0.785$, $b=2.416$, with alternative Zhang
calibration. Evaluated Lagrangian along 1-D edge streamlines.

### Tested Invariants

- **Newtonian vs non-Newtonian hemolysis**: non-Newtonian edge shear produces
  higher $\tau$ → higher per-edge $H$.

- **Accumulation law**: along a multi-edge path
  $H_{\mathrm{path}} = \sum_j C \tau_j^b t_j^a$ with strong nonlinearity from
  $b\approx 2.4$ (high-shear edges dominate).

- **Turn enhancement**: serpentine U-turn minor-loss shear adds local $H$
  spikes; scaling $\propto N_{\mathrm{turns}} \cdot \mathrm{De}$.

### API

```rust
use cfd_1d::hemolysis::{GiersiepenWurzinger, HemolysisPath};

let model = GiersiepenWurzinger::standard(); // C, a, b
let path = HemolysisPath::along_edges(&network, inlet_edge, &edge_shear_fn);
let hi = path.integrated_index(&model);
```

---

## Running the Suite

```bash
cargo run -p cfd-1d --example blood_bifurcation
cargo run -p cfd-1d --example tpms_blood_1d
cargo run -p cfd-1d --example cavitation_venturi_analysis
cargo run -p cfd-1d --example medical_millifluidic_screening
cargo run -p cfd-1d --example hemolysis_serpentine_analysis
cargo run -p cfd-1d --example microfluidic_chip
cargo run -p cfd-1d --example fda_shear_limit_screening
cargo run -p cfd-1d --example venturi_parallel_analysis
```

Each emits VTK for paraView plus JSON summary consumed by `comprehensive_validation_suite`.

---

## Convergence and Accuracy

The 1-D system uses one linear equation per edge (Poiseuille resistance × flowrate
equals pressure drop); grid-independence is immediate ($N_{\mathrm{edge}}\ge 1$
sufficient for per-edge mean values). Entrance-length effects, separation, and
vena contracta are not resolved at 1-D — errors from those are quantified in
[geometry fidelity ladder](geometry_and_meshing.md) ($L_2$ deviation $1$-D vs
$2$-D/$3$-D per edge).

Pseudo-transient steady-state solves converge quadratically when using Newton
iteration on the pressure-area coupling for compliant vessels; for rigid walls
the linear graph Laplacian solves directly.

---

## Atlas Status

| Surface | Atlas type | Status |
|---|---|---|
| `NetworkSolver` Krylov | `leto::NdArray<F, Ix1>` RHS + `leto::CsrMatrix` graph Laplacian | **in-progress** (legacy `ndarray` being replaced) |
| Rheology dispatch | `eunomia::RealField` + monomorphized Casson/Carreau impls | **complete** |
| Womersley profile | analytical, `leto::ComplexField` Bessel | **complete** |
| FL correction | `cfd-1d::microvascular` | **complete** |
| Screening / hemolysis | pure Rust, no float backend change needed | **complete** |

---

## Further Reading

- Smith et al., J. Biomech. 35, 747 (2002) — 1-D vascular network theory.
- `cfd-1d` source: `crates/cfd-1d/src/{network,rheology,womersley,microvascular,screening,hemolysis,poiseuille}.rs`.
- [Biomedical flows](biomedical_flows.md) for Womersley, Murray, Fahraeus, hemolysis.
- [Geometry chapter](geometry_and_meshing.md) for blueprint → solver binding.
- [Validation suite](crate_validation.md) for 1-D Poiseuille oracle.
