# Chapter 21 — 2-D and Schematic Examples

## Overview

This chapter covers schematic-driven 2-D CFD — the two-phase pattern of
designing a device in `cfd-schematics` and simulating it in `cfd-2d` or `cfd-1d`.
The schematic layer uses `NetworkBlueprint` — a JSON-serializable graph
description — to represent microfluidic chips, bifurcations, venturi networks,
and TPMS scaffold lattices. A downstream solver (`cfd-2d::physics` or
`cfd-1d::network`) ingests the blueprint, builds a structured or unstructured
grid, wires boundary conditions, and computes perfusion.

The two-phase split decouples geometry design (iterated rapidly by the user or
automated design-space searcher) from simulation (solver telemetry and
posterior computation), aligning with CFDrs's bounded-context isolation.

---

## Two-Phase Pattern

```text
Phase 1 — Design (cfd-schematics)
   Parametric spec → CSG primitives → NetworkBlueprint → SVG schematic PNG + JSON

Phase 2 — Simulation (cfd-2d or cfd-1d)
   Blueprint JSON → mesh generation → boundary wiring → network solve
   → flows, Δp, WSS, hemolysis → AnalysisOverlay PNG + posterior JSON
```

Phase 1 produces both a human-readable SVG showing the network topology and a
machine-readable JSON consumed by Phase 2. Phase 2 reports VTK fields, CSV of
edge quantities, and per-branch screening violation flags.

### Blueprint Structure

```rust
pub struct NetworkBlueprint<F> {
    pub nodes: Vec<NodeSpec<F>>,   // (x, y, type: Inlet/Outlet/Branch/Constriction)
    pub edges: Vec<EdgeSpec<F>>,   // (from_node, to_node, EdgeKind, params)
    pub global: GlobalParams<F>,   // fluid viscosity, density, pressure BC
}
```

`EdgeKind` variants: `Straight { length, diameter }`, `Venturi { throat_ratio,
cone_angle }`, `Bifurcation { angle }`, `Uturn { radius }`, `TpmsCell {
porosity, type: Gyroid|SchwarzP|Diamond }`.

---

## 1 — `bifurcation_schematic` — Murray's Law Bifurcation + Casson 2-D Solve

### Setup

Single bifurcation generated parametrically: parent diameter $5\,\mathrm{mm}$,
children $3.6\,\mathrm{mm}$ following Murray's law $r_p^3 = r_a^3 + r_b^3$.
Branches meshed with immersed boundary (IBM) at resolution $40 \times 80$ per
branch equivalent; Casson rheology $ \tau_y = 0.004$ Pa.

### Governing Equations

2-D incompressible N-S with non-Newtonian viscosity (see
[biomedical chapter](biomedical_flows.md) Eq.(1)–(10)) in IBM cut-cell form.
Boundary conditions: parabolic inlet $U_{\mathrm{avg}} = 0.2\,\mathrm{m/s}$,
Murray-optimal $Q$-split at child exits, no-slip walls with cut-cell volume
correction (Udaykumar et al. 2001).

### Tested Invariants

- **Murray flow split**: symmetric case $Q_a = Q_b = Q_{\mathrm{parent}}/2$ to
  machine epsilon (not numerical accuracy but structural invariant of equal
  hydraulic resistances).

- **Yield-stress core**: low-shear core where $\tau < \tau_y$ has
  $\dot{\gamma}=0$ and plug flow — width predicted by
  $y_{\mathrm{plug}} = 2\tau_y / (\Delta p / L)$.

- **Pressure drop vs Murray-optimal non-stenosed**: $\Delta p$ from 2-D solver
  vs 1-D Poiseuille prediction agrees within $5\%$ when entrance length is
  small relative to branch length.

- **SVG topology match**: rendered SVG edge count equals blueprint edge count;
  outlet arrows correctly labelled with $Q_a, Q_b$ values.

### Convergence

Second-order central IBM advection; pressure Poisson via multigrid ($V$-cycle).
Grid refinement: $40\times80 \to 80\times160$ reduces $L_\infty$ velocity error
by factor $\sim 4$ ($p\approx2$) in non-plug regions.

### Example Run

```bash
cargo run -p cfd-2d --example bifurcation_schematic
# Emits: examples/output/bifurcation_schematic.svg
#        examples/output/bifurcation_schematic_vtk/
```

---

## 2 — `venturi_schematic` — Venturi Topology + ISO 5167 2-D Pressure Solve

### Setup

Converge-diverge channel built from `VenturiSpec` with inlet diameter $6\,\mathrm{mm}$,
throat ratio $\beta = d/D = 0.5$, cone angles $10^\circ$ converging,
$7^\circ$ diverging per ISO 5167 venturi tube specification. Straight-channel
length before/after venturi $2L_h$ ensures validation slice is in fully
developed region.

### Governing Equations and ISO 5167 Comparison

Incompressible N-S plus Bernoulli-corrected energy for pressure drop.

ISO 5167 discharge coefficient $C(\beta, Re_D)$:

$$C = 0.9965 - 0.00653\sqrt{10^6/\mathrm{Re}_D} \cdot \beta^{2.5} \tag{1}$$

Mass flowrate $\dot{m} = C E \varepsilon (\pi/4)d^2 \sqrt{2\Delta p \rho_1}$
with expansibility factor $\varepsilon$ and velocity-of-approach factor
$E = 1/\sqrt{1-\beta^4}$.

CFDrs validates volumetric flowrate vs ISO 5167 prediction for the same
$\Delta p$:

$$Q_{\mathrm{ISO}} = C E A_t \sqrt{2\Delta p / \rho}, \quad Q_{\mathrm{CFD}} = \int_{A_t} u_z dA \tag{2}$$

with discrepancy $<2\%$ for $\beta=0.5$, $\mathrm{Re}_D=1000$–$5000$.

### Tested Invariants

- **Throat velocity scaling**: $U_t/U_{\mathrm{in}} = 1/\beta^2 = 4$ for
  incompressible continuity; deviation $<1\%$ in fully developed venturi
  away from reattachment zone.

- **Pressure recovery**: $\Delta p_{\mathrm{recovery}} / \Delta p_{\mathrm{total}} \approx 0.8$;
  consistent with single-phase venturi benchmark.

- **Separation bubble stress**: recirculation inside diverging section
  produces secondary low-shear zone; identified via Q-criterion posterior.

### Convergence

Second-order pressure, divergence residual $10^{-8}$.

---

## 3 — `serpentine_mixing_schematic` — Serpentine Mixer Schematic + Overlay

### Setup

Serpentine chip with $N=20$ turns, channel width $D_h = 200\,\mu\mathrm{m}$,
turning radius $R_c = 500\,\mu\mathrm{m}$, $\mathrm{Re}=10$, Dean number
$\mathrm{De}=\mathrm{Re}\sqrt{D_h/R_c} \approx 6.3$ (weak Dean regime).

### Governing and Mixing Metrics

Concentration advection-diffusion

$$\frac{\partial c}{\partial t} + \mathbf{u}\cdot\nabla c = D_m \nabla^2 c \tag{3}$$

with two inlet species (saline $c=0$, dyed fluid $c=1$). Mixing efficiency

$$M = 1 - \frac{\sigma_c}{\sigma_{c,\max}}, \quad \sigma_c^2 = \frac{1}{N}\sum (c_i - \bar{c})^2 \tag{4}$$

$M=0$ unmixed, $M=1$ perfectly mixed. Péclet $Pe = U D_h / D_m \sim 10^3$–$10^5$
(high advection dominance).

### Tested Invariants

- **Turn-count vs mixing**: $M$ grows $\propto \sqrt{N_{\mathrm{turns}}}$ in
  diffusion-dominant ($Pe < 100$) and $\propto N_{\mathrm{turns}}$ in
  advection-dominant ($Pe > 1000$) with Dean vortex enhancement.

- **Outlet mass conservation**: $\int c\, dA = 0.5$ (equal inlets) within
  $10^{-5}$.

- **SVG concentration striping**: downstream striping at U-turns visualized
  with blue-red colormap in overlay.

---

## 4 — `blood_venturi` — IBM + Carreau-Yasuda Blood, 50% Stenosis

### Setup

Single straight channel $6\,\mathrm{mm}$ inlet with symmetric 50% area stenosis
(venturi placed via IBM masking). Blood Carreau-Yasuda:

$$\mu(\dot{\gamma}) = \mu_\infty + (\mu_0-\mu_\infty)[1+(\lambda\dot{\gamma})^a]^{(n-1)/a} \tag{5}$$

with $\mu_\infty=0.00345$, $\mu_0=0.056$, $\lambda=3.313$, $n=0.3568$, $a=2.0$.

### Tested Invariants

- **Jet velocity amplification**: stenosed velocity max $U_{\max} = U_{\mathrm{in}} /
  (1-S)$ where $S=0.5$ area fraction → $U_{\max}=2U_{\mathrm{in}}$ (continuity).

- **Downstream recovery length**: $L_R \approx 8 D_h$ for 50% stenosis at
  $Re_D=500$; beyond $L_R$, parabolic profile with Newtonian $\mu_\infty$
  limit for fully-developed section.

- **WSS peak**: $\tau_w$ max at throat inner wall; non-Newtonian vs Newtonian
  $\tau_w$ ratio $\approx \mu_{\mathrm{eff}}/\mu_\infty$ at wall-averaged
  shear ~ 1.2× for blood.

- **Hemolysis flag**: $H = C\tau^b t^a$ accumulated along centerline vs wall
  — wall particles show $H$ orders of magnitude higher due to near-wall high
  $\tau$.

---

## 5 — `tpms_blood_2d` — TPMS Sinusoidal Wall Constrictions + Masked SIMPLE

### Setup

Channel with sinusoidal masking pattern approximating a TPMS unit cell cross-
section via $y_{\mathrm{w}}(x) = y_0 + a\sin(kx)$. Channel height varies
periodically → local porosity $\phi(x)$.

### Tested Invariants

- **Hydraulic resistance vs porosity**: modelled with Ergun-type correlation

  $$R_h \propto \frac{(1-\phi)^2}{\phi^3 D_c^2} \tag{6}$$

  where $D_c$ is cell characteristic length.

- **Vortex shedding periodic**: for $Re > 100$, shed vortices behind each
  constriction; Strouhal $St = f D_h / U \sim 0.1$.

- **Perfusion uniformity**: edge velocity variation along periodic pattern;
  compared against 1-D TPMS lattice network result for consistency.

---

## 6 — `serpentine_venturi_1d_vs_2d` — 1-D vs 2-D Fidelity Comparison

### Setup

Serpentine chip with venturi constrictions: the same `NetworkBlueprint` solved
both by `cfd-1d` network solver and `cfd-2d` IBM solver. Quantitative fidelity
comparison:

$$\epsilon_{1D/2D} = \frac{\|p_{1D} - \bar{p}_{2D}\|_{L_2}}{\|p_{2D}\|_{L_2}} \tag{7}$$

where $\bar{p}_{2D}$ is cross-section-averaged 2-D pressure.

### Findings / Tested Invariants

- **Entrance length error**: 1-D underpredicts per-edge pressure drop by
  $10$–$30\%$ in short edges $L/D < 5$ due to missing entrance/reattachment.

- **Separation error**: 1-D cannot represent recirculation; pressure recovery
  in U-turn post-stenosis is lower in 2-D.

- **Throat accuracy**: local Poiseuille-based per-edge resistance gives
  throat $\Delta p$ within $5\%$ of 2-D wedge when edge length $L/D > 2$.

---

## 7 — `schematic_demo_integration` and `geometry_integration_demo`

### Full-pipeline Glue

```rust
use cfd_schematics::blueprint::NetworkBlueprint;
use cfd_1d::network::NetworkSolver;

// Phase 1: design
let bp = NetworkBlueprint::from_json(include_str!("my_chip.json"))?;
let svg = bp.render_svg()?;        // human-readable

// Phase 2: simulate
let network = NetworkSolver::from_blueprint(&bp)?;
let solution = network.solve_steady()?;
let report = cfd_1d::screening::full_device_screening(&solution, &cfg)?;
report.write_json("screening.json")?;
report.render_overlay(&bp, "overlay.png")?;
```

These two examples demonstrate the full generate→solve→overlay→JSON pipeline
for `cfd-1d` (geometry_integration_demo) and `cfd-schematics`→`cfd-1d` bridge
(schematic_demo_integration).

### Tested Invariants

- **Blueprint roundtrip**: JSON serialize→deserialize preserves $R_j$,
  $L_j$, connectivity — regression test in CI.

- **Overlay coverage**: every edge with $Q_j > 0$ gets a non-zero color on
  the overlay (sanity for custom-coloured chips).

---

## Running the Suite

```bash
cargo run -p cfd-2d --example bifurcation_schematic
cargo run -p cfd-2d --example venturi_schematic
cargo run -p cfd-2d --example serpentine_mixing_schematic
cargo run -p cfd-2d --example blood_venturi
cargo run -p cfd-2d --example tpms_blood_2d
cargo run -p cfd-2d --example serpentine_venturi_1d_vs_2d
cargo run -p cfd-1d --example schematic_demo_integration
cargo run -p cfd-1d --example geometry_integration_demo
```

Emits SVG schematic + optional AnalysisOverlay PNG + JSON per example into
`examples/output/` (git-ignored). The `geometry_integration_demo` JSON is
consumed by `comprehensive_validation_suite` for blueprint roundtrip
regression.

---

## API Mapping

| Concern | Path |
|---|---|
| Blueprint | `cfd-schematics::blueprint::{NetworkBlueprint, NodeSpec, EdgeSpec}` |
| Bifurcation builder | `cfd-schematics::bifurcation::BifurcationBuilder` |
| Venturi spec | `cfd-schematics::venturi::VenturiSpec` + `cfd-2d::physics::venturi` |
| Serpentine builder | `cfd-schematics::serpentine::SerpentineBuilder` |
| TPMS scaffold | `cfd-schematics::tpms::{TpmsScaffold, TpmsKind}` |
| 2-D solver | `cfd-2d::physics::{vorticity_stream, ibm, muscl}` |
| IBM masking | `cfd-2d::ibm::ImmersedBoundaryMask` |
| 1-D network | `cfd-1d::network::NetworkSolver` |
| Screening | `cfd-1d::screening::{shear_violations, full_device_screening}` |
| SVG rendering | `cfd-schematics::render::SvgRenderer` |

---

## Further Reading

- `cfd-2d` source: `crates/cfd-2d/src/{grid,physics/stencil,ibm,muscl,adaptive,io}.rs`.
- `cfd-schematics` source: `crates/cfd-schematics/src/{blueprint,bifurcation,venturi,serpentine,tpms,render}.rs`.
- [Foundations](foundations.md) for BC taxonomy.
- [Geometry chapter](geometry_and_meshing.md) — CSG primitive cat and dimension lifting.
- [Biomedical chapter](biomedical_flows.md) — blood rheology and Fahraeus effect.
- [Core benchmarks](core_flows.md) — Poiseuille and entrance length effects.
- [1-D flows](crate_1d_flows.md) for 1-D network variant of these cases.
- Udaykumar et al., JCP 174, 345 (2001) — immersed boundary implementation.
- ISO 5167-1:2022 — Venturi tube discharge coefficients.
