# cfd-1d — Agent Reference

> **Role**: 1D lumped-network solver for microfluidic channel networks.  
> **Depends on**: `cfd-core`, `cfd-schematics`, `cfd-math`, `petgraph`

---

## Purpose

`cfd-1d` solves pressure-driven flow in microfluidic channel networks using the
**electrical circuit analogy**: each channel is a lumped hydraulic resistance,
and the network is solved as a sparse linear system (Kirchhoff's laws).

**Physical scope**: fully-developed laminar flow (Re ≪ 1), Hagen-Poiseuille or
Darcy-Weisbach resistance, Newtonian and Casson blood models.

**Cost**: microseconds — ideal for rapid design-space exploration before
committing to `cfd-2d` field simulations.

### Governing Physics

Hagen-Poiseuille resistance per edge:
```
Q = ΔP / R
  where R = 128·μ·L / (π·D⁴)          (circular cross-section)
            12·μ·L / (w·h³·φ(w/h))     (rectangular, φ = Bahrami shape factor)
```

Network system (Kirchhoff):
```
G · p = f    where G is the nodal conductance matrix (sparse, symmetric, PD)
             p = nodal pressures, f = flow-rate boundary source vector
```

---

## Module Structure

```
src/
  lib.rs                    pub mods; relationship docs vs cfd-schematics + cfd-2d

  network/
    mod.rs
    node.rs                 Node<T>: pressure node (maps from NodeSpec)
    edge.rs                 EdgeProperties<T>: resistance, geometry (maps from ChannelSpec)
    graph.rs                NetworkGraph<T>: petgraph DiGraph  wrapper
    builder.rs              NetworkBuilder: From<&NetworkBlueprint>
    component_type.rs       ComponentKind enum
    metadata.rs             NetworkMetadata: node/edge counts, total length
    wrapper.rs              Network<T>: validated graph + metadata

  solver/
    mod.rs
    problem.rs              NetworkProblem<T>: BCs + body forces
    state.rs                NetworkState<T>: p, Q, residual per iteration
    linear_system.rs        Assembles G·p = f from network
    matrix_assembly.rs      Nodal conductance matrix assembly (CSR)
    convergence.rs          KirchhoffResidual check
    geometry.rs             Cross-section area / perimeter helpers for assembly
    transient_composition.rs Transient species transport (Euler / backward Euler)
    transient_droplets.rs   Droplet tracking in time-dependent flow

  channel/
    mod.rs
    cross_section.rs        CrossSection<T>: area, perimeter, hydraulic diameter
    geometry.rs             ChannelGeometry<T>: length, elevation, curvature
    flow.rs                 ChannelFlowState<T>: Q, Re, τ_wall
    solver.rs               ChannelSolver: Q, τ_wall, dp for a single channel
    surface.rs              SurfaceProperties: roughness, wetting angle

  resistance/
    mod.rs
    traits.rs               HydraulicResistance<T>
    calculator.rs           ResistanceCalculator: dispatches to correct model
    factory.rs              ResistanceFactory: selects model from channel type
    geometry.rs             Geometry helpers (hydraulic radius, entry length)
    models/
      mod.rs
      hagen_poiseuille.rs   Circular laminar: R = 128μL/πD⁴
      rectangular.rs        Rectangular: Shah-London correlation
      darcy_weisbach.rs     Turbulent / general: R = f·L/(D·A²)·ρ/2
      serpentine.rs         Additional Dean-number correction for bends
      venturi.rs            Venturi throat: Bernoulli-based resistance
      entrance.rs           Entrance-length correction (hydrodynamic development)
      membrane.rs           Membrane resistance (filtration membranes)
      traits.rs
      tests.rs

  analysis/
    mod.rs
    flow.rs                 FlowAnalysis: Q per channel, Re, velocity profile
    pressure.rs             PressureAnalysis: nodal pressures, gradients
    resistance.rs           ResistanceAnalysis: dominant resistances, network Rₜₒₜ
    performance.rs          PerformanceAnalysis: uniformity index, efficiency
    blood_safety.rs         BloodShearLimits: FDA τ_wall < 150 Pa check
    analyzers/              Typed analyser traits (flow, pressure, resistance, performance)
    error.rs                AnalysisError

  junctions/
    mod.rs
    branching/
      mod.rs
      physics.rs            T / Y / bifurcation junction pressure balance
      solver.rs             Junction solver (iterative)
      validation.rs         Kirchhoff residual at junction nodes

  vascular/
    mod.rs
    murrays_law.rs          Murray's cube law: D₀³ = Σ Dᵢ³
    bifurcation.rs          Optimal bifurcation angle (Murray-optimised)
    womersley.rs            Womersley pulsatile flow profile

  cell_separation/
    mod.rs
    properties.rs           CellProperties (RBC, WBC, platelet sizes + densities)
    margination.rs          Margination model (near-wall cell depletion layer)
    separation_model.rs     DLD / pinched-flow fractionation approximation

  components/
    mod.rs
    channels.rs             Standard channel components
    pumps.rs                Pump boundary condition (pressure/flow source)
    valves.rs               Valve (adjustable resistance)
    mixers.rs               T/Y passive mixer
    membranes.rs            Membrane resistance component
    sensors.rs              Sensor tap (zero-flow measurement point)
    constants.rs            Standard component geometry presets
    factory.rs              ComponentFactory

  scheme_bridge/
    mod.rs
    converter.rs            SchematicsConverter: NodeSpec/ChannelSpec → Network<T>
    error.rs                BridgeError (missing cross-section, unsupported geometry)
```

---

## Key APIs

### Network Construction from Schematics

```rust
use cfd_1d::scheme_bridge::SchematicsConverter;
use cfd_schematics::domain::model::NetworkBlueprint;

let blueprint: NetworkBlueprint = /* from cfd-schematics preset */;
let network = SchematicsConverter::<f64>::convert(&blueprint)?;
```

### Solving

```rust
use cfd_1d::solver::{NetworkProblem, NetworkSolver};
use cfd_core::physics::boundary::BoundaryCondition;

let problem = NetworkProblem::new(network)
    .set_inlet_pressure(inlet_node_id, 200.0)   // Pa
    .set_outlet_pressure(outlet_node_id, 0.0);

let state = NetworkSolver::solve(&problem)?;
println!("Flow rate: {} µL/min", state.flow_rate(channel_id) * 1e9 * 60.0);
```

### Blood Safety Check

```rust
use cfd_1d::analysis::blood_safety::BloodShearLimits;

let safe = BloodShearLimits::check(&state);
// Returns Err if any channel τ_wall > 150 Pa (FDA haemolysis limit)
```

---

## Key Mathematical Theorems

| Module | Theorem |
|--------|---------|
| `resistance/models/hagen_poiseuille.rs` | Hagen-Poiseuille: Q = πD⁴ΔP / (128μL) |
| `resistance/models/rectangular.rs` | Shah-London shape factor φ(w/h) |
| `vascular/murrays_law.rs` | Murray's cube law: D₀³ = Σ Dᵢ³ |
| `vascular/womersley.rs` | Womersley pulsatile profile |
| `solver/matrix_assembly.rs` | Nodal conductance matrix is symmetric positive definite |
| `vascular/bifurcation.rs` | Murray-optimal bifurcation angle minimises viscous work |

---

## Relationship to Other Crates

| Crate | Relationship |
|-------|-------------|
| `cfd-schematics` | Authority for `NodeSpec` / `ChannelSpec`; `scheme_bridge` converts them |
| `cfd-2d` | `cfd-1d` nodal pressures seed inlet/outlet BCs for single-channel 2D simulation |
| `cfd-core` | Fluid properties (`Fluid`, blood models), physical constants, `ReynoldsNumber` |
| `cfd-math` | Sparse linear solver (GMRES/CG) for the conductance system |
| `cfd-validation` | `BloodFlow1DValidation`: Womersley and Poiseuille reference comparison |

---

## Prohibited Patterns

- No 2D/3D field computations — forward to `cfd-2d` / `cfd-3d`
- No schematic coordinate geometry — `cfd-schematics` owns x/y positions; `cfd-1d` only sees `length_m`
- No turbulence models — at Re ≪ 1, flow is always laminar
- `NetworkProblem` must be fully specified (no partial boundary conditions at solve time)
