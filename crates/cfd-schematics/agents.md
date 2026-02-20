# cfd-schematics — Agent Reference

> **Role**: Single authoritative source for microfluidic network topology, geometry generation,
> and 2D schematic visualisation. **Design-time only — no PDEs solved here.**  
> **Depends on**: `cfd-core`, `petgraph`, `plotters`

---

## Purpose

`cfd-schematics` is the **topology and geometry authority** for the CFDrs pipeline:

1. **Network blueprint** — `NodeSpec`, `ChannelSpec`, `CrossSectionSpec`, `NetworkBlueprint`  
2. **Geometry generation** — centreline placement, collision avoidance, 96-well-plate layout  
3. **2D schematic rendering** — SVG / PNG layout diagrams, `AnalysisOverlay`  
4. **Preset topologies** — bifurcation, trifurcation, serpentine, venturi factories  
5. **Graph persistence** — `PetgraphGraphSink` converts topology to `petgraph::DiGraph`  

Downstream crates consume the topology types:
- `cfd-1d` reads `NodeSpec` / `ChannelSpec` to build resistance networks  
- `cfd-2d` can use geometry for grid seeding and result overlays  

---

## Architecture: Clean Architecture Layers

```
Presentation
  interface/presets/          Preset factories (bifurcation, trifurcation, serpentine, venturi)

Application
  application/use_cases/      NetworkGenerationService

Infrastructure
  infrastructure/adapters/    PetgraphGraphSink, DesignGraph

Domain
  domain/model/               NodeSpec, ChannelSpec, CrossSectionSpec, NetworkBlueprint, ids
  domain/rules/               NetworkValidator (topology invariants)
```

Geometry and visualisation are *supporting* sub-systems, not domain concepts:
```
geometry/                     Builders, generator, strategies, collision detection, optimisation
state_management/             Parameter registry, bilateral symmetry, constraint satisfaction
visualizations/               SchematicRenderer, AnalysisOverlay, plotters backend
```

---

## Module Structure

```
src/
  lib.rs                      pub mods (config, geometry, state_management, visualizations,
                              application, domain, infrastructure, interface)
  config.rs                   ChannelConfig, GeometryConfig (global design parameters)
  config_constants.rs         Default channel width/height/spacing constants in mm
  error.rs                    SchematicsError

  domain/
    model/
      ids.rs                  NodeId, ChannelId (newtype wrappers)
      specs.rs                NodeSpec, ChannelSpec, CrossSectionSpec
      blueprint.rs            NetworkBlueprint (validated assembled topology)
      mod.rs
    rules/
      validator.rs            NetworkValidator: connectivity, cross-section completeness

  application/
    use_cases/
      generate_network.rs     NetworkGenerationService: orchestrates build + validate + render
    ports/
      graph_sink.rs           GraphSink<N, E> trait

  infrastructure/
    adapters/
      petgraph_graph_sink.rs  PetgraphGraphSink → petgraph::DiGraph<NodeSpec, ChannelSpec>

  geometry/
    mod.rs
    types.rs                  Point2D, Segment2D, ChannelSegment (chip plane, z=0)
    builders.rs               ChannelSystemBuilder (fluent)
    generator.rs              GeometryGenerator: node placement from blueprint
    strategies.rs             LayoutStrategy: grid, radial, custom
    collision_detection.rs    SegmentIntersectionTest (channel overlap)
    adaptive_collision.rs     Adaptive collision resolution (shift + re-test)
    optimization.rs           GeometryOptimizer: minimise total channel length
    metadata.rs               GeometryMetadata: total length, bounding box, scale
    state_integration.rs      Sync geometry with state_management registry

  state_management/
    mod.rs
    parameters.rs             ParameterSet: typed key-value store
    registry.rs               ParameterRegistry: named parameter collections
    constraints.rs            ParameterConstraint (bounds, relations)
    managers.rs               StateManager: atomic parameter updates
    bilateral_symmetry.rs     BilateralSymmetry: mirror-image layout helper
    symmetry_integration.rs   Applies symmetry to ChannelSystem geometry
    adaptive.rs               AdaptiveStateManager: dynamic reconfiguration
    validation.rs             StateValidator: consistency checks
    errors.rs                 StateError

  visualizations/
    mod.rs
    traits.rs                 Renderer<T> trait
    schematic.rs              SchematicRenderer: node+channel 2D diagram
    analysis_field.rs         AnalysisOverlay: pressure/velocity heat-map on schematic
    plotters_backend.rs       plotters (SVG/PNG) backend adapter
    shared_utilities.rs       Colour maps, label placement

  interface/
    presets/
      bifurcation.rs          BifurcationPreset: 1-to-2 Y-junction network
      trifurcation.rs         TrifurcationPreset: 1-to-3 network
      serpentine.rs           SerpentinePreset: N-bend serpentine
      venturi.rs              VenturiPreset: converging-diverging throat
```

---

## Key Domain Types

```rust
// Cross-section geometry consumed by cfd-1d resistance models
pub enum CrossSectionSpec {
    Circular  { diameter_m: f64 },
    Rectangular { width_m: f64, height_m: f64 },
}

// Single node in the schematic topology
pub struct NodeSpec {
    pub id:       NodeId,
    pub position: Point2D,          // layout plane (mm)
    pub kind:     NodeKind,         // Inlet | Junction | Outlet | Well
}

// Single channel edge
pub struct ChannelSpec {
    pub id:           ChannelId,
    pub from:         NodeId,
    pub to:           NodeId,
    pub length_m:     f64,
    pub cross_section: CrossSectionSpec,
    pub material:     ChannelMaterial,
}

// Fully validated, immutable network
pub struct NetworkBlueprint {
    pub nodes:    Vec<NodeSpec>,
    pub channels: Vec<ChannelSpec>,
}
impl NetworkBlueprint {
    /// Build and validate in one step — returns Err if topology is invalid.
    pub fn build(nodes: Vec<NodeSpec>, channels: Vec<ChannelSpec>) -> Result<Self, SchematicsError>
}
```

---

## Preset Topology Factories

```rust
use cfd_schematics::interface::presets::{BifurcationPreset, SerpentinePreset};

// Symmetric Y-bifurcation
let blueprint = BifurcationPreset::default()
    .parent_diameter(200e-6)
    .daughter_diameter(140e-6)      // Murray's law geometry — see cfd-1d::vascular::murrays_law
    .build()?;

// 6-bend serpentine (used in SDT device designs)
let blueprint = SerpentinePreset::new(n_bends: 6, width_m: 200e-6, height_m: 50e-6)
    .channel_length_m(5e-3)
    .build()?;
```

---

## 2D Schematic Coordinates vs. Simulation Domain

| Concept | `cfd-schematics` | `cfd-2d` |
|---------|-----------------|---------|
| Coordinate system | Layout plane x, y (mm) — chip diagram | Simulation domain u(x,y) — PDE fields |
| z component | Always 0 (flat schematic) | Not applicable (2D plane) |
| Purpose | topology map, visualisation | velocity/pressure field resolution |
| Consumer | design engineer, cfd-1d | cfd-2d solver |

The schematic coordinates are **layout** coordinates, not simulation domain coordinates.

---

## Visualisation

```rust
use cfd_schematics::visualizations::{SchematicRenderer, AnalysisOverlay};

// Render topology diagram
SchematicRenderer::new(&blueprint)
    .with_labels()
    .render_svg("output/schematic.svg")?;

// Overlay 1D analysis result
AnalysisOverlay::from_blueprint(&blueprint)
    .pressure_field(&node_pressures)
    .flow_rate_field(&channel_flow_rates)
    .render_svg("output/analysis.svg")?;
```

---

## Key Invariants (enforced by `NetworkValidator`)

1. Every `ChannelSpec.from` and `.to` references a valid `NodeId` in `nodes`
2. The topology graph is connected
3. At least one `Inlet` and one `Outlet` node exist
4. No degenerate channels (`length_m > 0`, cross-section area > 0)
5. For symmetric presets: bilateral symmetry is verified geometrically

---

## Prohibited Patterns

- No solver calls (no `cfd-1d`, `cfd-2d`, `cfd-3d` imports) — this is design-time only
- Do not store simulation results in `NodeSpec` / `ChannelSpec` — use `AnalysisOverlay`
- No mutable shared state in geometry generation — `GeometryGenerator` is pure-functional
- `NetworkBlueprint` is immutable once built; do not add mutation methods
