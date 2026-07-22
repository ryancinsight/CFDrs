# Chapter 6 — Geometry Construction and CFD Coupling

CFDrs separates solid-mesh operations from design-time microfluidic topology.
The two representations have different contracts and do not share a synthetic
cross-domain facade.

## Canonical Owners

| Concern | Owner | Public contract |
|---|---|---|
| Tessellated 3-D primitives and indexed meshes | `cfd-mesh` | `Cube`, `Cylinder`, `Cone`, `UvSphere`, `PrimitiveMesh`, `IndexedMesh` |
| Mesh Boolean operations | `cfd-mesh` | `csg_boolean(BooleanOp, &left, &right)` |
| Microfluidic network topology and layout | `cfd-schematics` | `NetworkBlueprint`, `NodeSpec`, `ChannelSpec`, preset factories |
| Physical boundary conditions | `cfd-core` | `BoundaryCondition<T>`, `BoundaryConditionSet<T>` |
| Flow discretization and solve | `cfd-1d`, `cfd-2d`, `cfd-3d` | solver-specific grid and mesh APIs |

## Solid Meshes and CSG

The root CSG examples construct `cfd-mesh` primitives through
`PrimitiveMesh::build`. They validate signed volume, watertightness, connected
components, Euler characteristic, and normal orientation. Boolean composition
uses one operation enum:

```rust
use cfd_mesh::application::csg::boolean::{csg_boolean, BooleanOp};

let union = csg_boolean(BooleanOp::Union, &left, &right)?;
let intersection = csg_boolean(BooleanOp::Intersection, &left, &right)?;
let difference = csg_boolean(BooleanOp::Difference, &left, &right)?;
```

These operations return `IndexedMesh` values. The examples optionally write
STL output after validating each result.

## Design-Time Network Blueprints

`cfd-schematics` owns planar channel topology, layout, and visualization. Its
preset factories—including `venturi_rect`, `serpentine_rect`,
`symmetric_bifurcation`, and `symmetric_trifurcation`—produce validated
`NetworkBlueprint` values. A blueprint carries nodes, channels, cross-section
specifications, topology metadata, and rendering hints. It is a design-time
network description, not a volume mesh or a 2-D PDE grid.

The `cfd-1d` solver consumes blueprint node and channel specifications for
lumped resistance-network models. Continuous `cfd-2d` and `cfd-3d` solvers
retain their own grid and mesh construction boundaries.

## Boundary Handoff

Physical conditions use
`cfd_core::physics::boundary::BoundaryCondition<T>`. Named regions are stored
in `BoundaryConditionSet<T>`; solver-specific setup code associates those
names with grid faces or mesh regions. Schematic branch metadata describes
topology and does not replace the physical boundary-condition contract.

## Examples Referenced by This Chapter

- [`csg_primitives_demo`](examples/csg_primitives_demo.md) validates primitive
  tessellation and topology.
- [`csg_operations`](examples/csg_operations.md) validates union,
  intersection, and difference results.
- [`csg_cfd_simulation`](examples/csg_cfd_simulation.md) demonstrates the
  existing CSG-to-solver example path.
- [`mesh_3d_integration`](examples/mesh_3d_integration.md) exercises a
  three-dimensional mesh integration path.
- [`dimension_scenarios_plots`](examples/dimension_scenarios_plots.md)
  compares the repository's fidelity-specific examples.
- The additional Part VI pages document the shipped schematic generator and
  preset examples with their exact Cargo target names.

## Further Reading

- [`cfd-schematics` source](../../crates/cfd-schematics/src/)
- [`cfd-core` boundary source](../../crates/cfd-core/src/physics/boundary/)
- [Biomedical and Specialized Flows](biomedical_flows.md)
- [Turbulence, Multiphase, and Cavitation](turbulence_multiphase.md)
