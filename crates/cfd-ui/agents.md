# cfd-ui Agent Instructions

## Purpose
CAD-like GUI for millifluidic CFD simulation. Built with wgpu for 3D mesh rendering
(offscreen render-to-texture) and standard Rust UI patterns for panels, toolbars, and
dialog-driven workflows.

## Architecture
Clean Architecture with four layers:
- **Domain**: Pure data models (scene graph, camera, drawing specs, document state).
  No dependencies on UI or GPU frameworks.
- **Application**: Use-case orchestrators (mesh creation commands, simulation runner,
  export commands). Depends on domain + CFDrs crates.
- **Infrastructure**: wgpu rendering, file dialogs, SVG/PDF drawing export.
  Depends on domain + external I/O crates.
- **Presentation**: UI views, panels, toolbars, dialogs, menus, actions.
  Depends on all layers.

## Key Integration Points
- `cfd-mesh`: `IndexedMesh<T>`, 38 primitive builders, CSG booleans, STL/OpenFOAM I/O,
  quality analysis, channel builders.
- `cfd-3d`: `FemConfig`, `FemSolver`, `StokesFlowProblem` for simulation.
- `cfd-schematics`: `NetworkBlueprint`, `ChannelSystem`, `SchematicRenderer`, presets.
- `cfd-core`: `FluidProperties`, `BoundaryCondition`, wgpu GPU compute patterns.

## Module Size
Hard limit: 500 lines production code. Advisory: 565 lines test files.

## Quality Gates
- Zero `TODO`, `FIXME`, `unimplemented!()`, `todo!()`.
- Zero build warnings, zero clippy warnings.
- Every non-trivial algorithm needs `# Theorem` + `Proof sketch` in rustdoc.
