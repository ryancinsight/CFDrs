# Chapter 6 — Geometry Construction and CFD Coupling

CFDrs geometry is built constructively: a **CSG primitive catalogue**
composes with a **mesh generator** and a **fluid-solid mapping**.  The
result is a single `Mesh` that any 1-D, 2-D, or 3-D solver consumes.

## CSG Primitives

`cfd-schematics::primitives` exposes the building blocks:

```rust
pub fn box_<F>(origin: Point3<F>, extents: Vector3<F>) -> Mesh;
pub fn cylinder<F>(origin: Point3<F>, axis: Vector3<F>, r: F, h: F) -> Mesh;
pub fn cone<F>(origin: Point3<F>, axis: Vector3<F>, r_base: F, r_tip: F, h: F) -> Mesh;
pub fn sphere_section<F>(origin: Point3<F>, r: F, theta_min: F, theta_max: F) -> Mesh;
```

Each primitive returns a `Mesh` whose vertices are `Leto::Point3<F>` and
whose elements are 2-D / 3-D references into the same `(i, j, k)` index
space as the solver's grid.

## CSG Operations

The three CSG ops are exposed through a single trait:

```rust
pub trait CsgOp {
    fn union(self, other: Mesh) -> Mesh;
    fn intersection(self, other: Mesh) -> Mesh;
    fn difference(self, other: Mesh) -> Mesh;
}
```

Operations are deferred; a CFG flag (`csgrs_backend`) selects whether
the CSG kernel runs on `csgrs` (entry-level), `brepkit` (curved-boundary),
or the Atlas [`hephaestus`]-backed GPU implementation.

## Mesh Generation

```rust
pub trait MeshGenerator<F: FloatElement> {
    fn generate(&self, geometry: &Mesh, density: &DensitySpec) -> Result<Mesh, CfdError>;
}
```

The default mesh generator produces **tetrahedral** P1 elements for 3-D
geometries and **triangular** P1 elements for 2-D.  Density is a
`DensitySpec` carrying per-region mesh-size hints.

## Fluid-Solid Mapping

A fluid-solid interface is a `Interface` with two sides — fluid and
solid — and an explicit boundary condition (slip / no-slip / partial-slip):

```rust
pub struct Interface<F: FloatElement> {
    fluid_id: RegionId,
    solid_id: RegionId,
    boundary_kind: BoundaryKind<F>,
}
```

The solver indexes its degree of freedom at the interface and applies
the chosen boundary condition.

## CSG-CFD Coupling

The full pipeline is:

```
  CSG primitive  -> CSG op  -> Mesh  -> Boundary wiring  -> Solver
```

CFDrs hides each step behind a trait so that the swap is local:

```rust
pub trait GeometryPipeline<F: FloatElement> {
    fn build_mesh(&self, primitive: &Mesh, density: &DensitySpec) -> Result<Mesh, CfdError>;
    fn wire_boundaries(&self, mesh: &Mesh) -> Result<Vec<Boundary>, CfdError>;
    fn submit(&self, mesh: &Mesh, boundaries: &[Boundary]) -> Result<SolverHandle, CfdError>;
}
```

## Examples Referenced by This Chapter

Part VI opens with five example chapters:

- [`csg_primitives_demo`](examples/csg_primitives_demo.md) — primitive
  catalogue walkthrough.
- [`csg_operations`](examples/csg_operations.md) — union / intersection /
  difference over primitives.
- [`csg_cfd_simulation`](examples/csg_cfd_simulation.md) — CSG-built
  channel integrated into a 2-D solver.
- [`mesh_3d_integration`](examples/mesh_3d_integration.md) — 3-D mesh
  generated from a CSG boundary.
- [`dimension_scenarios_plots`](examples/dimension_scenarios_plots.md) —
  1-D / 2-D / 3-D fidelity comparisons.

## Further Reading

- [`cfd-schematics` source](../../crates/cfd-schematics/src/)
- [Biomedical and Specialized Flows](biomedical_flows.md) for the
  geometries that ship pre-built.
- [Turbulence, Multiphase, and Cavitation](turbulence_multiphase.md) for
  the boundary conditions that the mesh serves.
