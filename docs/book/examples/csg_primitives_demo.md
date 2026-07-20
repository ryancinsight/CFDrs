# Example: csg_primitives_demo

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example csg_primitives_demo --features csg`  
**Source**: [`examples/csg_primitives_demo.rs`](../../examples/csg_primitives_demo.rs)

## What This Example Demonstrates

Tessellates all supported CSG primitives (Cube, UV sphere, Cylinder, Cone)
and validates each for volume accuracy, watertightness, Euler characteristic,
and outward-normal orientation — no boolean operations, pure primitive baseline.

| Primitive | Expected volume | Validation |
|---|---|---|
| Cube (side=2) | 8 m³ | Signed-volume check |
| UV sphere (r=1) | 4π/3 ≈ 4.189 m³ | Signed-volume check |
| Cylinder (r, h) | π r² h | Watertight + Euler χ = 2 |
| Cone (r, h) | π r² h / 3 | Normal orientation |

## Key Code Snippet

```rust
use cfd_mesh::{Cube, IndexedMesh, analyze_normals};
use cfd_mesh::application::watertight::check::check_watertight;
use cfd_mesh::domain::geometry::primitives::{Cone, Cylinder, PrimitiveMesh, UvSphere};

let mut mesh = Cube { origin: Point3r::new(-1.0,-1.0,-1.0),
                      width: 2.0, height: 2.0, depth: 2.0 }.build()?;
// Euler characteristic: V - E + F = 2 for any closed orientable surface
check_watertight(&mesh)?;
```

## Architecture

`cfd-mesh` (backed by the `gaia` Atlas crate) provides the topological substrate.
The `IndexedMesh` type owns vertices and cells; `check_watertight` and
`analyze_normals` operate through the `AdjacencyGraph` interface without
duplicating topology storage.

## Book Chapter

[← Geometry Construction and CFD Coupling](../geometry_and_meshing.md)

