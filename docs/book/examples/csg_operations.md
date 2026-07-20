# Example: csg_operations

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example csg_operations --features csg`  
**Source**: [`examples/csg_operations.rs`](../../examples/csg_operations.rs)

## What This Example Demonstrates

Union, intersection, and difference on three primitive pairs (Cube × Cylinder,
Sphere × Cube, Cylinder × Sphere) using the Mesh Arrangement pipeline, with
full topology validation and STL export to `outputs/csg/`.

| Operation | API | Check |
|---|---|---|
| Union | `csg_boolean(&a, &b, BooleanOp::Union)` | Euler χ = 2, watertight |
| Intersection | `csg_boolean(&a, &b, BooleanOp::Intersection)` | Volume ≤ min(Va,Vb) |
| Difference | `csg_boolean(&a, &b, BooleanOp::Difference)` | Volume = Va − overlap |

## Key Code Snippet

```rust
use cfd_mesh::application::csg::boolean::{csg_boolean, BooleanOp};
use cfd_mesh::{Cube, IndexedMesh};
use cfd_mesh::domain::geometry::primitives::{Cylinder, PrimitiveMesh};

let cube = Cube { origin: Point3r::new(-1.0,-1.0,-1.0),
                  width: 2.0, height: 2.0, depth: 2.0 }.build()?;
let cyl  = Cylinder { radius: 0.6, height: 3.0, segments: 24 }.build()?;

let result = csg_boolean(&cube, &cyl, BooleanOp::Difference)?;
// result written to outputs/csg/cube_minus_cylinder.stl
```

## Diagnostic Tracing

```bash
CSG_TRACE=1 cargo run --example csg_operations --features csg
```

Enables seam-repair trace output for investigating watertightness failures
after boolean operations.

## Book Chapter

[← Geometry Construction and CFD Coupling](../geometry_and_meshing.md)

