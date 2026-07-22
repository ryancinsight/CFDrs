# Chapter 11 — Leto: Geometry — Points, Vectors, Isometries

CFDrs geometry — cell centres, face normals, mesh indices, coordinate
transforms — migrates from `nalgebra::Point3`/`Vector3`/`Isometry3` to the
matching **Leto** types.  Leto's geometry types share the same in-memory
layout and trait surface as their `nalgebra` analogues, but encode compile-time
distinctions between scalar types (`f32` vs `f64`), precision wrappers, and
backend (CPU vs GPU) via `PhantomData` so that callers cannot accidentally
cross wires.

## The Core Types

```rust
pub struct Point3<T: RealField>       { data: [T; 3], _marker: PhantomData<()> }
pub struct Vector3<T: RealField>      { data: [T; 3], _marker: PhantomData<()> }
pub struct Rotation3<T: RealField>    { data: Quaternion<T>, _marker: PhantomData<()> }
pub struct Translation3<T: RealField> { data: Vector3<T>, _marker: PhantomData<()> }
pub struct Isometry3<T: RealField>    { data: Rotation3<T>, translation: Vector3<T> }
```

Each type:

- Stores three scalar fields contiguously (`size_of::<Point3<f64>>() == 24`).
- Carries a ZST `PhantomData` to force the **monomorphization** of every
  method per `(F, backend)` pair.
- Implements `Copy` when `T: Copy` (true for `f32`/`f64`), so passing
  geometry into a kernel is **zero-clone**.

## Migration From nalgebra

| Legacy | Atlas |
|---|---|
| `nalgebra::Point3<f64>` | `leto::Point3<f64>` |
| `nalgebra::Vector3<f64>` | `leto::Vector3<f64>` |
| `nalgebra::UnitQuaternion<f64>` | `leto::Rotation3<f64>` |
| `nalgebra::Translation3<f64>` | `leto::Translation3<f64>` |
| `nalgebra::Isometry3<f64>` | `leto::Isometry3<f64>` |
| `Matrix3::face_normal(...)` | `Vector3::cross(a, b).normalize()` |

The type-matching imports are intentionally identical so that a CFDrs file
porting from one to the other often requires only a `use` rewrite.

## Face-Normal and Cell-Centroid Kernels

```rust
use leto::Vector3;

#[inline]
pub fn face_normal<F: FloatElement>(a: Vector3<F>, b: Vector3<F>, c: Vector3<F>) -> Vector3<F> {
    let ab = b - a;
    let ac = c - a;
    ab.cross(&ac).normalize()
}
```

This kernel ships monomorphized for `(f32, Cpu)`, `(f64, Cpu)`, and (when
GPU is wired) `(f64, Wgpu)` separately.  No virtual dispatch — the compiler
generates the right SIMD width per build.

## CSG and Schematics

The `cfd-schematics` crate uses Leto geometry for all CSG primitives — boxes,
cylinders, cones, splines.  Lathe-able surfaces (venturi, serpentine, bifurcation)
compose `Point3`/`Vector3` sequences into `BezierCurve`/`MeshHandle` types that
the solver consumes.

## Validation Examples

- [`csg_operations`](examples/csg_operations.md) — CSG primitive composition
  over `Leto::Point3`/`Leto::Vector3`.
- [`schematic_demo_integration`](examples/schematic_demo_integration.md) —
  venturi + bifurcation + serpentine schematics, ported to Atlas geometry.
- [`mesh_3d_integration`](examples/mesh_3d_integration.md) — 3-D tetrahedral
  mesh generated from `Leto::Isometry3` boundaries.

## Further Reading

- [`leto` geometry module](../../../leto/crates/)
- [Leto: Arrays](migration_arrays.md)
- [Atlas Dependency Map](appendix_dependencies.md)
