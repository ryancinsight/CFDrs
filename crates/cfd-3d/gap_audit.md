# CFD-3D Gap Audit

## Current State
`cargo check -p cfd-3d` fails with numerous errors related to a recent revision in `cfd-mesh`. 

## Identified Gaps
- **Missing Methods on `IndexedMesh`**: `face()` and `vertex()` are no longer available. `boundary_label()` signature has changed.
- **Missing Methods on Builders**: `BranchingMeshBuilder`, `SerpentineMeshBuilder`, `VenturiMeshBuilder` no longer have `build()` methods available directly without traits or their signatures changed.
- **Trait Implementations**: `num_traits::Float` vs `simba::scalar::ComplexField` conflicts (multiple applicable items in scope) for `abs()`, `min()`, `max()`, `powi()`, `powf()`, `cos()`, `sin()`.
- **Iterator Type Inference**: `fold()` operations in `bifurcation/solver.rs` erroring out because of `&f64` vs `f64` mismatch.
- **Associated Constraints**: `DMatrix::<T: Scalar>` is no longer valid in `spectral/poisson.rs` and other areas.
