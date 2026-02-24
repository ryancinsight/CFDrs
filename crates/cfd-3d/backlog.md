# CFD-3D Update Backlog

## Strategy
Update `cfd-3d` to be compatible with the recent `cfd-mesh` revisions. The goal is a successful compile and passage of all tests for `cfd-3d`.

## Tasks
1. Fix trait resolution errors (`num_traits::Float` vs `ComplexField`)
2. Fix missing `face` and `vertex` methods on `IndexedMesh` and `builder.build()` methods.
3. Fix associated item constraints in `DMatrix::<T: Scalar>`.
4. Fix `&f64` to `f64` parameter mismatches in iterators.
