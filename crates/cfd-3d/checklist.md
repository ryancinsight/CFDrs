# CFD-3D Update Checklist

- [ ] Audit compilation errors
- [ ] Fix `num_traits::Float` vs `simba::scalar::ComplexField` trait ambiguities.
- [ ] Fix missing `vertex()` and `face()` methods from `IndexedMesh`.
- [ ] Fix missing `build()` method on `SerpentineMeshBuilder`.
- [ ] Fix `DMatrix` generic parameter constraints.
- [ ] Fix `&f64` type inference in `validation.rs` / `solver.rs`.
- [ ] Verify `cargo check -p cfd-3d`
- [ ] Verify `cargo test -p cfd-3d`
