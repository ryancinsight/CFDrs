# Leto: Arrays and Linear Algebra

The bulk of CFDrs's tensor arithmetic — state vectors, velocity fields,
pressure fields, sparse system matrices — migrates from `ndarray::Array*` and
`nalgebra::DMatrix` to **Leto's** unified `NdArray` and sparse CsrMatrix types.
Leto is the **single source of truth** for dense and sparse CPU storage.

## The NdArray Type

```rust
pub struct NdArray<T: RealField, D: Dimension> {
    storage: Vec<T>,
    shape:   D,
}
```

`NdArray<T, D>` is the only dense type Atlas consumers see.  It replaces:

| Legacy | Atlas |
|---|---|
| `ndarray::Array1<f64>` | `leto::NdArray<f64, Ix1>` |
| `ndarray::Array2<f64>` | `leto::NdArray<f64, Ix2>` |
| `ndarray::Array3<f64>` | `leto::NdArray<f64, Ix3>` |
| `ndarray::Array4<f64>` | `leto::NdArray<f64, Ix4>` |
| `ralgebra::DMatrix<f64>` | `leto::DMatrix<f64>` (alias for `NdArray<f64, Ix2>`) |

The `IxN` dimension types are const-generic-shaped where possible; the
storage backing is a `Vec<T>` allocated via [`mnemosyne::Arena`] so that
frequent allocations do not fragment the heap.

## CowArray: Zero-Copy Read-Only Views

```rust
pub struct CowArray<'a, T: RealField, D: Dimension> {
    inner: Cow<'a, NdArray<T, D>>,
}
```

`CowArray` returns **borrowed** data when the underlying storage can be shared
and **owned** data when a write requires a copy.  CFD adjoint solves and
matrix-free operator applications pass `CowArray` everywhere — the typical
write is rare and small, so the borrow case dominates and there is **zero copy**
on the hot path.

## CSR Sparse Storage

Spectral and finite-element stiffness matrices port to:

```rust
pub struct CsrMatrix<T: RealField> {
    rows:   usize,
    cols:   usize,
    indptr: Vec<usize>,
    indices: Vec<u32>,
    data:   Vec<T>,
}
```

`CsrMatrix<T>` is monomorphised per `T`, so `f32` and `f64` builds produce
distinct kernels.  `SpMV` and `SpMM` kernels route through
[`hermes-simd`] where SIMD-width permits.

## Migration Procedure (CFDrs)

A solver port typically follows:

1. **Identify the boundary.** Every solver has a `from_state` API that takes
   the "from-`ndarray`" type.  Replace with `&NdArray<F, Ix3>` (or higher).
2. **Replace the imports.** `use ndarray::Array3;` becomes
   `use leto::NdArray; use eunomia::RealField;`.
3. **Update the access pattern.** `arr[[i,j,k]]` becomes `arr.get([i,j,k])?`
   (returns `Option<&F>`); bulk passes use `arr.as_slice()`.
4. **Run validation.** Compare `NdArray`-ported vs legacy `.npy` outputs
   within the tolerance budget (see [Migration Validation](migration_validation.md)).

## Validation Examples

- [`simd_performance_benchmark`](examples/simd_performance_benchmark.md) —
  bench SpMV through both spokes.
- [`matrix_free_demo`](examples/matrix_free_demo.md) — matrix-free operator
  application over `NdArray` slices.
- [`spectral_3d_poisson`](examples/spectral_3d_poisson.md) — sparse 3-D
  Poisson solve, ported from `nalgebra::CsrMatrix` to `leto::CsrMatrix`.

## Further Reading

- [`leto` source](../../../leto/crates/)
- [Leto: Geometry](migration_geometry.md)
- [Leto: GAT Tiling](migration_gat_tiles.md)
