# Chapter 10 — Leto: Arrays and Linear Algebra

## Overview

The bulk of CFDrs's tensor arithmetic — state vectors, velocity fields,
pressure fields, sparse system matrices — migrates from `ndarray::Array*` and
`nalgebra::DMatrix` to **Leto's** unified `NdArray` and sparse `CsrMatrix`
types. Leto is the **single source of truth** for dense and sparse CPU
storage in Atlas: where legacy CFDrs depended on `ndarray` for tensor storage
and `nalgebra` for dense/sparse linear algebra, Atlas provides both under one
`eunomia::FloatElement`-generic roof with Atlas-cache-friendly allocation
(`mnemosyne::Arena`) and SIMD dispatch (`hermes-simd`) wired in.

---

## The NdArray Type

```rust
pub struct NdArray<T: RealField, D: Dimension> {
    storage: Vec<T>,   // contiguous T, row-major (C-order)
    shape: D,
    strides: D::Stride, // computed from shape
}
```

`NdArray<T, D>` is the only dense type Atlas consumers see.  It replaces:

| Legacy | Atlas | Dimension notation |
|---|---|---|
| `ndarray::Array1<f64>` | `leto::NdArray<f64, Ix1>` | 1-D vector |
| `ndarray::Array2<f64>` | `leto::NdArray<f64, Ix2>` | 2-D matrix |
| `ndarray::Array3<f64>` | `leto::NdArray<f64, Ix3>` | 3-D field (solver state) |
| `ndarray::Array4<f64>` | `leto::NdArray<f64, Ix4>` | space-time slab |
| `nalgebra::DMatrix<f64>` | `leto::DMatrix<f64>` alias | `NdArray<F, Ix2>` with BLAS dispatch |

The `IxN` dimension types are const-generic-shaped where possible; the backing
storage may optionally be arena-allocated via `mnemosyne` for hot-path fields.

Layout semantics (matching `ndarray` C-order):

- **Row-major**: in `Ix2`, stride is `(ncols, 1)` (last index varies fastest).
- **Contiguous**: `storage.len() == product(shape)`. No ghost cells — ghost
  buffers are separate arrays.
- **Indexing**: `arr.get([i,j,k])` returns `Option<&T>` (bounds-checked);
  `arr.as_slice()` / `arr.as_slice_mut()` gives flat `&[T]` / `&mut [T]`
  for SIMD loops.

Allocation pattern in the solver:

```rust
use leto::{NdArray, Ix3};

let mut u = NdArray::<f64, Ix3>::zeros((nx, ny, nz));   // owned
u.fill(1.0);
let v = &mut u;   // &mut NdArray — no view type needed for mut refs
```

### Dimension Types and Const Generics

`Ix1..Ix4` expose `const fn` shape accessors and `strides` computed via const
eval. Higher-dim `u128`-wide shape available for batch-field slabs. For
solver-code that previously generic'd over dimension count, the trait
`Dimension: Copy + Clone + PartialEq + Debug` replaces manual scalar
repeating. Const-generic specializations exist for `Ix2::eye()` and
`Ix2::zeros_square(n)`.

---

## CowArray — Zero-Copy Read-Only Views

```rust
pub struct CowArray<'a, T: RealField, D: Dimension> {
    inner: Cow<'a, NdArray<T, D>>,
}
```

`CowArray` returns **borrowed** data when the underlying storage can be shared
and **owned** data when a write requires a copy.  CFD adjoint solves and
matrix-free operator applications pass `CowArray` everywhere — the typical
adjoint write is rare and small, so the borrow case dominates and there is
**zero copy** on the hot path.

Pattern in adjoint gradient accumulation:

```rust
fn accumulate_gradient<'a>(field: CowArray<'a, f64, Ix2>, adj: &mut NdArray<f64, Ix2>) {
    // if `field` borrowed (most cases): no allocation
    // if `field` owned (write-dirty adj.): field.into_owned() triggers clone
    for (g, &f) in adj.iter_mut().zip(field.iter()) { *g += f; }
}
```

Vs legacy `ndarray::CowArray`: leto's `CowArray` additionally tracks whether
the borrowed data lives in a `mnemosyne::Arena` so that lifetime violations
produce a clear error message (`Arena lifetime mismatch: borrowed view outlives
arena region`) rather than a dangling-pointer UB.

---

## CSR Sparse Storage

Spectral and finite-element stiffness matrices port to:

```rust
pub struct CsrMatrix<T: RealField> {
    pub rows: usize,
    pub cols: usize,
    pub indptr: Vec<usize>,   // length rows+1, indptr[0]=0, indptr[rows]=nnz
    pub indices: Vec<u32>,    // column indices, length nnz
    pub data: Vec<T>,         // values, length nnz
}
```

Invariant: `indptr` non-decreasing, `indices[i] < cols` for all $i$,
`indptr[rows] == data.len() == indices.len() = \text{nnz}$. Debug-validated at
construction; release ignores row/column duplicates but accumulation via
`add_into(row, col, value)` deduplicates on next sort.

Atlas `CsrMatrix` has built-in kernels:

- **`apply(x, y)`** — SpMV $y \leftarrow A x$ via `hermes-simd` bands per row.
- **`apply_add(x, y)`** — $y \leftarrow y + A x$ (accumulate).
- **`diag()`** — extract diagonal $D$ for Jacobi preconditioner.
- **`nnz_per_row_distribution()`** — load-balance metric for domain decomposition.

```rust
use leto::CsrMatrix;

let mut csr = CsrMatrix::<f64>::from_triplets(
    &[(0,0,2.0),(0,1,-1.0),(1,0,-1.0),(1,1,2.0)],
    rows=2, cols=2
);

let x = NdArray::<f64, Ix1>::from_vec(vec![1.0, 2.0]);
let mut y = NdArray::<f64, Ix1>::zeros(Ix1(2));
csr.apply(&x, &mut y);   // y = [0, 3]
```

Legacy replacement:

| Legacy | Atlas | Notes |
|---|---|---|
| `nalgebra_sparse::CsrMatrix` | `leto::CsrMatrix` | Atlas has `.apply()` routed through hermes; legacy needed manual iteration |
| Hand-rolled `csr_sum`, `csr_transpose` | `CsrMatrix::transpose()` (optional CSC intermediary) | Row→col transpose via counting sort $O(\text{nnz}+cols)$ |
| `nalgebra::DMatrix` dense matvec | `NdArray<F,Ix2> @ NdArray<F,Ix1>` via `DMatrix::gemv_scalar` | BLAS dispatch on `f32`/`f64` |
| `nalgebra::Matrix3::mul_vec` | `leto::Vector3::dot/matmul` not via mat-vec but per-cell cross product pattern | Avoid temporary 3x3 allocation — SIMD batched |

---

## DMatrix Alias and BLAS Dispatch

`DMatrix<F>` is a type alias over `NdArray<F, Ix2>` with an additional stipulation
that allocation goes through a `BLAS`-compatible layout (aligned to SIMD width,
contiguous columns within a cache block). Operations:

```rust
pub type DMatrix<F> = NdArray<F, Ix2>;

impl<F: RealField> DMatrix<F> {
    pub fn gemm(&self, other: &DMatrix<F>, out: &mut DMatrix<F>, alpha: F, beta: F);
    // out = alpha * self * other + beta * out
    pub fn gemv(&self, x: &NdArray<F, Ix1>, y: &mut NdArray<F, Ix1>, alpha: F, beta: F);
}
```

- `f32`, `f64` → dispatch to `blas_src` (`cblas_dgemm` / `sgemm` internally)
  when shape $> 64$.
- `f16` → upcast to `f32` in registers, compute, downcast (when Atlas `f16` feature
  on); `bf16` similarly.
- Generic `F: FloatElement` → scalar fallback (no BLAS), iterates inner loop.

Supersedes `nalgebra::DMatrix` where the DMatrix was used as a dense 2-D field
inside a solver — `let m: DMatrix<f64> = DMatrix::zeros(n, n)` etc.

---

## Migration Procedure (CFDrs)

A solver port typically follows:

### Step 1 — Identify the boundary

Every solver has a `from_state` API that takes the "from-`ndarray`" type.
Replace with `&NdArray<F, Ix3>` (or higher) as primary boundary — keep legacy
boundary behind `#[cfg(feature="legacy_ndarray")]` shim during intermediate
steps.

### Step 2 — Replace imports

```rust
// Before
use ndarray::{Array2, Array3, Array1, s, CowArray};
use nalgebra::{DMatrix, CsrMatrix, Point3};

// After
use leto::{NdArray, CowArray, Ix1, Ix2, Ix3, DMatrix, CsrMatrix, Point3};
use eunomia::FloatElement;
```

### Step 3 — Replace access pattern

- `arr[[i,j,k]]` → `arr.get([i,j,k])` (Option-based indexing) for safety, or
  `arr[i * ny * nz + j * nz + k]` index computation for hot inner loop where
  bounds already validated and SIMD loop uses flat slices.
- `arr.slice(s![.., j, ..])` (ndarray macro) → `arr.slice_along(Axis(1), j)` /
  direct flat-slice with stride math.
- `ndarray::Zip::indexed` → `arr.as_slice()` flat zip for SIMD.
- `.iter()` over `ndarray::Array3` → `NdArray::as_slice().iter()` still
  row-major.

Example translation for a stencil interior loop:

```rust
// Legacy ndarray
for (i,j,k) in itertools::iproduct!(1..nx-1, 1..ny-1, 1..nz-1) {
    laplacian[[i,j,k]] = (u[[i+1,j,k]] + u[[i-1,j,k]] + u[[i,j+1,k]] + u[[i,j-1,k]] + u[[i,j,k+1]] + u[[i,j,k-1]] - 6.0*u[[i,j,k]]) / (dx*dx);
}

// Atlas NdArray via flat slices + stride
let slice = u.as_slice();
let mut out = laplacian.as_slice_mut();
for idx in 0..(nx*ny*nz) {
    let (i,j,k) = unravel(idx, (nx,ny,nz));
    if i==0||i==nx-1||j==0||j==ny-1||k==0||k==nz-1 { continue; }
    let nbs = [ idx+x_stride, idx-x_stride, idx+y_stride, idx-y_stride, idx+z_stride, idx-z_stride];
    out[idx] = (nbs.iter().map(|&nb| slice[nb]).sum::<f64>() - 6.0*slice[idx])/(dx*dx);
}
```

### Step 4 — Run validation

Compare `NdArray`-ported vs legacy `.npy` outputs within the tolerance budget (see
[Migration Validation](migration_validation.md)):

```bash
cargo run -p cfd-validation --example richardson_convergence
pytest tests/ -k "test_ndarray_vs_leto_parity" -x
```

Convergence order $p \approx 2$ must be preserved; any deviation $> 10^{-8}$
relative flags a layout bug.

---

## Validation Examples

- [Example: simd_performance_benchmark](examples/simd_performance_benchmark.md) — bench SpMV through both spokes (legacy ndarray vs leto).
- [Example: matrix_free_demo](examples/matrix_free_demo.md) — matrix-free operator application over `NdArray` slices with tile streaming.
- [Example: spectral_3d_poisson](examples/spectral_3d_poisson.md) — sparse 3-D Poisson solve, ported from `nalgebra::CsrMatrix` to `leto::CsrMatrix`.

## Further Reading

- `leto` source: `D:\atlas\crates\leto\crates\` — `NdArray`, `CsrMatrix`, `DMatrix`, `Point3/Vector3`.
- [Leto: Geometry](migration_geometry.md) — point/vector geometry types.
- [Leto: GAT Tiling](migration_gat_tiles.md) — GAT-based lending-iterator tiling over `NdArray`.
- [Hermes: SIMD](migration_simd.md) — `hermes-simd` lane dispatch over `NdArray` slices.
- [Migration Overview](migration_overview.md) — Atlas motivation.
- [Migration Notes](appendix_migration.md) — full type alias table.
