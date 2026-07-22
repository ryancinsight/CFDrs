# Appendix B — Migration Notes: ndarray/nalgebra/burn → leto/hephaestus/coeus

This appendix is the **type-mapping quick reference** for the CFDrs
Atlas migration.  Use it when porting a specific call site — the
table below tells you which Atlas crate owns the surface and how the
import rewrite looks.

## Type Map

| Legacy | Atlas | Crate |
|---|---|---|
| `ndarray::Array1<f64>` | `leto::NdArray<f64, Ix1>` | `leto` |
| `ndarray::Array2<f64>` | `leto::NdArray<f64, Ix2>` | `leto` |
| `ndarray::Array3<f64>` | `leto::NdArray<f64, Ix3>` | `leto` |
| `nalgebra::DMatrix<f64>` | `leto::DMatrix<f64>` (= `NdArray<f64, Ix2>`) | `leto` |
| `nalgebra::CsrMatrix<f64>` | `leto::CsrMatrix<f64>` | `leto` |
| `nalgebra::Point3<f64>` | `leto::Point3<f64>` | `leto` |
| `nalgebra::Vector3<f64>` | `leto::Vector3<f64>` | `leto` |
| `nalgebra::Isometry3<f64>` | `leto::Isometry3<f64>` | `leto` |
| `num_traits::Float` (bound) | `eunomia::RealField` (bound) | `eunomia` |
| `num_complex::Complex64` | `eunomia::ComplexField` (`f64` impl) | `eunomia` |
| `rustfft::FftPlanner` | `apollo::FftPlan` | `apollo` |
| `realfft::RealFftPlanner` | `apollo::FftPlan::forward_real` | `apollo` |
| `rayon::ParIter` | `moirai::Executor::scope` | `moirai` |
| `rayon::join(a, b)` | `moirai::TaskGraph::add_edge` | `moirai` |
| `tokio::spawn` | `moirai::Executor::spawn` | `moirai` |
| `mpsc::channel` | `moirai::channel` | `moirai` |
| `Vec<T>::with_capacity(N)` | `mnemosyne::Arena::alloc_slice<T>(N)` | `mnemosyne` |
| `jemalloc` / `mimalloc` | `mnemosyne::Arena::with_capacity(bytes)` per subsystem | `mnemosyne` |
| `std::thread::spawn` | `moirai::Executor::spawn` | `moirai` |
| `packed_simd::f64x4` | `hermes_simd::Avx2Lane<f64>` | `hermes-simd` |
| `Iterator<Item = &Tile<'a>>` | `LendingIterator<Item<'a> = &Tile<'a>>` | `leto` |

## Migration Recipe (CFDrs)

A typical CFDrs solver port follows this recipe:

1. **Identify call sites.** Start at the solver entry point and walk
   down to leaf kernels.  Tag every callsite that touches one of the
   legacy crates in the type-map above.
2. **Replace imports.** `use ndarray::Array3` → `use leto::NdArray;
   use eunomia::RealField;`.
3. **Update bounds.** Anywhere a generic solver declared `T: Float`,
   replace the bound with `T: FloatElement` (or the tighter trait that
   your concrete kernel needs).
4. **Re-route storage.** Replace `Vec<T>` accumulator buffers with
   `&mut Arena<T>` sub-arenas; drop the per-element `drop` cost.
5. **Re-thread parallelism.** Replace `rayon::par_iter` with a
   `moirai::Executor::scope` block.  Re-think oversubscription; the
   executor already limits it.
6. **Run the parity harness.** Compare against the legacy baseline
   under [`comprehensive_validation_suite`](examples/comprehensive_validation_suite.md).

## Deprecated Paradigms

When porting, **delete**:

- `unsafe` SIMD intrinsics. Replace with [`hermes-simd`] lanes.
- Manual `csr` column indexing. Replace with [`leto::CsrMatrix`].
- Ad hoc `Arc<Mutex<Vec<T>>>` channels. Replace with [`moirai::channel`].
- `unsafe { Vec::from_raw_parts(...)` } for arena access. Replace with
  [`mnemosyne::Arena::alloc`].
- `Box<dyn FnMut(...)>` callbacks in the solver. Replace with
  trait-bound closures per the [`eunomia`] trait frontier.

## Validation Checklist

Before a migration can graduate to *cleanup* phase:

- [ ] All call sites have Atlas-compatible types (no `ndarray::Array*`,
      no `nalgebra::DMatrix`, no `Vec<f64>` accumulators).
- [ ] All parallelism routes through [`moirai::Executor`].
- [ ] All FFT calls go through [`apollo::FftPlan`].
- [ ] All numeric bounds are [`eunomia`] traits.
- [ ] Parity test passes within the documented tolerance (see
      [Migration Validation](migration_validation.md)).
- [ ] `cargo tree | grep -E '(ndarray|nalgebra|rayon|tokio|rustfft)'` is empty.

## Further Reading

- [Atlas Dependency Map](appendix_dependencies.md) — full Atlas crate map.
- [Migration Overview](migration_overview.md) — Atlas stack rationale.
- [BOOK_ORGANIZATION.md](BOOK_ORGANIZATION.md) — forward roadmap.
