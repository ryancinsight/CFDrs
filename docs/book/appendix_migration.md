# Appendix B — Migration Notes: ndarray/nalgebra/burn → leto/hephaestus/coeus

## Overview

This appendix is the **type-mapping quick reference** for the CFDrs Atlas
migration. Use it when porting a specific call site — the table below tells
you which Atlas crate owns the surface and how the import rewrite looks. For
the per-crate dependency map and migration order, see
[Atlas Crate Dependency Map](appendix_dependencies.md).

Atlas is a family replacement for an earlier mixture of third-party crates:
`ndarray` / `nalgebra` / `nalgebra-sparse` for storage and linear algebra,
`num-traits` / `num-complex` for trait bounds, `packed_simd` and raw
`std::arch` intrinsics for SIMD, `rayon` / `tokio` for concurrency,
`rustfft` / `realfft` for Fourier transforms, and `rand` crates for sampling.
Every one of those legacy crates is replaced by a single Atlas crate that
shares the same `eunomia::RealField` trait frontier.

---

## Type Map — Legacy → Atlas

| Legacy type | Atlas equivalent | Crate | Notes / semantic change |
|---|---|---|---|
| `ndarray::Array1<f64>` | `leto::NdArray<f64, Ix1>` | `leto` | Same 1-D semantics, `FloatElement`-generic |
| `ndarray::Array2<f64>` | `leto::NdArray<f64, Ix2>` | `leto` | 2-D field; 0-indexed row-major same as ndarray |
| `ndarray::Array3<f64>` | `leto::NdArray<f64, Ix3>` | `leto` | 3-D field — solver state containers |
| `ndarray::ArrayBase` | `leto::NdArray` (view or owned) or `CowArray<F,D>` | `leto` | `CowArray` replaces ndarray's view/owned union |
| `ndarray::CowArray<'a, f64, D>` | `leto::CowArray<'a, f64, D>` | `leto` | Read-only views with no clone |
| `ndarray::ArrayView2<f64>` | `let view: &NdArray<f64, Ix2>` | `leto` | Use reference rather than separate view type |
| `ndarray::ArrayViewMut2<f64>` | `let view: &mut NdArray<f64, Ix2>` | `leto` | Mut reference; alias-safe borrowing verified by type system |
| `nalgebra::DMatrix<f64>` | `leto::DMatrix<f64>` (= `NdArray<f64, Ix2>`) alias | `leto` | Atlas `DMatrix` is BLAS dispatching wrapper over `NdArray<F,Ix2>` |
| `nalgebra::CsrMatrix<f64>` | `leto::CsrMatrix<f64>` | `leto` | Atlas Csr with builtin hermes-simd `apply`, `diag` method |
| `nalgebra::DVector<f64>` | `leto::NdArray<f64, Ix1>` | `leto` | 1-D array — no separate vector type needed |
| `nalgebra::Point3<f64>` | `leto::Point3<f64>` | `leto` | `Copy` where `F: Copy`; dist method renamed `.distance(&other)` |
| `nalgebra::Vector3<f64>` | `leto::Vector3<f64>` | `leto` | `dot`, `cross`, `normalize` have same names |
| `nalgebra::Vector2<f64>` | `leto::Vector2<f64>` | `leto` | Same for 2-D solvers (advection schemes) |
| `nalgebra::Isometry3<f64>` | `leto::Isometry3<f64>` | `leto` | `*` operator composes isometries; `inverse()` |
| `nalgebra::Isometry2<f64>` | `leto::Isometry2<f64>` | `leto` | 2-D transform for grid orientation |
| `num_traits::Float` (bound) | `eunomia::RealField` (bound) | `eunomia` | `RealField = Num + PartialOrd + Copy + ...` with trig, powers |
| `num_traits::Zero + One` | already in `RealField` | `eunomia` | No longer listed separately — implied by `RealField` |
| `num_traits::Num` | already in `RealField` | `eunomia` | Same |
| `num_complex::Complex64` | `eunomia::ComplexField<F=RealField>` (`f64` impl) | `eunomia` | Spectral path complex wrapper |
| `num_complex::Complex<Float>` | `eunomia::Complex<F>` generic | `eunomia` | `F: RealField` underlying |
| `rustfft::FftPlanner<f64>` | `apollo::FftPlan` | `apollo` | Builder-pattern plan; `.forward(&mut [Complex<F>])` method |
| `realfft::RealFftPlanner<f64>` | `apollo::FftPlan::forward_real(len)` | `apollo` | Real-to-complex plan |
| `rayon::ParIter` / `par_iter()` | `moirai::Executor::scope(|s| s.spawn(...))` | `moirai` | Scope-based; nested task graph — no par_iter adapter sugar |
| `rayon::join(a, b)` | `moirai::TaskGraph::add_edge(a_node, b_node)` or `scope.join(a,b)` | `moirai` | Task graph edges or simple scope join |
| `rayon::scope(fn)` | `moirai::Executor::scope(fn)` | `moirai` | Same semantics — spawned tasks join at scope end |
| `tokio::spawn(async {...})` | `moirai::Executor::spawn(async { ... })` | `moirai` | Unified async executor |
| `mpsc::channel::<T>(cap)` | `moirai::channel::<T>(cap)` | `moirai` | Bounded typed channel |
| `std::thread::spawn(fn)` | `moirai::Executor::spawn_blocking(fn)` | `moirai` | Blocking task on ded. pool |
| `Vec<T>::with_capacity(N)` | `mnemosyne::Arena::<T>::with_capacity(N)` or direct `Bump::alloc_slice::<T>(N)` | `mnemosyne` | Arena scratch region — zero ed/zero dealloc |
| `Vec<T>::new()` accumulation | `arena.alloc(T)` per element | `mnemosyne` | Writes into sub-arena; lifetime tied to scope |
| `jemalloc` global allocator | `mnemosyne::Arena::with_capacity(bytes)` region per subsystem | `mnemosyne` | Region pinning per solver phase — placement hint via themis |
| `mimalloc` replace | Same: arena region per subsystem | `mnemosyne` | Same |
| `packed_simd::f64x4` | `hermes_simd::Avx2Lane<f64>` or `NeonLane<f64>` | `hermes-simd` | Per-target lane; scalar fallback automatic elsewhere |
| `std::arch::x86_64::_mm256_fmadd_pd` | `hermes_simd::LaneRef<f64>::mul_add(a,b,c)` | `hermes-simd` | Dispatches to target ISA |
| `Iterator<Item = &Tile<'a>>` | `LendingIterator<Item<'a> = &Tile<'a>>` GAT | `leto` (GAT feature) | Lifetime-generic associated type — allows borrowed sub-tiles |
| `rand::rngs::StdRng` + `SliceRandom` etc. | `tyche_core::RngCore + SobolSeq / LhsSampler` | `tyche-core` (Atlas rng crate) | Quasi-random + standard sampling — used by cfd-optim |
| `Iterator<Item = (i,j,u)>` stencil scan | `moirai::par_chunks` + `hermes-simd` lane for inner loop | `moirai` + `hermes-simd` | Parallel outer scan + SIMD inner loop rather than scalar Iterator |
| `Box<dyn FnMut(&mut f64)>` callbacks | Trait-bound closure `F: FnMut(&mut NdArray<F,D>)` per `eunomia` trait frontier | `eunomia` | Type-erased callbacks now trait-bound for monomorphization |
| `Arc&lt;Mutex&lt;Vec&lt;T&gt;&gt;&gt;` channel | `moirai::channel::Bounded<T>` | `moirai` | Typed bounded channel; no `Arc&lt;Mutex&lt;...&gt;&gt;` wrapping |

---

## Import Rewrite Cheat-Sheet

```rust
// Before
use ndarray::{Array2, Array3, Array1, CowArray};
use nalgebra::{DMatrix, CsrMatrix, Point3, Vector3, Isometry3};
use num_traits::{Float, Zero, One};
use num_complex::Complex64;
use rustfft::{FftPlanner, num_complex::Complex};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};

// After
use leto::{NdArray, CowArray, Ix1, Ix2, Ix3, DMatrix, CsrMatrix,
            Point3, Vector3, Isometry3};
use eunomia::{RealField, FloatElement, ComplexField};
use apollo::FftPlan;
use moirai::Executor;
use mnemosyne::Arena;
```

---

## Migration Recipe (CFDrs Solver Port)

A typical CFDrs solver port follows this recipe, tested on `cfd-2d::stencil`
(Section 2 of [numerics chapter](numerics_and_solvers.md)) and
`cfd-3d::fem::p1`:

### Step 1 — Identify Call Sites

Start at the solver entry point and walk down to leaf kernels. Tag every call
site that touches one of the legacy crates in the type-map above. Automated
assist:

```bash
grep -R "ndarray::\|nalgebra::\|rustfft::\|rayon::\|packed_simd" crates/cfd-*/src/
grep -R "std::arch::\|_mm256_\|_mm512_\|neon\|core::arch" crates/cfd-*/src/
```

Covered call sites are tagged per subtype in the `BOOK_ORGANIZATION.md`
per-crate tracker table.

### Step 2 — Replace Imports

`use ndarray::Array3` → `use leto::{NdArray, Ix3}; use eunomia::FloatElement;`.
Leave the legacy import commented or removed — legacy path gone.

### Step 3 — Update Generic Bounds

Anywhere a generic solver declared `T: Float + Clone`, replace bound with
`T: FloatElement` (or the tighter trait that the concrete kernel needs:
`RealField + Clone + Default + Debug + Send + Sync` for threaded storages).
`FloatElement` is `eunomia` alias that bundles `RealField + Default + Debug +
'static`.

Before:

```rust
fn solve<T: Float + Clone + Send + Sync>(...)
```

After:

```rust
fn solve<F: FloatElement + Send + Sync>(...) // FloatElement includes RealField + Default + Debug
```

### Step 4 — Re-Route Storage

Replace `Vec<T>` accumulator buffers with `&mut Arena<T>` sub-arenas; drop the
per-element `drop` cost. Rationale: inside a steady-state loop, a temporary
`Vec<f64>` holding per-cell residual is heap-allocated each iteration and
freed on scope exit — `mnemosyne::Arena` allocated once outside the loop
amortizes to $O(1)$ bulk dealloc.

Before:

```rust
let mut resid = vec![0.0; n];
for i in 0..n { resid[i] = rhs[i] - A.dot(i, &x); }
```

After:

```rust
let mut arena = Arena::with_capacity(n * size_of::<F>());
let resid = arena.alloc_slice_with(n, |i| rhs[i] - A.dot(i, &x));
```

Variation: for short-lived closure buffers, explicit `ScratchArena<F>`.

For IO buffers (`VtkWriter` output), plain `Vec<F>` is acceptable — those are
not hot path.

### Step 5 — Re-Thread Parallelism

Replace `rayon::par_iter` with a `moirai::Executor::scope` block. Re-think
oversubscription; the executor already limits it.

Before (rayon):

```rust
arr.par_iter_mut().for_each(|cell| { cell.update(); });
rayon::join(|| a.compute(), || b.compute());
```

After (moirai):

```rust
executor.scope(|scope| {
    for chunk in arr.chunks_mut(chunk_size) {
        scope.spawn(|_| { for cell in chunk { cell.update(); } });
    }
});
executor.join(a_task, b_task);
```

Pitfall: `rayon::current_num_threads()` → `moirai::Executor::num_threads()`.

### Step 6 — Convert SIMD

Replace `packed_simd::f64x4` or raw `core::arch::x86_64::_mm256_*` with
`hermes-simd` lanes.

Before:

```rust
use packed_simd::f64x4;
let v = f64x4::new(a0,a1,a2,a3);
let w = v + f64x4::splat(b);
```

After:

```rust
use hermes_simd::{Avx2Lane, SimdLane};
type Lane = Avx2Lane<f64>;
let v = Lane::from_slice(&[a0,a1,a2,a3]);   // 4-wide
let w = v + Lane::splat(b);
```

Fallback: when `Avx2Lane` is not available (e.g. WASM/ARM), `SimdLane` trait
provides portable ops via scalar emulation or NEON automatically.

### Step 7 — Convert FFT

Before:

```rust
use rustfft::{FftPlanner, num_complex::Complex};
let mut planner = FftPlanner::new();
let fft = planner.plan_fft_forward(n);
fft.process(&mut buffer);
```

After:

```rust
use apollo::FftPlan;
let plan = FftPlan::forward(n);
plan.execute(&mut buffer);  // buffer: &mut [ComplexField<F>]
```

For real-FFT:

```rust
let plan = FftPlan::forward_real(n);
let (real_out, imag_out) = plan.execute_with_buffers(&input, &mut scratch);
```

### Step 8 — Convert RNG / Sampling

For `cfd-optim`:

Before:

```rust
use rand::seq::SliceRandom;
let mut rng = rand::rng();
samples.shuffle(&mut rng);
```

After:

```rust
use tyche_core::{SobolSeq, LhsSampler};
let sampler = LhsSampler::new(d=4, n=100, seed=42);
let pts = sampler.sample(); // &[f64] in [0,1]^d
```

### Step 9 — Run the Parity Harness

Compare against the legacy baseline under
[`comprehensive_validation_suite`](examples/comprehensive_validation_suite.md):

```bash
cargo run -p cfd-validation --example comprehensive_validation_suite -- --atlas-ref legacy_ref.json
cargo run -p cfd-validation --example richardson_convergence
curl --max-time 30 -s <DOC_URL> || echo "skipped URL validation (offline)"
```

Atlas vs legacy agreement tolerance is $5\times10^{-8}$ relative on dense
arrays and $10^{-6}$ on eigenvalues (condition-number amplifies) — see
[migration validation](migration_validation.md).

---

## Deprecated Paradigms — Delete When Spotted

When porting, **delete entirely** (no compat shim — compatibility soup
prohibition applies):

- `unsafe { core::arch::x86_64::_mm256_fmadd_pd(...) }` intrinsics → replace
  with `hermes-simd::Lane::mul_add()`.

- `unsafe { Vec::from_raw_parts(ptr, len, cap)` } for arena access → replace
  with `mnemosyne::Arena::alloc_slice::<T>(N)`.

- `Box<dyn FnMut(...)>` callbacks in the solver → trait-bound closures per the
  `eunomia` trait frontier (monomorphized, no vtable indirect on hot path).

- Manual `csr` column indexing (raw `ptr::offset` loops) → replace with
  `leto::CsrMatrix {-row_ptr, -col_ind, -data (Slice)}` and `.apply()`.

- `Arc&lt;Mutex&lt;Vec&lt;T&gt;&gt;&gt;` channels → `moirai::channel::Bounded<T>`.

- `BitSlice` truncation tricks in `moirai` port → now `moirai::Executor`
  with proper barrier.

- `nightly`-only `generic_const_exprs` in hand-rolled specialization →
  replaced by `const fn` in `eunomia::RealField` default methods.

---

## Validation Checklist — Before Graduation to Cleanup Phase

Before a migration can move from **bulk** to **cleanup**:

- [ ] All call sites have Atlas-compatible types (no `ndarray::Array*`,
      no `nalgebra::DMatrix`, no `Vec<f64>` transient accumulators,
      no `packed_simd`).
- [ ] All parallelism routes through `moirai::Executor` (no remaining
      `rayon::join`, `rayon::par_iter`, `tokio::spawn`, or raw
      `std::thread::spawn` for numkernels/worker tasks).
- [ ] All FFT calls go through `apollo::FftPlan` — `rustfft` / `realfft`
      no longer in that crate's `Cargo.lock` closure.
- [ ] All SIMD goes via `hermes-simd` — no bare `_mm256_` / `_mm512_` /
      `_mm256_fmadd` / ARM `vaddq` intrinsics.
- [ ] All numeric bounds are `eunomia::RealField` / `FloatElement` traits —
      no `num_traits::Float` / `num-traits::Zero+One` combo surviving.
- [ ] Parity test passes within documented tolerance (see
      [migration validation](migration_validation.md)).
- [ ] `cargo tree | grep -E '(ndarray|nalgebra|rayon|tokio|rustfft|num-traits|packed_simd)'` is empty in that crate's dependency closure.
- [ ] `cargo test -p <crate> --all-features` green including downstream
      consumer crate tests that import this crate.
- [ ] Bench (`criterion`) shows Atlas path no slower than legacy (tolerate
      <1% slowd on legacy-matched machine).
- [ ] `validation_reports/migration_<crate>.md` side-by-side written with
      before/after metrics on canonical cases.
- [ ] No `#[cfg(legacy)]` / `#[cfg(atlas)]` feature gates introduced —
      old path removed directly; branch is isolation mechanism per integrity
      policy.

---

## Further Reading

- [Atlas Dependency Map](appendix_dependencies.md) — per-crate consumption + migration order.
- [Migration Overview](migration_overview.md) — Atlas motivation, layers, and forward roadmap.
- [Leto Arrays](migration_arrays.md) — `NdArray` migration detail including view semantics.
- [Hermes SIMD](migration_simd.md) — vectorization ladder and ISA dispatch.
- [Leto GAT Tiling](migration_gat_tiles.md) — lending-iterator GAT design.
- [Migration Validation](migration_validation.md) — parity harness tolerance methodology.
- [Performance and Atlas Integration Overview](performance_and_atlas.md) — performance contract.
- [`BOOK_ORGANIZATION.md`](BOOK_ORGANIZATION.md) — forward roadmap and tracker.
- Atlas book top-level: `D:\atlas\README.md` — Atlas monolith rationale.
