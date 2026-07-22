# Appendix C — Atlas Glossary

A consolidated glossary spanning every Atlas crate that CFDrs, Helios, and
Kwavers consume.  Where a term has an Atlas-crate-specific meaning,
the relevant crate is named in parentheses.

## A

**Alloc, Arena Allocation** (`mnemosyne`).  A type-erased memory region
that owns its chunks and recycles them in O(1) on `reset()`.  The Atlas
default for transient solver storage.

**Apollo** (`apollo`).  The Atlas forward-only FFT crate.  Apollo
exposes `FftPlan::forward_real` and `forward_complex`.  Inverse FFT is
done by forward + complex-conjugate.

**Atlas Stack**.  The unified first-party mathematics, memory, SIMD,
and concurrency stack that replaces `ndarray`/`nalgebra`/`burn`,
`tokio`/`rayon`, `packed_simd`, `rustfft`, etc.  Atlas crates span
`eunomia`, `leto`, `hephaestus`, `coeus`, `apollo`, `hermes-simd`,
`mnemosyne`, `themis`, `moirai`, `ritk`, and `consus`.

**Autodiff / Autograd** (`coeus`).  Reverse-mode automatic
differentiation.  Coeus tracks an autodiff tape so that any
`Tensor<T>` participates in a backward graph.

## B

**Backend (Compute)** (`hephaestus`).  A `hephaestus::Backend` is one
of `Cpu`, `Wgpu`, or `Cuda`.  Atlas consumers carry a backend as a
type parameter so the same source compiles for any backend.

**BoundaryType** (`cfd-core`, `helios-domain`).  A typed enum of
boundary conditions (no-slip wall, lid, slip, partial-slip, etc.).

## C

**Carreau / Carreau-Yasuda / Casson / Power-law**.  Non-Newtonian
blood rheology models in `cfd-1d`.

**Coeus** (`coeus`).  Atlas tensor + autodiff crate.  Equivalent to
PyTorch/JAX/Burn from the Atlas perspective.

**Compile-time dispatch**.  The Atlas technique of using
`#[cfg(target_arch)]` + `PhantomData` so that the same source
compiles to specialized code per (`Scalar`, `Backend`) pair, **no
virtual dispatch**.

**ComputeBackend**.  A `eunomia::Backend` is `Cpu`, `Wgpu`, or `Cuda`.

**CowArray** (`leto`).  A `Cow<'a, NdArray<T, D>>` used for
read-then-write slicing without cloning the underlying storage.

**Cow (Copy-on-write)**.  Atlas pattern for views that
preferentially borrow; only copy when a write requires owning.

## D

**DICOM** (`ritk`).  Digital Imaging and Communications in Medicine.
The standard medical-imaging format.  Ritk parses DICOM into typed
`DicomElement`s.

**DIF, DIP**.  Dependency-Inversion Principle: depend on traits
(`eunomia::RealField`, `moirai::Executor`), not impls.

**DRY**.  Don't Repeat Yourself.  Atlas enforces DRY at the
trait frontier — there is one `RealField`, not five.

## E

**Eunomia** (`eunomia`).  Atlas numeric-trait crate that unifies
`num-traits` + `num-complex` + custom precision types.

**Executor** (`moirai`).  The unified async + parallel worker pool
that replaces `tokio::Runtime` and `rayon::ThreadPool`.

## F

**FftPlan** (`apollo`).  A pre-computed FFT plan for a fixed
shape.  Used by every Atlas consumer running spectral methods.

**FloatElement** (`eunomia`).  Atlas trait bound for floating-point
scalars that want Atlas-wide opt-in to SIMD / autodiff / MPI.

**Foundations** (`kwavers`/`cfd-core`/`helios-core`).  The lowest
chapter / layer: type definitions, validation contracts, dimension
checking.

## G

**GAT — Generic Associated Type** (Rust).  A type-level parameter
on a trait.  Atlas uses GATs for `LendingIterator::Item<'a>` and
`TileStreaming::LendingIter`.

**GpuArray / GpuBackend** (`hephaestus`).  The GPU-side counterpart
of `leto::NdArray` and `eunomia::ComputeBackend`.

## H

**Hephaestus** (`hephaestus`).  The Atlas GPU crate — wgpu (cross-
platform) and CUDA backends.

**Hermes** (`hermes-simd`).  The Atlas SIMD crate with portable
SSE2 / AVX2 / AVX-512 / NEON / WASM lane dispatch.

**Host–Device Sync**.  The boundary between `leto` (CPU) and
`hephaestus` (GPU).  Sync is mediated by `hephaestus::sync_*`
functions; Atlas arrays can flow across it with explicit
`to_device()` / `to_host()`.

## I

**Isometry3** (`leto`).  A 3-D rigid transformation (rotation +
translation).  Equivalent to `nalgebra::Isometry3`, but uses
Leto's quaternion representation.

**Ix1..IxN** (`leto`).  Const-generic dimension types for
multi-dimensional arrays.

## K

**KernelCache** (`hephaestus`).  Pre-compiled shader / CUDA
kernels, keyed by `(shape, dtype)`.  Avoids recompilation on every
launch.

**Krylov solver** (`cfd-math`, `coeus`).  An iterative linear solver
(BiCGSTAB, GMRES, IDR(s)) that needs only `apply(A, x)` for the
matrix-free path.

## L

**Lending Iterator**.  An Atlas iterator that yields a value
whose lifetime is **shorter** than `&self` — encoded via
`type Item<'a>`.  See Rust RFC 1598.

**Leto** (`leto`).  Atlas CPU dense/sparse storage crate —
`NdArray<T, D>`, `CowArray<'a, T, D>`, `CsrMatrix<T>`,
`DMatrix<T>`, and the geometry types `Point3<T>`, `Vector3<T>`,
`Isometry3<T>`.

**LtoArray / LtoView** (`leto`).  Aliases for `NdArray` and
`CowArray` respectively, used in adjoint passes.

## M

**MatrixFree** (`cfd-math`, `apollo`).  An operator-application
pattern where the linear operator is a trait and the solver
allocates only the action `A·x`.

**Migration Reference**.  The Atlas-stack migration chapter in each
of the three books (CFDrs Part VII, Helios Part VIII, Kwavers
Part VI).

**Mnemosyne** (`mnemosyne`).  Atlas arena allocator crate —
`Arena`, `ScratchArena`, with chunked growth.

**Moirai** (`moirai`).  Atlas async + parallel crate.  Replaces
`tokio::Runtime` + `rayon::ThreadPool`.

**Moiraiexecutor** (`moirai`).  The runtime type that spawns futures
and dispatches parallel work.

## N

**NdArray** (`leto`).  Atlas dense n-dimensional array.  Replaces
`ndarray::Array*` and `nalgebra::DMatrix/Tensor`.

**NUMA Placement** (`themis`).  Per-`Arena` binding to a physical
core or NUMA node — see `Placement::bind_to`,
`Placement::numa_aware`.

## P

**PhantomData**.  Zero-sized marker type that forces the Rust
compiler to honour capability requirements.  Atlas uses
`PhantomData<F>`, `PhantomData<B>`, `PhantomData<(F, B)>` to tie
scalar type and backend to every struct.

**PyO3 boundary** (`cfd-python`, `helios-python`, `kwavers-python`).
  Atlas Rust ↔ Python boundary.  All cross-boundary types must
implement `IntoPyObject` and `FromPyObject`.

## R

**Rayleigh–Plesset** (`cfd-3d::cavitation`).  Single-bubble
cavitation ODE.  Atlas' default cavitation model.

**RealField** (`eunomia`).  Atlas trait bound for any
floating-point-real-type.  The minimum surface that numeric
Atlas consumers must implement.

**Ritk** (`ritk`).  Atlas image-toolkit crate — DICOM, PNG, NIfTI.

## S

**SIMD Lane** (`hermes-simd`).  A `SimdLane`-implementing type
(`Avx2Lane<F>`, `Avx512Lane<F>`, …).  Used in vectorized kernels.

**Spectral Tile** (`apollo`, `leto`).  A frequency-domain
partition struct used by tiled spectral solvers.

**SRP**.  Single Responsibility Principle — each Atlas crate owns
one well-defined surface.

**SSOT**.  Single Source of Truth — every Atlas capability lives
in exactly one crate.

**Spectral FFT** (`apollo`, `cfd-3d::spectral`).  The Atlas FFT
pipeline (Apollo forward FFT + Atlas autodiff tape).

## T

**TaskGraph** (`moirai`).  A DAG of dependent compute tasks,
executed by an `Executor`.  Replaces `rayon::join`, `tokio::join`,
and ad-hoc thread pools.

**Tensor** (`coeus`).  An autodiff-aware n-dimensional array.

**Themis** (`themis`).  Atlas NUMA / physical-core placement
crate.  Bind Arena chunks to NUMA nodes so that solver reads
land cache-hot.

**TileStreaming** (`leto`).  A trait with a GAT-based iterator
that yields `&Tile<'a>` over an Atlas volume.

**Typed Boundary**. A CFD boundary condition carried by
`cfd_core::physics::boundary::BoundaryCondition<T>` variants such as
`VelocityInlet`, `PressureOutlet`, and `Wall`, rather than by runtime
string tags.

## Z

**Zero-Copy**.  Atlas pattern: every read-only consumer takes
`Cow<'_, T>` instead of `T`, so the underlying storage is shared
when no write is needed.

**Zero-Cost Abstraction**.  Atlas pattern: every abstraction
(ZSTs, GATs, `PhantomData`, const generics) compiles away to
zero cycles per call site.

**ZST**.  Zero-Sized Type.  Atlas uses `PhantomData<F>` and
similar ZSTs to encode capability at compile time without runtime
cost.
