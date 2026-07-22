# Appendix A — Atlas Crate Dependency Map

A CFDrs-centric view of the Atlas stack — which crate owns which surface, which
CFDrs crate consumes it, and the migration phase each consumer is in.

## Crate Map

| Atlas crate | Surface | CFDrs consumer | Status |
|---|---|---|---|
| `eunomia` | `RealField`, `ComplexField`, `FloatElement`, `IntElement` | every CFDrs crate | COMPLETE |
| `leto` | `NdArray<T,D>`, `CowArray`, geometry types | `cfd-math`, `cfd-1d`, `cfd-2d`, `cfd-3d`, `cfd-schematics` | in-progress |
| `leto` GAT | `LendingIterator`, `TileStreaming` | `cfd-2d::stencil`, `cfd-3d::stencil` | in-progress |
| `hermes-simd` | `SimdLane`, vectorized kernels | `cfd-2d::stencil`, `cfd-3d::fem`, `cfd-1d::kernel` | in-progress |
| `mnemosyne` | `Arena`, `ScratchArena` | every CFDrs crate | in-progress |
| `themis` | `Placement`, NUMA bindings | `cfd-3d`, `cfd-2d` | in-progress |
| `moirai` | `Executor`, `TaskGraph`, channels | `cfd-2d::parallel`, `cfd-3d::parallel`, `cfd-optim` | in-progress |
| `apollo` | `FftPlan` | `cfd-3d::spectral`, `cfd-2d::spectral` | graduating |
| `coeus` | `Tensor`, autodiff | n/a — no autodiff in CFDrs | n/a |
| `hephaestus` | `Backend`, GPU kernels | opt-in via `cfd-3d` GPU feature | partial |
| `ritk` | image I/O | n/a — no image data in CFDrs | n/a |
| `consus` | persistent storage | `cfd-io` (planned) | planned |
| `apollo` | spectral autodiff bridge | (planned) | planned |

## Per-Crate Migration Phase

| CFDrs crate | Migration phase | Atlas deps |
|---|---|---|
| `cfd-core` | bulk — finalize eunomia trait frontier | `eunomia`, `mnemosyne` |
| `cfd-math` | graduation — finalize leto migration | `eunomia`, `leto`, `hermes-simd`, `apollo` |
| `cfd-1d` | bulk — mid-migration | `eunomia`, `leto`, `mnemosyne`, `moirai` |
| `cfd-2d` | bulk — mid-migration | `eunomia`, `leto`, `mnemosyne`, `moirai`, `hermes-simd`, `apollo` |
| `cfd-3d` | bulk — mid-migration | `eunomia`, `leto`, `mnemosyne`, `themis`, `moirai`, `hermes-simd`, `apollo`, (optional) `hephaestus` |
| `cfd-validation` | graduation — eunomia harness complete | `eunomia`, `leto` (output) |
| `cfd-schematics` | bulk — atomic operations on point/vector typed mesh | `leto`, `hephaestus` (GPU opt-in) |
| `cfd-io` | legacy | (planned `consus`) |
| `cfd-optim` | legacy / audit | (planned `moirai`) |
| `cfd-python` | legacy | PyO3 boundary |

## Migration Order (recommended)

1. Complete `cfd-math` `leto::NdArray` migration first — it is the
   dependency for all other crates.
2. Migrate `cfd-1d` next — it has the smallest surface and exercises
   the most-used Atlas APIs.
3. Migrate `cfd-2d` and `cfd-3d` together — they share solver patterns.
4. Migrate `cfd-schematics` last (after solver migrations are stable) —
   it depends on the new geometry types becoming consistent.
5. After all solvers graduate, run the parity harness and prune the
   legacy crates.

## Cleanup Gate

Markers that the migration is complete and legacy crates can be removed:

- [ ] `cargo tree | grep -E '(ndarray|nalgebra|rayon|tokio|rustfft)'` is empty.
- [ ] `cfd-validation` parity passes on every canonical scenario.
- [ ] `validation_reports/` is updated with `migration_<crate>.md` for
      every migrated crate.
- [ ] No `unsafe` SIMD intrinsics remain outside `hermes-simd`.
- [ ] No `unsafe { Vec::from_raw_parts ... }` remains outside
      `mnemosyne`.

## Further Reading

- [Atlas Book](../../../README.md)
- [Migration Overview](migration_overview.md)
- [Leto Arrays](migration_arrays.md)
- [Leto Tiling](migration_gat_tiles.md)
- [Apollo FFT](migration_fft.md)
