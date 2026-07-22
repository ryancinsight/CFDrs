# Chapter 16 — Leto: GAT-Based Tiling and Lending Iterators

CFDrs iterates over volumetric data in **tiles** (sliding windows, boundary
extrusions, ghost-cell haloes) at multiple granularities.  The legacy approach
allocates a `Vec<Tile>` per call and pays for the heap traffic.  Atlas
**Leto** uses **generic associated types (GATs)** to encode tile iteration
without an internal allocator.

## The GAT Tile Iterator

```rust
pub trait TileStreaming<'a> {
    type Item;
    type LendingIter: LendingIterator<Item = &'a Self::Item> + 'a;
    fn tiles(&'a self) -> Self::LendingIter;
}
```

`TileStreaming::tiles` returns a **lending iterator** that borrows from the
source.  Each `next()` yields a `&Tile<'a>` that lives only as long as the
iterator state — the type system enforces no `'static` clones, so the
hot path is **zero-copy and zero-allocator**.

## The LendingIterator Trait

```rust
pub trait LendingIterator {
    type Item<'a> where Self: 'a;
    fn next(&mut self) -> Option<Self::Item<'_>>;
}
```

`LendingIterator::Item` is a **GAT** — its lifetime parameter lets the
iterator yield references whose lifetime is **shorter** than `&self`, which
is impossible with `Iterator` (whose `Item` is fixed at the trait level).

## Migration Procedure

| Legacy | Atlas |
|---|---|
| `Vec<Tile>` allocation per pass | `TileStreaming::tiles(&volume)` |
| `impl Iterator<Item = &Tile<'a>>` | `LendingIterator<Item<'a> = &Tile<'a>>` |
| per-tile `clone()` of the inner buffer | `&'a Tile<'a>` borrow |
| ghost-cell halo allocation | `halo_lazy(window_offset, Self::slice)` |

A typical CFDrs port:

```rust
for tile in volume.tiles() {
    // tile: &Tile<'_>, borrowed from `volume`
    apply_stencil(tile, weights)?;
}
// No tile allocated, no per-step clone.
```

## How Hermes Composes

When `TileStreaming::tiles` is iterated over an `NdArray<F, Ix3>`, the tile
inner slices are `&[F]` — flat enough for [`hermes-simd`] to apply
vectorized stencils without per-call conversion.  The Leto + Hermes
combination keeps the kernel-loop short and the SIMD throughput high.

## Validation Examples

- [`matrix_free_demo`](examples/matrix_free_demo.md) — tile-by-tile
  matrix-free operator, demonstrating GAT streaming.
- [`compact_serpentine_mixing`](examples/serpentine_mixing_comprehensive.md) —
  serpentine mixing analysis on tile-streamed volume.
- [`csg_cfd_simulation`](examples/csg_cfd_simulation.md) — CSG-derived
  boundary tiles fed into the solver.
- [`spectral_3d_poisson`](examples/spectral_3d_poisson.md) — spectral
  Poisson using per-tile FFT.

## Further Reading

- [`leto` GAT module](../../../leto/crates/)
- [Leto: Arrays](migration_arrays.md)
- [Migration Overview](migration_overview.md)
