# Chapter 13 — Mnemosyne and Themis: Memory

CFDrs migrates from `std::alloc` (jemalloc/mimalloc on production) to
**Mnemosyne** arenas plus **Themis** NUMA placement.  Two crates compose:
Mnemosyne reuses memory in type-erased arenas; Themis binds the arenas to
the physical / logical NUMA node so that a solver thread reads its state
cache-hot.

## Mnemosyne Arenas

```rust
pub struct Arena {
    chunks: Mutex<Vec<ArenaChunk>>,
    free_list: Mutex<Vec<NonNull<u8>>>,
}

impl Arena {
    pub fn with_capacity(bytes: usize) -> Self;
    pub fn alloc<T>(&self, value: T) -> &mut T;
    pub fn try_alloc<T>(&self, value: T) -> Option<&mut T>;
    pub fn reset(&mut self);
}
```

Mnemosyne's `Arena` is what every Atlas consumer sees.  The default CFDrs
solver-side workflow:

1. **Allocate once at construction.** `state: Arena<State>` instead of
   `Vec<Field>`; the arena outlives the solver.
2. **Reset between cases.** `arena.reset()` after each test/scene, freeing
   everything in O(1) without per-element `drop`.
3. **Sub-arenas for transient kernels.** Bandwidth-bound stencil application
   uses a `ScratchArena` for intermediate buffers, then drops the entire
   arena at scope exit — no `drop` per buffer.

A solver that previously alloc'd millions of temporaries per timestep sees
**zero fragmentation** and a **predictable hot-cache footprint** after the
port.

## Themis NUMA Placement

```rust
pub struct PhysicalCore(pub u32);
pub struct NumaNode(pub u32);

pub trait Placement {
    fn bind_to(core: PhysicalCore) -> Self;
    fn numa_aware() -> Self;
    fn current() -> Self;
}
```

`themis::Placement` exposes the host topology (cores, NUMA nodes, LLC
partitions).  CFDrs binds its solver pools to NUMA-aware locations:

```rust
let placement = Placement::numa_aware();
let pool = MoiraiPool::new(placement.clone(), num_workers);
let state_arena = Arena::with_capacity(state_bytes).bind(placement);
```

The binding does not move memory at runtime — it influences which
**arena chunk** the allocator serves next, so consecutive allocations land
on the same NUMA node and the solver reads cache-hot data.

## Migration Procedure

| Legacy | Atlas |
|---|---|
| global `jemalloc`/`mimalloc` init | `Arena::with_capacity(...)` per subsystem |
| `Vec<T>::with_capacity(N)` | `Arena::alloc_slice::<T>(N)` |
| per-step `drop` cost | `arena.reset()` once per case |
| `sched_setaffinity` (Linux) | `Placement::bind_to(core)` (portable) |

## Validation Examples

- [`simd_performance_benchmark`](examples/simd_performance_benchmark.md) —
  measures cache-hot throughput with `Arena`-backed state.
- [`matrix_free_demo`](examples/matrix_free_demo.md) — `ScratchArena`
  lifecycle in a matrix-free operator application.
- [`gpu_detection`](examples/gpu_detection.md) — when Hekaistos is wired,
  numbers shown here are host-side only.

## Further Reading

- [`mnemosyne` source](../../../mnemosyne/crates/)
- [`themis` source](../../../themis/src/)
