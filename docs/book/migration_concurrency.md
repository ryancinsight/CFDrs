# Chapter 14 — Moirai: Concurrency

CFDrs migrates parallel execution from `rayon::ParIter` and bespoke
`std::thread::spawn` calls to **Moirai**'s unified executor.  Moirai
unifies what tokio does for async and what rayon does for CPU-bound
parallel iteration, so CFD workloads share one runtime.

## Executor Surface

```rust
pub struct Executor {
    inner: Arc<ExecutorInner>,
}

impl Executor {
    pub fn new(num_workers: usize) -> Self;
    pub fn num_cpus() -> usize;
    pub fn block_on<F: Future>(&self, fut: F) -> F::Output;
    pub fn spawn<F: Future + Send + 'static>(&self, fut: F) -> JoinHandle<F::Output>;
}
```

The runtime presents an `Executor` whose `.spawn` and `.block_on` accept
**any** future — async I/O and parallel numerical kernels share the same
worker pool.

## Task Graph

```rust
pub struct TaskGraph {
    nodes: Vec<TaskNode>,
    edges: Vec<(TaskId, TaskId)>,
}

impl TaskGraph {
    pub fn add_task(&mut self, ...) -> TaskId;
    pub fn add_edge(&mut self, from: TaskId, to: TaskId);
    pub fn execute(&self, exec: &Executor) -> Vec<TaskResult>;
}
```

`TaskGraph` is the migration surface for CFDrs when a workflow has
explicit dependencies (e.g. assemble → solve → post-process).

## Channels

```rust
pub fn channel<T: Send>(cap: usize) -> (Sender<T>, Receiver<T>);
```

Moirai channels are MPMC and back-pressure-aware.  CFDrs uses them for
streaming boundary data into a solver loop without `block_on`.

## Migration Procedure

| Legacy | Atlas |
|---|---|
| `rayon::par_iter().for_each(...)` | `exec.scope(|s| s.spawn_many(...))` |
| `rayon::join(a, b)` | `TaskGraph::add_edge` + `execute` |
| `std::thread::spawn` | `Executor::spawn` |
| `tokio::spawn` | `Executor::spawn` |
| `mpsc::channel` | `moirai::channel<T>(cap)` |

A typical CFDrs parallel region becomes:

```rust
exec.scope(|scope| {
    for tile in tiles {
        scope.spawn(async move { tile.step(&state) });
    }
});
```

The executor decides worker count, schedules across cores, and avoids
oversubscription — three wins the legacy code reproduced per-iteration.

## Validation Examples

- [`simd_performance_benchmark`](examples/simd_performance_benchmark.md) —
  parallel stencil benchmarks under Moirai.
- [`cavitation_damage_simulation`](examples/cavitation_damage_simulation.md)
  — independent prototype evaluations inside a `TaskGraph`.
- [`venturi_parallel_analysis`](examples/venturi_parallel_analysis.md) —
  fan-out sweep, ported from rayon.
- [`hemolysis_serpentine_analysis`](examples/hemolysis_serpentine_analysis.md)
  — large-scale parallel run, scoped to one Moirai executor.

## Further Reading

- [`moirai` source](../../../moirai/crates/)
- [Mnemosyne and Themis: Memory](migration_memory.md)
