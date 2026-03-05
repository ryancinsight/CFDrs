# CFD-1D Gap Audit
- **Started:** 2026-03-05
- **Goal:** Identify numerical, physical, and performance gaps in the `cfd-1d` implementation.

## 1. Mathematical Purity & Physical Laws
- [x] Check for hardcoded constants instead of physical parameters.
- [x] Check for correct application of Hagen-Poiseuille equations.
- [x] Ensure non-Newtonian flow edge cases (if any) are explicitly prohibited or handled mathematically.
- [x] Invariant bounding: Maximum Re for strictly laminar assumption, maximum Mach for incompressibility.
  - **Resolution**: Replaced `unwrap_or(T::zero)` fallbacks in `NetworkBuilder` with explicit strict positive bounds checking for all `ChannelSpec` lengths, radii, and widths.

## 2. Performance and Memory
- [x] `Network` data structure memory locality.
- [x] Assembly of sparse matrix for linear solver.
  - **Resolution**: Refactored `cfd-1d/src/solver/matrix_assembly.rs` to replace O(N) `HashMap<usize, T>` lookups during network parallel edge iteration with O(1) continuous memory `Vec<Option<T>>` indexed arrays.
- [x] Preallocation and re-use of solver buffers.
- [x] Redundant validations taking O(N) when O(1) is possible.

## 4. CFD-2D Purity and Performance
- [x] Analyze `cfd-2d` simple algorithm implementations for memory efficiency.
  - **Resolution**: Hoisted `SparseMatrixBuilder<T>` allocation from inside `SimpleAlgorithm` `simple_iteration` hot loop to a cached generic pointer, resolving a structural $O(N)$ allocation penalty executed identically 50+ times per step.
- [x] Investigate float-injected heuristic numbers in continuous math systems.
  - **Resolution**: Replaced empirical `unwrap_or_else(|| T::one())` constant heuristics mathematically via algebraic type-safe constant equivalents (`T::one() + T::one()`) within the PISO corrector Rhie-Chow iterations.
