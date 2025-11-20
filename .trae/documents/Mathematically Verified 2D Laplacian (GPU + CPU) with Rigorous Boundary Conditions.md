## Observations and Gaps
- The 2D Laplacian WGSL implements Dirichlet/Neumann/Periodic with correct stencils; interior uses 3‑point (second‑order), periodic uses 5‑point (fourth‑order) leading to mixed accuracy (laplacian_2d.wgsl:28‑66,104‑128).
- Boundary type is single global `bc_type`, preventing mixed boundaries per axis; proposal: per‑axis `bc_x`, `bc_y` (laplacian_2d.wgsl:7‑16 and gpu_operators.rs:95‑105,147‑156).
- CPU fallback in `GpuLaplacianOperator2D::apply_cpu` only implements Neumann‑like mirroring; it ignores the operator’s `bc_type` and periodic 5‑point accuracy (gpu_operators.rs:399‑446). This breaks mathematical consistency.
- Heuristic `dx_inv2`/`dy_inv2` derivation uses grid‑count sniffing; replace with exact `1/(dx^2)` and `1/(dy^2)` derived from inputs to avoid hidden approximations (gpu_operators.rs:113‑139,302‑324).
- Momentum and Poisson shaders use simpler/uneven boundary handling; they should share a consistent Laplacian core or documented differences (momentum_2d.wgsl:52‑76, poisson_3d.wgsl:45‑90).

## Mathematical Specification
- Discrete Laplacian ∇² on Cartesian grid with spacing `dx`, `dy`:
  - Interior central 3‑point second derivative (order h²): `(u_{i-1} − 2u_i + u_{i+1})/dx²`, likewise for y.
  - 5‑point periodic second derivative (order h⁴): `(-1, 16, -30, 16, -1) / (12·dx²)` stencil; same for y.
  - Neumann at boundary (du/dn=0): second‑order one‑sided approximation: `∂²u/∂x²|_{i=0} ≈ (2u_0 − 5u_1 + 4u_2 − u_3)/dx²`, symmetric at `i=nx-1`; analog in y.
  - Dirichlet (u=0 at ghost): ghost point enforces zero; boundary second derivative uses `u_{-1}=0` or `u_{nx}=0` consistently.
- Operator properties and invariants:
  - Dirichlet: symmetric negative‑definite (SPD) matrix; unique solution.
  - Neumann: symmetric negative‑semidefinite with 1‑D nullspace span{constant}; zero mean constraint for solvability.
  - Periodic: symmetric negative‑semidefinite; Fourier eigenpairs with discrete eigenvalues consistent with chosen stencil.

## Implementation Changes
- WGSL shader `laplacian_2d.wgsl`:
  - Add per‑axis boundary controls: `bc_x`, `bc_y` with same coding (0: Dirichlet, 1: Neumann, 2: Periodic).
  - Use accuracy‑consistent stencils:
    - For periodic: 5‑point in both axes (already present) across full domain.
    - For Dirichlet/Neumann: use 5‑point in interior where `i∈[2..nx‑3]`, fallback to 3‑point near boundaries; same for y.
  - Retain FMA for stability; keep `@workgroup_size(16,16)` and contiguous memory access.
- Rust `GpuLaplacianOperator2D` (gpu_operators.rs):
  - Extend `BoundaryType` or add `BoundaryAxis` pair `(bc_x, bc_y)`; update uniform `Params` struct to include both.
  - Remove heuristic `dx_inv2`/`dy_inv2` detection; compute exact `1/(dx^2)`, `1/(dy^2)` in `f32`.
  - Ensure GPU `apply_gpu` packs buffers with both boundary codes; update bind group to supply only required bindings (1=input, 3=output).
  - Implement `apply_cpu` with identical stencil logic as WGSL, honoring `bc_x`, `bc_y` and 5/3‑point transitions; ensure sign conventions match.
- Momentum/Poisson consistency:
  - Option A: Refactor to reuse a shared Laplacian core (CPU and WGSL) with boundary handling; compose operators from Laplacian.
  - Option B: Document and confine differences; add feature‑gated alignment later.

## Testing Strategy (Mathematical Validation)
- Manufactured solutions:
  - Dirichlet: `u(x,y)=x²+y²` on `[0,1]^2` with `u=0` at boundary ⇒ interior ∇²u=4; check interior error < 1e‑3 (existing test; expand to mixed axes).
  - Neumann: constant `u=1` ⇒ ∇²u=0 everywhere; verify max |∇²u| < 1e‑6 (existing test).
  - Periodic: `u=sin(2πx)sin(2πy)` ⇒ ∇²u=−8π²u; verify max error < 2e‑3 (existing test; confirm with 5‑point discrete eigenvalue comparison).
- Operator invariants:
  - Symmetry: for random fields `a,b`, check `<a, L b> == <b, L a>` to machine epsilon on CPU path; spot‑check GPU vs CPU.
  - Negative (semi)definiteness: verify `<u, L u> ≤ 0` and detect nullspace (Neumann/Periodic: constant).
  - Mixed boundary axes: periodic in x, Dirichlet in y; verify correctness via MMS tailored to axis conditions.
- Property‑based tests:
  - Grid‑size variations (nx,ny∈{8,16,32,64,128}); stability near boundaries with 5→3‑point fallback.
  - dx,dy scaling invariance: `L(α·u)` scales linearly; `L` respects `1/dx²`/`1/dy²` scaling.
- GPU/CPU consistency:
  - Bitwise tolerance comparisons after shader run vs CPU for same inputs across boundary modes; assert max relative error thresholds per mode.

## Performance and Numerics
- Maintain coalesced reads: X‑neighbors contiguous; Y‑neighbors stride `nx`. Use local temporaries to reduce repeated loads.
- Exploit FMA in WGSL for 3‑point and ghost computations to reduce cancellation.
- Keep `@workgroup_size(16,16)`; dispatch `(ceil(nx/16), ceil(ny/16),1)`.
- Avoid conditional divergence hotspots by grouping checks; early return outside bounds.

## Deliverables
- Updated WGSL `laplacian_2d.wgsl` with `bc_x`, `bc_y` and unified stencil logic.
- Updated Rust `GpuLaplacianOperator2D` API and uniform packing honoring per‑axis boundaries; exact spacing uniforms.
- CPU fallback implementing identical mathematics as WGSL.
- Expanded test suite: MMS, invariants, property‑based, GPU/CPU consistency across modes and mixed axes.
- Inline module docs for invariants and boundary condition semantics; consistency note for momentum/poisson.
