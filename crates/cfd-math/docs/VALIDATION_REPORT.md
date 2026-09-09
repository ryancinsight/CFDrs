# cfd-math Validation Report

## Mathematical Models and Algorithms
- Discrete 2D Laplacian `∇²` on Cartesian grids with spacing `dx`, `dy`.
  - Interior Dirichlet/Neumann: central 3-point second derivative `d²u/dx² ≈ (u_{i-1} − 2u_i + u_{i+1})/dx²` and analog for `y`.
  - Periodic axes: fourth-order 5-point second derivative `(-1, 16, -30, 16, -1)/(12·dx²)` and analog for `y`.
  - Neumann boundaries: one-sided second-order second derivative at boundary `∂²u/∂x²|_{i=0} ≈ (2u_0 − 5u_1 + 4u_2 − u_3)/dx²`.
  - Dirichlet boundaries: ghost value `0` enforced at boundary.
- Iterative solvers: Conjugate Gradient (CG), GMRES, BiCGSTAB; preconditioners and multigrid components for acceleration.

## Physical and Mathematical Invariants
- Dirichlet: discrete operator is symmetric negative definite; unique solution.
- Neumann/Periodic: symmetric negative semidefinite; nullspace includes constants; solvability requires zero-mean RHS.
- Linearity: `L(αu+βv)=αLu+βLv` verified.
- Symmetry: `<u, Lv> = <v, Lu>` verified under Dirichlet interiors.

## Test Suite Overview
- Manufactured solutions:
  - Dirichlet: `u=x²+y²` interior ⇒ `∇²u=4` with max error `<1e-3`.
  - Neumann: constant `u=1` ⇒ `∇²u=0` with max abs `<1e-6`.
  - Periodic: `u=sin(2πx)sin(2πy)` ⇒ `∇²u=−8π²u` with max error `<2e-3`.
  - Mixed: periodic `x`, Dirichlet `y` product solution; interior error `<5e-3`.
- Property tests:
  - Linearity and symmetry validation for Dirichlet interiors.
- Iterative methods:
  - CG solves SPD systems; residual convergence and theoretical bounds checked.

## Floating-Point Precision & Stability
- Use of fused multiply-add (FMA) in WGSL stencils minimizes cancellation.
- Exact uniforms `dx_inv2=1/(dx²)` and `dy_inv2=1/(dy²)` avoid heuristic spacing inference.
- Convergence and residual validation in CG with tolerance-based stopping.

## Dimensional Analysis
- Discrete Laplacian has units `1/length² × field`; uniforms `dx_inv2`, `dy_inv2` validated and used consistently.

## Performance Benchmarks
- Criterion benchmarks:
  - Laplacian CPU apply across sizes `{64,128,256}`.
  - CG solve on tridiagonal SPD for sizes `{128,256,512}`.

## Assumptions and Limitations
- Neumann/Periodic modes exhibit constant-field nullspace; solvers require appropriate constraints.
- Boundary-adjacent points use 3-point stencils for Dirichlet/Neumann to preserve stability.

## Coverage and Verification
- Core algorithms (Laplacian CPU, CG) covered by unit and property tests.
- GPU paths validated functionally via manufactured solutions when feature `gpu` is enabled.

## References
- Saad, Barrett, Golub & Van Loan; see inline rustdoc in solver modules for detailed theorems and bounds.
