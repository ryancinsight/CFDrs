# Example: spectral_poisson_3d (cfd-3d)

**Crate**: `cfd-3d`
**Run**: `cargo run -p cfd-3d --example spectral_poisson_3d`
**Source**: [`crates/cfd-3d/examples/spectral_poisson_3d.rs`](../../../crates/cfd-3d/examples/spectral_poisson_3d.rs)

## What This Example Demonstrates

The 3D spectral Chebyshev solver on a Poisson problem with a manufactured
analytical solution, verifying exponential convergence of the L2 and L∞ error
as the polynomial-mode count `N` increases through `{4, 6, 8, 10, 12}`.

| Aspect | Reference |
|---|---|
| PDE | `−∇²u = f` on `[0,1]³`, Dirichlet `u = 0` on all faces |
| Source term | `f = 3·π²·sin(πx)·sin(πy)·sin(πz)` |
| Exact solution | `u = sin(πx)·sin(πy)·sin(πz)` |
| Discretization | Chebyshev-Gauss-Lobatto nodes, mapped `[−1,1] → [0,1]` |
| Reference | Canuto et al. (2006), *Spectral Methods: Fundamentals in Single Domains*, Springer, Ch. 3 |

## Key Code Snippet

```rust
use cfd_3d::spectral::poisson::PoissonBoundaryCondition;
use cfd_3d::spectral::solver::PoissonProblem;
use cfd_3d::spectral::{SpectralConfig, SpectralSolver};
use leto::Array1;

let source = |x: f64, y: f64, z: f64| -> f64 {
    3.0 * std::f64::consts::PI.powi(2) * (std::f64::consts::PI * x).sin()
        * (std::f64::consts::PI * y).sin() * (std::f64::consts::PI * z).sin()
};

for &n in &[4usize, 6, 8, 10, 12] {
    let config = SpectralConfig::<f64>::new(n, n, n);
    let solver = SpectralSolver::new(config);
    let problem = PoissonProblem {
        source_term: /* covered source */,
        bc_x: PoissonBoundaryCondition::Dirichlet(Box::new(|_, _| 0.0)),
        bc_y: PoissonBoundaryCondition::Dirichlet(Box::new(|_, _| 0.0)),
        bc_z: PoissonBoundaryCondition::Dirichlet(Box::new(|_, _| 0.0)),
    };
    let sol = solver.solve(&problem)?;
    let (l2, linf) = sol.errors_against(&|x, y, z| (PI*x).sin() * (PI*y).sin() * (PI*z).sin());
    println!("N = {:>2} | L2 = {:.3e} | Linf = {:.3e}", n, l2, linf);
}
```

## Physics Background

A spectral Chebyshev method discretizes a PDE on `N+1` Chebyshev-Gauss-Lobatto
points per axis; for smooth (analytic) solutions the L2 error decays
*exponentially* in `N`, in contrast to the algebraic `O(h^k)` decay of
finite-difference schemes. This example constructs a Poisson problem whose
exact answer is a single trigonometric mode, then sweeps `N` from 4 to 12 and
verifies that the L2 and L∞ errors both decrease monotonically. A runtime
invariant — `assert!(l2_err < prev_l2)` for every `N > 4` — encodes the
monotonic-convergence expectation that the published `cfd-3d` spectral solver
contract guarantees; a violation surfaces as a panic at the offending `N`.

## Book Chapter

[← 3-D Flows](../crate_3d_flows.md)
