//! Multi-stage cascade 3D FEM solver for CIF (Cascade Inertial Focusing) networks.
//!
//! This module orchestrates independent 3D FEM Navier-Stokes solves for each
//! channel segment in a multi-channel device.  The caller supplies per-channel
//! geometry and pre-computed flow rates (from a 1D Kirchhoff network solver),
//! and this module runs a 3D solve on each channel independently.
//!
//! # Theorem — Independent Channel Decomposition
//!
//! For fully-developed flow in long rectangular ducts ($L/D_h \gg 1$),
//! entrance effects are confined to the first $\sim 0.06 \cdot Re \cdot D_h$
//! of each channel.  Outside this region, each channel segment can be solved
//! independently with prescribed inlet flow rate from the 1D network solution.
//!
//! **Proof sketch**: The parabolic nature of the downstream momentum equation
//! ensures that information propagates only downstream.  Once the 1D solver
//! provides the global flow distribution, each channel's velocity field is
//! determined by its own geometry and the assigned inlet flow rate.  The
//! inter-channel coupling is fully captured by the 1D pressure-flow balance.
//!
//! # Architecture
//!
//! `CascadeSolver3D` does not depend on `cfd-schematics` or `cfd-1d`.  It
//! accepts a `Vec<CascadeChannelSpec>` with pre-computed geometry and flow
//! rates.  This keeps the crate dependency graph acyclic and allows callers
//! in `cfd-optim` to integrate 1D → 3D coupling externally.
//!
//! # Theorem — Picard Iteration Convergence for Generalised-Newtonian Stokes
//!
//! Given the viscosity function $\mu: \mathbb{R}^+ \to [\mu_\infty, \mu_0]$ is
//! Lipschitz-continuous with constant $L_\mu$, the Picard iteration
//!
//! ```text
//! μ^{(k+1)} = μ(γ̇(u^{(k)}))
//! ```
//!
//! converges to a fixed point in $H^1$ norm provided $L_\mu < \rho / C_K$
//! where $C_K$ is the Korn inequality constant and $\rho$ is the coercivity
//! constant of the bilinear form $a(\cdot, \cdot)$.
//!
//! **Proof sketch.** The mapping $T: \mu \mapsto \mu(\dot{\gamma}(\mathbf{u}(\mu)))$
//! is a contraction on $[\mu_\infty, \mu_0]$ under the stated bound. The Stokes
//! operator with bounded viscosity is coercive on $H^1_0$, so the linear solve
//! at each step is well-posed. Convergence follows from the Banach fixed-point
//! theorem.
//!
//! **Reference:** Hirn, A. (2013). "Finite element approximation of singular
//! power-law systems." *Math. Comp.* 82:1247–1268.

mod solver;
mod types;

pub use solver::CascadeSolver3D;
pub use types::{CascadeChannelSpec, CascadeConfig3D, CascadeResult3D, ChannelResult3D};
