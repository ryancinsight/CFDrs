//! Numerical integration / quadrature methods for CFD applications.
//!
//! This module provides Gauss quadrature rules, composite integration, and
//! tensor-product 3D quadrature optimized for FEM assembly.
//!
//! ## Theorem — Gauss Quadrature Exactness
//!
//! **Theorem (Gauss–Legendre Exactness)**: An N-point Gauss–Legendre rule
//! integrates polynomials of degree ≤ 2N-1 exactly:
//!
//! ```text
//! ∫_{-1}^{1} p(x) dx = Σᵢ wᵢ p(xᵢ)   for all p ∈ P_{2N-1}
//! ```
//!
//! where xᵢ are the roots of the N-th Legendre polynomial Pₙ(x) and wᵢ are
//! the corresponding Christoffel weights wᵢ = 2/((1-xᵢ²)[Pₙ'(xᵢ)]²).
//!
//! **Proof sketch**: The N-point rule has 2N free parameters (N nodes + N weights).
//! Requiring exactness for P_{2N-1} (2N conditions) exactly determines the system.
//! Choosing nodes at Legendre polynomial roots maximises the degree of precision
//! by a classical optimality argument (Gaussian quadrature theorem, Hildebrand 1956).
//!
//! ## Theorem — Composite Rule Error Bound
//!
//! **Theorem (Composite p-point Rule)**: For f ∈ C^{2p}([a,b]), the composite
//! p-point Gauss rule on a uniform partition of step h satisfies:
//!
//! ```text
//! |∫_a^b f dx - Q_h[f]| ≤ C_{p,2p} · (b-a) · h^{2p} · max_{ξ}|f^{(2p)}(ξ)|
//! ```
//!
//! where C_{p,2p} depends only on the rule order.
//!
//! **Invariant**: For finite element assembly, quadrature order 2p must satisfy
//! 2p ≥ 2k-1 where k is the polynomial degree of the basis functions to achieve
//! optimal FEM convergence rates.
//!
//! ## Theorem — Fubini–Tonelli for Tensor Quadrature
//!
//! **Theorem**: For a separable integrand f(x,y,z) = g(x)·h(y)·k(z), the
//! tensor-product rule achieves exact integration if the 1D rule is exact for
//! each factor individually. For non-separable integrands, the error is bounded
//! by the supremum of cross-partial derivatives.
//!
//! ## Invariants
//! - Quadrature weights are strictly positive (wᵢ > 0) for all Gauss families.
//! - All weights sum to the integration domain measure: Σ wᵢ = 2 on [-1,1].
//! - Nodes are symmetric on [-1,1]: if xᵢ is a node then -xᵢ is a node.

pub mod composite;
pub mod quadrature;
pub mod quadrature_3d;
pub mod tensor;
pub mod traits;
pub mod utils;
pub mod variable;

// Re-export main types for convenience
pub use composite::CompositeQuadrature;
pub use quadrature::{GaussQuadrature, SimpsonsRule, TrapezoidalRule};
pub use quadrature_3d::TetrahedronQuadrature;
pub use tensor::TensorProductQuadrature;
pub use traits::{Quadrature, Quadrature3D};
pub use utils::IntegrationUtils;
pub use variable::VariableQuadrature;
