//! Generalized Minimal Residual (GMRES) method for non-symmetric linear systems
//!
//! GMRES is a Krylov subspace method for solving non-symmetric systems Ax = b.
//! This implementation follows Saad & Schultz (1986) and Saad (2003) §6.5.
//!
//! ## Algorithm Overview
//!
//! GMRES(m) builds an orthonormal basis for the Krylov subspace K_m(A, r0) using
//! the Arnoldi process, then minimizes the residual over this subspace.
//!
//! Key features:
//! - Modified Gram-Schmidt (MGS) orthogonalization for numerical stability
//! - Givens rotations for incremental least-squares solution
//! - Restart mechanism GMRES(m) to control memory usage
//! - Preconditioner support for improved convergence
//!
//! ## References
//!
//! - Saad, Y., & Schultz, M. H. (1986). GMRES: A generalized minimal residual
//!   algorithm for solving nonsymmetric linear systems. SIAM Journal on
//!   Scientific and Statistical Computing, 7(3), 856-869.
//! - Saad, Y. (2003). Iterative Methods for Sparse Linear Systems (2nd ed.).
//!   SIAM, Philadelphia, §6.5.

mod solver;

pub use solver::GMRES;
