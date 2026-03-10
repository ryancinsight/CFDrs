//! Advanced pressure-velocity coupling algorithms
//!
//! This module provides SIMPLEC and PIMPLE algorithms for improved
//! convergence in incompressible flow simulations.
//!
//! ## References
//!
//! - Van Doormaal, J. P., & Raithby, G. D. (1984). Enhancements of the SIMPLE method for predicting incompressible fluid flows.
//! - OpenFOAM PIMPLE implementation
//!
//! # Theorem (SIMPLEC Improved Convergence — Van Doormaal & Raithby 1984)
//!
//! SIMPLEC (SIMPLE-Consistent) permits $\alpha_p = 1$ (no pressure under-relaxation)
//! by replacing $1/A_P$ with $1/(A_P - \sum_{nb} A_{nb})$ in the pressure correction,
//! yielding faster convergence than standard SIMPLE.
//!
//! **Proof sketch**:
//! In SIMPLE, the velocity correction $\mathbf{u}' = -\mathbf{d}\,\nabla p'$ neglects
//! the $\sum_{nb} A_{nb} \mathbf{u}'_{nb}$ neighbour terms. SIMPLEC approximates
//! $\mathbf{u}'_{nb} \approx \mathbf{u}'_P$ (consistent approximation), yielding
//! $\hat{d} = V/(A_P - \sum A_{nb})$. This larger correction coefficient compensates
//! for the neighbour neglection, removing the need for pressure under-relaxation and
//! reducing the spectral radius of the outer iteration operator.

mod algorithms;
pub mod config;
mod diagnostics;
mod interpolation;
pub mod solver;

pub use config::{AlgorithmType, SimplecPimpleConfig};
pub use solver::SimplecPimpleSolver;
