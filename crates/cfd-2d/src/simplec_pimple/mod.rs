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
//! For diagonally dominant momentum stencils, SIMPLEC reduces the splitting
//! error relative to SIMPLE by using the consistent coefficient
//! $1/(A_P - \sum_{nb} A_{nb})$ in the pressure correction. The method often
//! tolerates less pressure under-relaxation, but $\alpha_p = 1$ is not a universal
//! guarantee.
//!
//! **Proof sketch**:
//! In SIMPLE, the velocity correction $\mathbf{u}' = -\mathbf{d}\,\nabla p'$ retains
//! the $\sum_{nb} A_{nb} \mathbf{u}'_{nb}$ neighbour contribution implicitly. SIMPLEC
//! closes that term with the consistent single-cell relation
//! $\mathbf{u}'_{nb} \approx \mathbf{u}'_P$, yielding
//! $\hat{d} = V/(A_P - \sum A_{nb})$. This larger correction coefficient compensates
//! for neighbour coupling and reduces the spectral radius of the outer iteration
//! when the underlying linearization remains stable.

mod algorithms;
pub mod config;
mod diagnostics;
mod interpolation;
pub mod solver;

pub use config::{AlgorithmType, SimplecPimpleConfig};
pub use solver::SimplecPimpleSolver;
