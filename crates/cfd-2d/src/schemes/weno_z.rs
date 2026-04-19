//! Fifth-order WENO-Z reconstruction.
//!
//! # Theorem
//! WENO-Z replaces the classical Jiang-Shu nonlinear weights with a global
//! smoothness correction, reducing dissipation and recovering the optimal
//! fifth-order interface accuracy more robustly near smooth critical points.
//!
//! **Proof sketch**:
//! The nonlinear weights are constructed from the standard WENO-5 smoothness
//! indicators `β_k`, together with the global indicator `τ_5 = |β_0 - β_2|`.
//! In smooth regions all `β_k = O(Δx^2)` and `τ_5 = O(Δx^5)`, so the weights
//! converge to the optimal linear weights while preserving the non-oscillatory
//! stencil selection near discontinuities.
//!
//! ## References
//!
//! - Borges, M., Carmona, M., Costa, B., & Don, W. S. (2008).
//!   "An improved weighted essentially non-oscillatory scheme for hyperbolic
//!   conservation laws." *Journal of Computational Physics*, 227(6), 3191-3211.
//! - Recent WENO reviews continue to recommend WENO-Z as the low-dissipation
//!   baseline for shock-capturing finite-difference and finite-volume schemes.

use super::weno_helpers::{weno5_candidate_fluxes, weno5_smoothness_indicators, weno5_z_weights};
use super::{constants, Grid2D, SpatialDiscretization};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Fifth-order WENO-Z scheme.
#[derive(Debug, Clone)]
pub struct WENOZ5<T: RealField + Copy> {
    epsilon: T,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive> Default for WENOZ5<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive> WENOZ5<T> {
    /// Create a WENO-Z5 scheme.
    pub fn new() -> Self {
        Self {
            epsilon: T::from_f64(constants::WENO_EPSILON).expect("analytical constant conversion"),
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> SpatialDiscretization<T> for WENOZ5<T> {
    fn compute_derivative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        let v = [
            grid.data[(i - 2, j)],
            grid.data[(i - 1, j)],
            grid.data[(i, j)],
            grid.data[(i + 1, j)],
            grid.data[(i + 2, j)],
        ];

        let beta = weno5_smoothness_indicators(&v);
        let w = weno5_z_weights(self.epsilon, &beta);
        let flux = weno5_candidate_fluxes(&v);

        (w[0] * flux[0] + w[1] * flux[1] + w[2] * flux[2]) / grid.dx
    }

    fn order(&self) -> usize {
        5
    }

    fn is_conservative(&self) -> bool {
        true
    }

    fn cfl_limit(&self) -> f64 {
        0.5
    }
}
