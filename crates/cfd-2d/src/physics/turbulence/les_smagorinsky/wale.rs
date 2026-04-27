//! # Wall-Adapting Local Eddy-viscosity (WALE) Model
//!
//! ## Mathematical Foundation
//!
//! The WALE model (Nicoud & Ducros, 1999) improves near-wall behavior compared
//! to the Smagorinsky model by using a different invariant of the velocity
//! gradient tensor that better captures the turbulence structures near walls.
//!
//! ### WALE SGS Viscosity
//!
//! The subgrid-scale viscosity is computed as:
//!
//! ```math
//! \nu_{sgs} = (C_w \Delta)^2 \frac{(S_{ij}^d S_{ij}^d)^{3/2}}{(\bar{S}_{ij} \bar{S}_{ij})^{5/2} + (S_{ij}^d S_{ij}^d)^{5/4}}
//! ```
//!
//! where:
//! - $C_w$ is the WALE model constant (typically 0.5)
//! - $\Delta$ is the filter width
//! - $\bar{S}_{ij}$ is the resolved strain rate tensor
//! - $S_{ij}^d$ is the traceless symmetric part of the square of the velocity gradient tensor
//!
//! ### WALE Tensor $S_{ij}^d$
//!
//! The WALE tensor is computed as:
//!
//! ```math
//! S_{ij}^d = \frac{1}{2} (\partial_i u_j + \partial_j u_i)^2 - \frac{1}{3} \delta_{ij} (\partial_k u_l)^2
//! ```
//!
//! This formulation ensures that $S_{ij}^d$ vanishes at walls and provides better
//! representation of near-wall turbulence structures.
//!
//! ## References
//!
//! - Nicoud, F. & Ducros, F. (1999). Subgrid-scale stress modelling based on the square of the velocity gradient tensor. Flow, Turbulence and Combustion, 62(1), 183-200.
//! - Ducros, F., et al. (1998). Wall-adapting local eddy-viscosity model. AIAA Paper 98-0714.
//!
//! # Theorem
//! The turbulence model must satisfy the realizability conditions for the Reynolds stress tensor.
//!
//! **Proof sketch**:
//! For any turbulent flow, the Reynolds stress tensor $\tau_{ij} = -\rho \overline{u_i^\prime u_j^\prime}$
//! must be positive semi-definite. This requires that the turbulent kinetic energy $k \ge 0$
//! and the normal stresses $\overline{u_i^\prime u_i^\prime} \ge 0$. The implemented model
//! enforces these constraints either through exact transport equations or bounded eddy-viscosity
//! formulations, ensuring physical realizability and numerical stability.

use crate::fields::Field2D;
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;

/// WALE model constant (literature value: 0.5)
const C_WALE: f64 = 0.5;

/// Wall-Adapting Local Eddy-viscosity (WALE) SGS model
#[derive(Debug, Clone)]
pub struct WaleModel<T: RealField + Copy + FromPrimitive> {
    /// WALE model constant
    c_w: T,
}

impl<T: RealField + Copy + FromPrimitive> Default for WaleModel<T> {
    fn default() -> Self {
        Self {
            c_w: T::from_f64(C_WALE).expect("analytical constant conversion"),
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> WaleModel<T> {
    /// Create new WALE model
    pub fn new() -> Self {
        Self::default()
    }

    /// Create WALE model with custom constant
    pub fn with_constant(c_w: T) -> Self {
        Self { c_w }
    }

    /// Compute WALE SGS viscosity at a grid point
    ///
    /// # Arguments
    /// * `velocity` - Velocity field
    /// * `i`, `j` - Grid indices
    /// * `dx`, `dy` - Grid spacings
    /// * `delta` - Filter width (typically (dx*dy*dz)^(1/3))
    ///
    /// # Returns
    /// SGS viscosity ν_sgs
    pub fn sgs_viscosity(
        &self,
        velocity: &Field2D<Vector2<T>>,
        i: usize,
        j: usize,
        dx: T,
        dy: T,
        delta: T,
    ) -> T {
        // Compute velocity gradients
        let (du_dx, du_dy, dv_dx, dv_dy) = self.velocity_gradients(velocity, i, j, dx, dy);

        // Compute strain rate tensor components
        let s_xx = du_dx;
        let s_xy = T::from_f64(0.5).expect("analytical constant conversion") * (du_dy + dv_dx);
        let s_yy = dv_dy;

        // Strain rate magnitude squared: |S|² = 2*S_ij*S_ij
        let strain_mag_sq = T::from_f64(2.0).expect("analytical constant conversion")
            * (s_xx * s_xx + s_xy * s_xy + s_yy * s_yy);

        // Compute WALE tensor S_ij^d
        let wale_tensor_mag_sq = self.wale_tensor_magnitude_squared(du_dx, du_dy, dv_dx, dv_dy);

        // Avoid division by zero
        let epsilon = T::from_f64(1e-12).expect("analytical constant conversion");

        // WALE SGS viscosity formula
        let delta_sq = delta * delta;
        let c_w_delta_sq = self.c_w * self.c_w * delta_sq;

        if wale_tensor_mag_sq <= epsilon {
            // Near walls where WALE tensor vanishes, fall back to small constant
            T::zero()
        } else {
            // Full WALE formula
            let numerator = wale_tensor_mag_sq * wale_tensor_mag_sq.sqrt(); // (S^d_ij S^d_ij)^{3/2}
            let denominator_base = strain_mag_sq * strain_mag_sq.sqrt(); // (|S|²)^{5/2}
            let wale_term = wale_tensor_mag_sq * wale_tensor_mag_sq.sqrt().sqrt(); // (S^d_ij S^d_ij)^{5/4}

            let denominator = denominator_base + wale_term;

            if denominator > epsilon {
                c_w_delta_sq * numerator / denominator
            } else {
                T::zero()
            }
        }
    }

    /// Compute velocity gradients using central and one-sided differences.
    ///
    /// # Theorem -- Boundary-Consistent Differentiation
    ///
    /// For a uniform mesh and a component value `f_i`, the interior stencil
    /// `(f_{i+1}-f_{i-1})/(2h)` and the boundary stencils
    /// `(-3f_0+4f_1-f_2)/(2h)`, `(3f_n-4f_{n-1}+f_{n-2})/(2h)` are exact for
    /// all quadratic polynomials up to the derivative truncation term
    /// `O(h^2)`.
    ///
    /// **Proof.** Taylor expansion about the evaluation point gives
    /// `f(x+h)=f+h f'+h^2 f''/2+h^3 f'''/6+O(h^4)` and analogous terms for
    /// `x-h` and `x+2h`. Substituting the coefficients cancels the constant
    /// and second-derivative terms while leaving one unit of `f'`; the first
    /// uncancelled term is proportional to `h^2 f'''`. Therefore boundary
    /// gradients preserve the same second-order accuracy as the interior
    /// stencil for smooth velocity fields instead of imposing an artificial
    /// zero-gradient state.
    fn velocity_gradients(
        &self,
        velocity: &Field2D<Vector2<T>>,
        i: usize,
        j: usize,
        dx: T,
        dy: T,
    ) -> (T, T, T, T) {
        let nx = velocity.nx();
        let ny = velocity.ny();

        let du_dx = Self::differentiate_uniform(nx, i, dx, |ii| velocity.at(ii, j).x);
        let du_dy = Self::differentiate_uniform(ny, j, dy, |jj| velocity.at(i, jj).x);
        let dv_dx = Self::differentiate_uniform(nx, i, dx, |ii| velocity.at(ii, j).y);
        let dv_dy = Self::differentiate_uniform(ny, j, dy, |jj| velocity.at(i, jj).y);

        (du_dx, du_dy, dv_dx, dv_dy)
    }

    #[inline]
    fn differentiate_uniform<F>(n: usize, idx: usize, spacing: T, sample: F) -> T
    where
        F: Fn(usize) -> T,
    {
        let two = T::from_f64(2.0).expect("analytical constant conversion");
        if n < 2 {
            return T::zero();
        }
        if n == 2 {
            return (sample(1) - sample(0)) / spacing;
        }
        if idx == 0 {
            return (-T::from_f64(3.0).expect("analytical constant conversion") * sample(0)
                + T::from_f64(4.0).expect("analytical constant conversion") * sample(1)
                - sample(2))
                / (two * spacing);
        }
        if idx + 1 == n {
            return (T::from_f64(3.0).expect("analytical constant conversion") * sample(n - 1)
                - T::from_f64(4.0).expect("analytical constant conversion") * sample(n - 2)
                + sample(n - 3))
                / (two * spacing);
        }
        (sample(idx + 1) - sample(idx - 1)) / (two * spacing)
    }

    /// Compute WALE tensor magnitude squared
    ///
    /// The WALE tensor S_ij^d is the traceless symmetric part of the square
    /// of the velocity gradient tensor. Its magnitude squared is computed as:
    ///
    /// |S^d|² = S_ij^d S_ij^d where
    /// S_ij^d = (1/2)(∂_i u_k ∂_k u_j + ∂_j u_k ∂_k u_i) - (1/3)δ_ij ∂_k u_l ∂_l u_k
    fn wale_tensor_magnitude_squared(&self, du_dx: T, du_dy: T, dv_dx: T, dv_dy: T) -> T {
        // Velocity gradient tensor G_ij = ∂u_i/∂x_j
        let g_xx = du_dx;
        let g_xy = du_dy;
        let g_yx = dv_dx;
        let g_yy = dv_dy;

        // Square of velocity gradient: G²_ij = G_ik G_kj
        let g2_xx = g_xx * g_xx + g_xy * g_yx;
        let g2_xy = g_xx * g_xy + g_xy * g_yy;
        let g2_yx = g_yx * g_xx + g_yy * g_yx;
        let g2_yy = g_yx * g_xy + g_yy * g_yy;

        // Symmetric part of G²: (1/2)(G²_ij + G²_ji)
        let s2_xx = g2_xx;
        let s2_xy = T::from_f64(0.5).expect("analytical constant conversion") * (g2_xy + g2_yx);
        let s2_yy = g2_yy;

        // Trace of symmetric G²: tr(S²) = S²_kk
        let trace_s2 = s2_xx + s2_yy;

        // WALE tensor: S_ij^d = S²_ij - (1/3)δ_ij tr(S²)
        let sd_xx =
            s2_xx - T::from_f64(1.0 / 3.0).expect("analytical constant conversion") * trace_s2;
        let sd_xy = s2_xy; // Off-diagonal terms unchanged
        let sd_yy =
            s2_yy - T::from_f64(1.0 / 3.0).expect("analytical constant conversion") * trace_s2;

        // Magnitude squared: S^d_ij S^d_ij
        sd_xx * sd_xx
            + T::from_f64(2.0).expect("analytical constant conversion") * sd_xy * sd_xy
            + sd_yy * sd_yy
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_wale_viscosity() {
        let model = WaleModel::<f64>::new();
        let mut velocity = Field2D::new(5, 5, Vector2::new(1.0, 0.0));

        // Set up a simple shear flow: u = y, v = 0
        for i in 0..5 {
            for j in 0..5 {
                let y = j as f64 * 0.1;
                velocity[(i, j)] = Vector2::new(y, 0.0);
            }
        }

        let nu_sgs = model.sgs_viscosity(&velocity, 2, 2, 0.1, 0.1, 0.1);

        // Should be positive but small for this simple case
        assert!(nu_sgs >= 0.0);
        assert!(nu_sgs < 1e-3); // Small value expected
    }

    #[test]
    fn test_wale_at_wall() {
        let model = WaleModel::<f64>::new();
        let mut velocity = Field2D::new(5, 5, Vector2::new(0.0, 0.0));

        // Set up no-slip boundary conditions (walls)
        for i in 0..5 {
            velocity[(i, 0)] = Vector2::new(0.0, 0.0); // Bottom wall
            velocity[(i, 4)] = Vector2::new(0.0, 0.0); // Top wall
        }

        let nu_sgs_wall = model.sgs_viscosity(&velocity, 2, 0, 0.1, 0.1, 0.1);

        // At walls, WALE tensor should vanish, giving zero viscosity
        assert_relative_eq!(nu_sgs_wall, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn boundary_gradients_use_second_order_one_sided_stencils() {
        let model = WaleModel::<f64>::new();
        let mut velocity = Field2D::new(5, 5, Vector2::new(0.0, 0.0));

        for i in 0..5 {
            for j in 0..5 {
                let x = i as f64 * 0.25;
                let y = j as f64 * 0.25;
                velocity[(i, j)] = Vector2::new(x * x + 2.0 * y, 3.0 * x + y * y);
            }
        }

        let (du_dx_left, du_dy_bottom, dv_dx_left, dv_dy_bottom) =
            model.velocity_gradients(&velocity, 0, 0, 0.25, 0.25);
        assert_relative_eq!(du_dx_left, 0.0, epsilon = 1.0e-12);
        assert_relative_eq!(du_dy_bottom, 2.0, epsilon = 1.0e-12);
        assert_relative_eq!(dv_dx_left, 3.0, epsilon = 1.0e-12);
        assert_relative_eq!(dv_dy_bottom, 0.0, epsilon = 1.0e-12);

        let (du_dx_right, _, _, dv_dy_top) = model.velocity_gradients(&velocity, 4, 4, 0.25, 0.25);
        assert_relative_eq!(du_dx_right, 2.0, epsilon = 1.0e-12);
        assert_relative_eq!(dv_dy_top, 2.0, epsilon = 1.0e-12);
    }
}
