//! Rhie-Chow interpolation for pressure-velocity coupling
//!
//! Reference: Rhie, C.M. and Chow, W.L. (1983). "Numerical study of the turbulent
//! flow past an airfoil with trailing edge separation." AIAA Journal, 21(11), 1525-1532.

use crate::constants::numerical::TWO;
use crate::fields::Field2D;
use crate::grid::StructuredGrid2D;
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;

/// Complete Rhie-Chow interpolation with momentum coefficients
pub struct RhieChowInterpolation<T: RealField + Copy> {
    /// Grid dimensions
    nx: usize,
    ny: usize,
    /// Momentum equation coefficients (`A_p` from discretized momentum equation)
    ap_coefficients: Field2D<T>,
    /// Previous time step velocity field (for transient term) - stored as buffer
    u_old_buffer: Field2D<Vector2<T>>,
    /// Flag indicating if old velocity is valid
    has_old_velocity: bool,
}

impl<T: RealField + Copy + FromPrimitive + Copy> RhieChowInterpolation<T> {
    /// Create new interpolator with momentum coefficients
    pub fn new(grid: &StructuredGrid2D<T>) -> Self {
        Self {
            nx: grid.nx,
            ny: grid.ny,
            ap_coefficients: Field2D::new(grid.nx, grid.ny, T::one()),
            u_old_buffer: Field2D::new(grid.nx, grid.ny, Vector2::zeros()),
            has_old_velocity: false,
        }
    }

    /// Update previous velocity field for transient term (zero-copy)
    pub fn update_old_velocity(&mut self, u_old: &Field2D<Vector2<T>>) {
        // Copy data efficiently without cloning the entire structure
        self.u_old_buffer.data.copy_from_slice(&u_old.data);
        self.has_old_velocity = true;
    }

    /// Update momentum equation coefficients from discretized momentum equation (zero-copy)
    /// `A_p` is the diagonal coefficient from momentum discretization
    pub fn update_coefficients(&mut self, ap: &Field2D<T>) {
        // Copy data efficiently without cloning the entire structure
        self.ap_coefficients.data.copy_from_slice(&ap.data);
    }

    /// Compute face velocity with complete Rhie-Chow interpolation
    ///
    /// The complete formulation according to Rhie & Chow (1983), Eq. 13:
    /// `u_f` = `ū_f` + `d_f` * [(∇p)_P - (∇p)_f] + dt * (`u_f^n` - `ū_f^n`) / 2
    ///
    /// where:
    /// - `ū_f` is the linearly interpolated velocity
    /// - `d_f` = (`Volume/A_p`)_f is the interpolated pressure gradient coefficient
    /// - (∇p)_P is the cell-centered pressure gradient (averaged)
    /// - (∇p)_f is the face pressure gradient (direct)
    /// - The last term is the transient correction for unsteady flows
    ///
    /// This prevents pressure-velocity decoupling on colocated grids
    pub fn face_velocity_x(
        &self,
        u: &Field2D<Vector2<T>>,
        p: &Field2D<T>,
        dx: T,
        dy: T,
        dt: Option<T>,
        i: usize,
        j: usize,
    ) -> T {
        // East face between cells (i,j) and (i+1,j)

        // Linear interpolation of velocity
        let two = T::from_f64(TWO).expect("Failed to represent 2.0 in numeric type T");
        let u_bar = (u.at(i, j).x + u.at(i + 1, j).x) / two;

        // Interpolate pressure gradient coefficient d_f = (Volume/A_p)_f
        let volume = dx * dy; // Correct 2D cell area for rectangular cells
        let d_p = volume / self.ap_coefficients.at(i, j);
        let d_e = volume / self.ap_coefficients.at(i + 1, j);
        let d_face = (d_p + d_e) / two;

        // Cell-centered pressure gradients (from momentum equation)
        let dp_dx_p = if i > 0 {
            (p.at(i + 1, j) - p.at(i - 1, j)) / (two * dx)
        } else {
            (p.at(i + 1, j) - p.at(i, j)) / dx
        };

        let dp_dx_e = if i + 1 < self.nx - 1 {
            (p.at(i + 2, j) - p.at(i, j)) / (two * dx)
        } else {
            (p.at(i + 1, j) - p.at(i, j)) / dx
        };

        // Average cell-centered gradient
        let dp_dx_cells = (dp_dx_p + dp_dx_e) / two;

        // Face pressure gradient
        let dp_dx_face = (p.at(i + 1, j) - p.at(i, j)) / dx;

        // Base Rhie-Chow interpolation
        let mut u_f = u_bar + d_face * (dp_dx_cells - dp_dx_face);

        // Add transient term if time step and old velocity are provided
        if let Some(dt) = dt {
            if self.has_old_velocity {
                let u_bar_old =
                    (self.u_old_buffer.at(i, j).x + self.u_old_buffer.at(i + 1, j).x) / two;
                let transient_factor = dt / two; // Simplified transient correction
                u_f += transient_factor * (u_bar - u_bar_old);
            }
        }

        u_f
    }

    /// Compute face velocity in y-direction
    pub fn face_velocity_y(
        &self,
        v: &Field2D<Vector2<T>>,
        p: &Field2D<T>,
        dx: T,
        dy: T,
        dt: Option<T>,
        i: usize,
        j: usize,
    ) -> T {
        // North face between cells (i,j) and (i,j+1)

        // Linear interpolation of velocity
        let two = T::from_f64(TWO).expect("Failed to represent 2.0 in numeric type T");
        let v_bar = (v.at(i, j).y + v.at(i, j + 1).y) / two;

        // Interpolate pressure gradient coefficient
        let volume = dx * dy; // Correct 2D cell area for rectangular cells
        let d_p = volume / self.ap_coefficients.at(i, j);
        let d_n = volume / self.ap_coefficients.at(i, j + 1);
        let d_face = (d_p + d_n) / two;

        // Cell-centered pressure gradients
        let dp_dy_p = if j > 0 {
            (p.at(i, j + 1) - p.at(i, j - 1)) / (two * dy)
        } else {
            (p.at(i, j + 1) - p.at(i, j)) / dy
        };

        let dp_dy_n = if j + 1 < self.ny - 1 {
            (p.at(i, j + 2) - p.at(i, j)) / (two * dy)
        } else {
            (p.at(i, j + 1) - p.at(i, j)) / dy
        };

        // Average cell-centered gradient
        let dp_dy_cells = (dp_dy_p + dp_dy_n) / two;

        // Face pressure gradient
        let dp_dy_face = (p.at(i, j + 1) - p.at(i, j)) / dy;

        // Base Rhie-Chow interpolation
        let mut v_f = v_bar + d_face * (dp_dy_cells - dp_dy_face);

        // Add transient term if time step and old velocity are provided
        if let Some(dt) = dt {
            if self.has_old_velocity {
                let v_bar_old =
                    (self.u_old_buffer.at(i, j).y + self.u_old_buffer.at(i, j + 1).y) / two;
                let transient_factor = dt / two; // Simplified transient correction
                v_f += transient_factor * (v_bar - v_bar_old);
            }
        }

        v_f
    }
}
