//! Rhie-Chow interpolation for pressure-velocity coupling
//!
//! Reference: Rhie, C.M. and Chow, W.L. (1983). "Numerical study of the turbulent
//! flow past an airfoil with trailing edge separation." AIAA Journal, 21(11), 1525-1532.

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
}
impl<T: RealField + Copy + FromPrimitive + Copy> RhieChowInterpolation<T> {
    /// Create new interpolator with momentum coefficients
    pub fn new(grid: &StructuredGrid2D<T>) -> Self {
        Self {
            nx: grid.nx,
            ny: grid.ny,
            ap_coefficients: Field2D::new(grid.nx, grid.ny, T::one()),
        }
    }
    /// Update momentum equation coefficients from discretized momentum equation
    /// `A_p` is the diagonal coefficient from momentum discretization
    pub fn update_coefficients(&mut self, ap: &Field2D<T>) {
        self.ap_coefficients = ap.clone();
    /// Compute face velocity with complete Rhie-Chow interpolation
    ///
    /// The complete formulation according to Rhie & Chow (1983), Eq. 13:
    /// `u_f` = `ū_f` + `d_f` * [(∇p)_P - (∇p)_f] + dt * (`u_f^n` - `ū_f^n`) / 2
    /// where:
    /// - `ū_f` is the linearly interpolated velocity
    /// - `d_f` = (`Volume/A_p`)_f is the interpolated pressure gradient coefficient
    /// - (∇p)_P is the cell-centered pressure gradient (averaged)
    /// - (∇p)_f is the face pressure gradient (direct)
    /// - The last term is the transient correction for unsteady flows
    /// This prevents pressure-velocity decoupling on colocated grids
    pub fn face_velocity_x(
        &self,
        u: &Field2D<Vector2<T>>,
        p: &Field2D<T>,
        dx: T,
        i: usize,
        j: usize,
    ) -> T {
        // East face between cells (i,j) and (i+1,j)
        // Linear interpolation of velocity
        let u_bar = (u.at(i, j).x + u.at(i + 1, j).x) / T::from_f64(2.0).unwrap_or_else(T::zero);
        // Interpolate pressure gradient coefficient d_f = (Volume/A_p)_f
        let volume = dx * dx; // 2D cell volume (assuming square cells)
        let d_p = volume / self.ap_coefficients.at(i, j);
        let d_e = volume / self.ap_coefficients.at(i + 1, j);
        let d_face = (d_p + d_e) / T::from_f64(2.0).unwrap_or_else(T::zero);
        // Cell-centered pressure gradients (from momentum equation)
        let dp_dx_p = if i > 0 {
            (p.at(i + 1, j) - p.at(i - 1, j)) / (T::from_f64(2.0).unwrap_or_else(T::zero) * dx)
        } else {
            (p.at(i + 1, j) - p.at(i, j)) / dx
        };
        let dp_dx_e = if i + 1 < self.nx - 1 {
            (p.at(i + 2, j) - p.at(i, j)) / (T::from_f64(2.0).unwrap_or_else(T::zero) * dx)
        // Average cell-centered gradient
        let dp_dx_cells = (dp_dx_p + dp_dx_e) / T::from_f64(2.0).unwrap_or_else(T::zero);
        // Face pressure gradient
        let dp_dx_face = (p.at(i + 1, j) - p.at(i, j)) / dx;
        // Complete Rhie-Chow interpolation
        u_bar + d_face * (dp_dx_cells - dp_dx_face)
    /// Compute face velocity in y-direction
    pub fn face_velocity_y(
        v: &Field2D<Vector2<T>>,
        dy: T,
        // North face between cells (i,j) and (i,j+1)
        let v_bar = (v.at(i, j).y + v.at(i, j + 1).y) / T::from_f64(2.0).unwrap_or_else(T::zero);
        // Interpolate pressure gradient coefficient
        let volume = dy * dy;
        let d_n = volume / self.ap_coefficients.at(i, j + 1);
        let d_face = (d_p + d_n) / T::from_f64(2.0).unwrap_or_else(T::zero);
        // Cell-centered pressure gradients
        let dp_dy_p = if j > 0 {
            (p.at(i, j + 1) - p.at(i, j - 1)) / (T::from_f64(2.0).unwrap_or_else(T::zero) * dy)
            (p.at(i, j + 1) - p.at(i, j)) / dy
        let dp_dy_n = if j + 1 < self.ny - 1 {
            (p.at(i, j + 2) - p.at(i, j)) / (T::from_f64(2.0).unwrap_or_else(T::zero) * dy)
        let dp_dy_cells = (dp_dy_p + dp_dy_n) / T::from_f64(2.0).unwrap_or_else(T::zero);
        let dp_dy_face = (p.at(i, j + 1) - p.at(i, j)) / dy;
        v_bar + d_face * (dp_dy_cells - dp_dy_face)
