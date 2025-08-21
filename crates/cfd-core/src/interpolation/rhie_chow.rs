//! Rhie-Chow momentum interpolation
//!
//! Reference: Rhie, C.M. and Chow, W.L. (1983). "Numerical study of the turbulent 
//! flow past an airfoil with trailing edge separation". AIAA Journal, 21(11), 1525-1532.

use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Rhie-Chow interpolation for pressure-velocity coupling
/// 
/// This interpolation prevents checkerboard pressure oscillations in collocated grids
/// by adding a pressure-gradient-based correction to the interpolated velocity.
pub struct RhieChowInterpolation<T: RealField> {
    /// Grid spacing in x-direction
    dx: T,
    /// Grid spacing in y-direction
    dy: T,
    /// Relaxation factor
    alpha: T,
}

impl<T: RealField + FromPrimitive + Copy> RhieChowInterpolation<T> {
    /// Create new Rhie-Chow interpolator
    pub fn new(dx: T, dy: T) -> Self {
        Self {
            dx,
            dy,
            alpha: T::from_f64(0.8).unwrap_or_else(T::one),
        }
    }
    
    /// Interpolate velocity to face with pressure correction
    /// 
    /// u_f = ū_f - D_f * (∇p)_f
    /// 
    /// where:
    /// - u_f is the face velocity
    /// - ū_f is the linearly interpolated velocity
    /// - D_f is the face diffusion coefficient
    /// - (∇p)_f is the pressure gradient at the face
    pub fn interpolate_u_face(
        &self,
        u_p: T,      // velocity at cell P
        u_e: T,      // velocity at cell E
        p_p: T,      // pressure at cell P
        p_e: T,      // pressure at cell E
        a_p: T,      // diagonal coefficient for cell P
        a_e: T,      // diagonal coefficient for cell E
        dt: T,       // time step
    ) -> T {
        // Linear interpolation of velocity
        let u_bar = (u_p + u_e) / T::from_f64(2.0).unwrap();
        
        // Face diffusion coefficient (harmonic mean)
        let d_f = dt * T::from_f64(2.0).unwrap() / (a_p + a_e);
        
        // Pressure gradient at face
        let dp_dx = (p_e - p_p) / self.dx;
        
        // Rhie-Chow correction
        u_bar - d_f * dp_dx
    }
    
    /// Interpolate v-velocity to face with pressure correction
    pub fn interpolate_v_face(
        &self,
        v_p: T,      // velocity at cell P
        v_n: T,      // velocity at cell N
        p_p: T,      // pressure at cell P
        p_n: T,      // pressure at cell N
        a_p: T,      // diagonal coefficient for cell P
        a_n: T,      // diagonal coefficient for cell N
        dt: T,       // time step
    ) -> T {
        // Linear interpolation of velocity
        let v_bar = (v_p + v_n) / T::from_f64(2.0).unwrap();
        
        // Face diffusion coefficient
        let d_f = dt * T::from_f64(2.0).unwrap() / (a_p + a_n);
        
        // Pressure gradient at face
        let dp_dy = (p_n - p_p) / self.dy;
        
        // Rhie-Chow correction
        v_bar - d_f * dp_dy
    }
}