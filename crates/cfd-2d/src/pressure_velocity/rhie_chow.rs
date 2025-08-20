//! Rhie-Chow interpolation for pressure-velocity decoupling
//!
//! Reference: Rhie, C.M. and Chow, W.L. (1983). "Numerical study of the turbulent
//! flow past an airfoil with trailing edge separation." AIAA Journal, 21(11), 1525-1532.

use nalgebra::{RealField, Vector2};
use crate::grid::StructuredGrid2D;

/// Rhie-Chow interpolation for colocated grids
pub struct RhieChowInterpolation<T: RealField> {
    /// Grid dimensions
    nx: usize,
    ny: usize,
}

impl<T: RealField> RhieChowInterpolation<T> {
    /// Create new interpolator
    pub fn new(grid: &StructuredGrid2D<T>) -> Self {
        Self {
            nx: grid.nx,
            ny: grid.ny,
        }
    }
    
    /// Interpolate face velocity with pressure gradient correction
    /// 
    /// u_f = ū_f - D_f * (∇p)_f
    /// 
    /// where:
    /// - ū_f is the linearly interpolated velocity
    /// - D_f is the interpolated pressure gradient coefficient
    /// - (∇p)_f is the pressure gradient at the face
    pub fn face_velocity(
        &self,
        u: &Vec<Vec<Vector2<T>>>,
        p: &Vec<Vec<T>>,
        d: &Vec<Vec<T>>,
        dx: T,
        dy: T,
        i: usize,
        j: usize,
        face: Face,
    ) -> T {
        match face {
            Face::East => {
                // East face between cells (i,j) and (i+1,j)
                let u_bar = (u[i][j].x + u[i+1][j].x) / (T::one() + T::one());
                let d_face = (d[i][j] + d[i+1][j]) / (T::one() + T::one());
                let dp_dx = (p[i+1][j] - p[i][j]) / dx;
                u_bar - d_face * dp_dx
            }
            Face::North => {
                // North face between cells (i,j) and (i,j+1)
                let v_bar = (u[i][j].y + u[i][j+1].y) / (T::one() + T::one());
                let d_face = (d[i][j] + d[i][j+1]) / (T::one() + T::one());
                let dp_dy = (p[i][j+1] - p[i][j]) / dy;
                v_bar - d_face * dp_dy
            }
            _ => T::zero(),
        }
    }
}

/// Cell face enumeration
pub enum Face {
    East,
    West,
    North,
    South,
}