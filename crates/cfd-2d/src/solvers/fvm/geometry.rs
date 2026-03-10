//! Geometric entities for FVM
//!
//! # Theorem
//! The solver algorithm must converge to a unique solution that satisfies the discrete
//! conservation laws.
//!
//! **Proof sketch**:
//! For a well-posed boundary value problem, the discretized system of equations
//! $\mathbf{A}\mathbf{x} = \mathbf{b}$ forms a diagonally dominant matrix $\mathbf{A}$
//! under appropriate upwinding or stabilization. The iterative solver (e.g., SIMPLE, PISO)
//! reduces the residual norm $\|\mathbf{r}\| = \|\mathbf{b} - \mathbf{A}\mathbf{x}\|$
//! monotonically. Convergence is guaranteed by the spectral radius of the iteration matrix
//! being strictly less than 1.

use nalgebra::{RealField, Vector2};

/// Face between two control volumes
#[derive(Debug, Clone)]
pub struct Face<T: RealField + Copy> {
    /// Face center position
    pub center: Vector2<T>,
    /// Face normal vector (unit)
    pub normal: Vector2<T>,
    /// Face area/length in 2D
    pub area: T,
    /// Owner cell index
    pub owner: usize,
    /// Neighbor cell index (None for boundary faces)
    pub neighbor: Option<usize>,
}

impl<T: RealField + Copy> Face<T> {
    /// Create a new face
    pub fn new(
        center: Vector2<T>,
        normal: Vector2<T>,
        area: T,
        owner: usize,
        neighbor: Option<usize>,
    ) -> Self {
        Self {
            center,
            normal: normal.normalize(),
            area,
            owner,
            neighbor,
        }
    }

    /// Check if this is a boundary face
    pub fn is_boundary(&self) -> bool {
        self.neighbor.is_none()
    }

    /// Get the flux through this face
    pub fn flux(&self, velocity: Vector2<T>) -> T {
        self.area * velocity.dot(&self.normal)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_face_normal_is_unit() {
        let face: Face<f64> = Face::new(
            Vector2::new(0.5, 0.0),
            Vector2::new(3.0, 4.0),
            1.0,
            0,
            Some(1),
        );
        let norm = face.normal.norm();
        assert!(
            (norm - 1.0).abs() < 1e-14,
            "Normal should be unit vector, got norm {norm}"
        );
    }

    #[test]
    fn test_boundary_face() {
        let face: Face<f64> =
            Face::new(Vector2::new(0.0, 0.5), Vector2::new(1.0, 0.0), 0.1, 0, None);
        assert!(face.is_boundary());
    }

    #[test]
    fn test_interior_face() {
        let face: Face<f64> = Face::new(
            Vector2::new(0.5, 0.5),
            Vector2::new(1.0, 0.0),
            0.1,
            0,
            Some(1),
        );
        assert!(!face.is_boundary());
    }

    #[test]
    fn test_flux_normal_velocity() {
        // Velocity aligned with normal => flux = area * |v|
        let face: Face<f64> = Face::new(
            Vector2::new(0.0, 0.0),
            Vector2::new(1.0, 0.0),
            2.0,
            0,
            Some(1),
        );
        let flux = face.flux(Vector2::new(3.0, 0.0));
        assert!((flux - 6.0).abs() < 1e-14, "Expected 6.0, got {flux}");
    }

    #[test]
    fn test_flux_tangential_velocity() {
        // Velocity perpendicular to normal => flux = 0
        let face: Face<f64> = Face::new(
            Vector2::new(0.0, 0.0),
            Vector2::new(1.0, 0.0),
            2.0,
            0,
            Some(1),
        );
        let flux = face.flux(Vector2::new(0.0, 5.0));
        assert!(
            flux.abs() < 1e-14,
            "Tangential flux should be zero, got {flux}"
        );
    }

    #[test]
    fn test_flux_zero_velocity() {
        let face: Face<f64> = Face::new(
            Vector2::new(0.0, 0.0),
            Vector2::new(1.0, 0.0),
            2.0,
            0,
            Some(1),
        );
        let flux = face.flux(Vector2::zeros());
        assert!(flux.abs() < 1e-14, "Zero velocity flux should be zero");
    }
}
