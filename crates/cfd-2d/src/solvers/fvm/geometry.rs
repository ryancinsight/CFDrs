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

use cfd_core::error::Error;
use eunomia::{NumericElement, RealField};
use leto::geometry::Vector2;

/// Face between two control volumes
#[derive(Debug, Clone)]
pub struct Face<T: NumericElement> {
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

impl<T: NumericElement + RealField> Face<T> {
    /// Create a new face.
    ///
    /// # Panics
    /// Panics if `center` or `normal` has any non-finite component, `area`
    /// is non-finite or non-positive, or `normal` has zero magnitude (see
    /// [`Self::try_new`]).
    pub fn new(
        center: Vector2<T>,
        normal: Vector2<T>,
        area: T,
        owner: usize,
        neighbor: Option<usize>,
    ) -> Self {
        Self::try_new(center, normal, area, owner, neighbor).unwrap_or_else(|error| {
            panic!("Face::new called with invalid inputs: {error}");
        })
    }

    /// Create a new face with validation.
    ///
    /// # Errors
    /// Returns `Error::InvalidConfiguration` if:
    /// - any component of `center` is non-finite,
    /// - any component of `normal` is non-finite,
    /// - `area` is non-finite or non-positive,
    /// - `normal` has zero magnitude (cannot be normalized).
    pub fn try_new(
        center: Vector2<T>,
        normal: Vector2<T>,
        area: T,
        owner: usize,
        neighbor: Option<usize>,
    ) -> cfd_core::error::Result<Self> {
        if !<T as NumericElement>::is_finite(center.x)
            || !<T as NumericElement>::is_finite(center.y)
        {
            return Err(Error::InvalidConfiguration(format!(
                "Face::try_new: center components must be finite, got ({:?}, {:?})",
                center.x, center.y
            )));
        }
        if !<T as NumericElement>::is_finite(normal.x)
            || !<T as NumericElement>::is_finite(normal.y)
        {
            return Err(Error::InvalidConfiguration(format!(
                "Face::try_new: normal components must be finite, got ({:?}, {:?})",
                normal.x, normal.y
            )));
        }
        let n_mag_sq = normal.x * normal.x + normal.y * normal.y;
        if n_mag_sq <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidConfiguration(format!(
                "Face::try_new: normal must have positive magnitude, got |n|² = {n_mag_sq:?}"
            )));
        }
        if !<T as NumericElement>::is_finite(area) || area <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidConfiguration(format!(
                "Face::try_new: area must be finite and positive, got {area:?}"
            )));
        }
        Ok(Self {
            center,
            normal: normal.normalize(),
            area,
            owner,
            neighbor,
        })
    }

    /// Get the flux through this face
    pub fn flux(&self, velocity: Vector2<T>) -> T {
        self.area * velocity.dot(self.normal)
    }
}

impl<T: NumericElement> Face<T> {
    /// Check if this is a boundary face
    pub fn is_boundary(&self) -> bool {
        self.neighbor.is_none()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn abs(value: f64) -> f64 {
        <f64 as NumericElement>::abs(value)
    }

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
            abs(norm - 1.0) < 1e-14,
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
        assert!(abs(flux - 6.0) < 1e-14, "Expected 6.0, got {flux}");
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
            abs(flux) < 1e-14,
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
        assert!(abs(flux) < 1e-14, "Zero velocity flux should be zero");
    }

    /// **Positive**: `try_new` accepts valid inputs.
    #[test]
    fn face_try_new_accepts_valid_inputs() {
        let face = Face::<f64>::try_new(
            Vector2::new(0.5, 0.5),
            Vector2::new(1.0, 0.0),
            0.1,
            0,
            Some(1),
        )
        .expect("valid inputs must succeed");
        assert!(abs(face.area - 0.1) < 1e-14);
        assert!(face.normal.norm() - 1.0 < 1e-14);
    }

    /// **Adversarial**: non-finite `center` is rejected.
    #[test]
    fn face_try_new_rejects_non_finite_center() {
        match Face::<f64>::try_new(
            Vector2::new(f64::NAN, 0.5),
            Vector2::new(1.0, 0.0),
            0.1,
            0,
            None,
        ) {
            Err(e) => assert!(
                e.to_string().contains("center"),
                "error must mention center: {e}"
            ),
            Ok(_) => panic!("NaN center must be rejected"),
        }
    }

    /// **Adversarial**: non-finite `normal` is rejected.
    #[test]
    fn face_try_new_rejects_non_finite_normal() {
        match Face::<f64>::try_new(
            Vector2::new(0.5, 0.5),
            Vector2::new(f64::NAN, 0.0),
            0.1,
            0,
            None,
        ) {
            Err(e) => assert!(
                e.to_string().contains("normal"),
                "error must mention normal: {e}"
            ),
            Ok(_) => panic!("NaN normal must be rejected"),
        }
    }

    /// **Adversarial**: zero-magnitude `normal` is rejected.
    #[test]
    fn face_try_new_rejects_zero_normal() {
        match Face::<f64>::try_new(Vector2::new(0.5, 0.5), Vector2::new(0.0, 0.0), 0.1, 0, None) {
            Err(e) => assert!(
                e.to_string().contains("magnitude"),
                "error must mention magnitude: {e}"
            ),
            Ok(_) => panic!("zero normal must be rejected"),
        }
    }

    /// **Adversarial**: non-finite or non-positive `area` is rejected.
    #[test]
    fn face_try_new_rejects_invalid_area() {
        match Face::<f64>::try_new(Vector2::new(0.5, 0.5), Vector2::new(1.0, 0.0), 0.0, 0, None) {
            Err(e) => assert!(
                e.to_string().contains("area"),
                "error must mention area: {e}"
            ),
            Ok(_) => panic!("zero area must be rejected"),
        }
        match Face::<f64>::try_new(
            Vector2::new(0.5, 0.5),
            Vector2::new(1.0, 0.0),
            -0.1,
            0,
            None,
        ) {
            Err(e) => assert!(
                e.to_string().contains("area"),
                "error must mention area: {e}"
            ),
            Ok(_) => panic!("negative area must be rejected"),
        }
        match Face::<f64>::try_new(
            Vector2::new(0.5, 0.5),
            Vector2::new(1.0, 0.0),
            f64::NAN,
            0,
            None,
        ) {
            Err(e) => assert!(
                e.to_string().contains("area"),
                "error must mention area: {e}"
            ),
            Ok(_) => panic!("NaN area must be rejected"),
        }
    }

    /// **Boundary**: `new` panics on invalid `area` (thin wrapper contract).
    #[test]
    #[should_panic(expected = "area")]
    fn face_new_panics_on_invalid_area() {
        let _ = Face::<f64>::new(Vector2::new(0.5, 0.5), Vector2::new(1.0, 0.0), 0.0, 0, None);
    }
}
