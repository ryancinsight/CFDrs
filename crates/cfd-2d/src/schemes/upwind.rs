//! Upwind discretization schemes
//!
//! These schemes implement proper upwind discretization based on the velocity field,
//! not the scalar field being transported. This is critical for physical correctness.

use super::{FaceReconstruction, Grid2D, SpatialDiscretization};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// First-order upwind scheme
pub struct FirstOrderUpwind<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive + Copy> Default for FirstOrderUpwind<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> FirstOrderUpwind<T> {
    /// Create new first-order upwind scheme
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> FaceReconstruction<T> for FirstOrderUpwind<T> {
    fn reconstruct_face_value_x(
        &self,
        phi: &Grid2D<T>,
        velocity_at_face: T,
        i: usize,
        j: usize,
    ) -> T {
        // Check boundaries to prevent panics
        let nx = phi.data.ncols();
        if i >= nx - 1 {
            panic!(
                "Cannot apply upwind scheme at domain boundary i={}, nx={}",
                i, nx
            );
        }

        if velocity_at_face > T::zero() {
            // Positive velocity: flow is from left (i) to right (i+1)
            // The upwind value is from the cell on the left
            phi.data[(i, j)]
        } else {
            // Negative velocity: flow is from right (i+1) to left (i)
            // The upwind value is from the cell on the right
            phi.data[(i + 1, j)]
        }
    }

    fn reconstruct_face_value_y(
        &self,
        phi: &Grid2D<T>,
        velocity_at_face: T,
        i: usize,
        j: usize,
    ) -> T {
        // Check boundaries to prevent panics
        let ny = phi.data.nrows();
        if j >= ny - 1 {
            panic!(
                "Cannot apply upwind scheme at domain boundary j={}, ny={}",
                j, ny
            );
        }

        if velocity_at_face > T::zero() {
            // Positive velocity: flow is from bottom (j) to top (j+1)
            // The upwind value is from the cell on the bottom
            phi.data[(i, j)]
        } else {
            // Negative velocity: flow is from top (j+1) to bottom (j)
            // The upwind value is from the cell on the top
            phi.data[(i, j + 1)]
        }
    }

    fn order(&self) -> usize {
        1
    }
}

// Keep backward compatibility with old interface (deprecated)
impl<T: RealField + Copy + FromPrimitive + Copy> SpatialDiscretization<T> for FirstOrderUpwind<T> {
    fn compute_derivative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        // This old interface is fundamentally flawed as it doesn't have access to velocity
        // We mark it as deprecated and panic if used
        panic!("The compute_derivative interface is physically incorrect for upwind schemes. Use FaceReconstruction trait instead.");
    }

    fn order(&self) -> usize {
        1
    }

    fn is_conservative(&self) -> bool {
        true
    }
}

/// Second-order upwind scheme
pub struct SecondOrderUpwind<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive + Copy> Default for SecondOrderUpwind<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> SecondOrderUpwind<T> {
    /// Create new second-order upwind scheme
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> FaceReconstruction<T> for SecondOrderUpwind<T> {
    fn reconstruct_face_value_x(
        &self,
        phi: &Grid2D<T>,
        velocity_at_face: T,
        i: usize,
        j: usize,
    ) -> T {
        let half = T::from_f64(0.5).expect("Failed to represent 0.5 in numeric type T");
        let three_half = T::from_f64(1.5).expect("Failed to represent 1.5 in numeric type T");

        if velocity_at_face > T::zero() {
            // Positive velocity: use backward-biased stencil
            // Check boundaries
            if i < 1 {
                // Fall back to first-order at boundary
                return phi.data[(i, j)];
            }
            // Second-order reconstruction: phi_face = 3/2*phi_i - 1/2*phi_{i-1}
            three_half * phi.data[(i, j)] - half * phi.data[(i - 1, j)]
        } else {
            // Negative velocity: use forward-biased stencil
            // Check boundaries
            let nx = phi.data.ncols();
            if i >= nx - 2 {
                // Fall back to first-order at boundary
                return phi.data[(i + 1, j)];
            }
            // Second-order reconstruction: phi_face = 3/2*phi_{i+1} - 1/2*phi_{i+2}
            three_half * phi.data[(i + 1, j)] - half * phi.data[(i + 2, j)]
        }
    }

    fn reconstruct_face_value_y(
        &self,
        phi: &Grid2D<T>,
        velocity_at_face: T,
        i: usize,
        j: usize,
    ) -> T {
        let half = T::from_f64(0.5).expect("Failed to represent 0.5 in numeric type T");
        let three_half = T::from_f64(1.5).expect("Failed to represent 1.5 in numeric type T");

        if velocity_at_face > T::zero() {
            // Positive velocity: use backward-biased stencil
            // Check boundaries
            if j < 1 {
                // Fall back to first-order at boundary
                return phi.data[(i, j)];
            }
            // Second-order reconstruction: phi_face = 3/2*phi_j - 1/2*phi_{j-1}
            three_half * phi.data[(i, j)] - half * phi.data[(i, j - 1)]
        } else {
            // Negative velocity: use forward-biased stencil
            // Check boundaries
            let ny = phi.data.nrows();
            if j >= ny - 2 {
                // Fall back to first-order at boundary
                return phi.data[(i, j + 1)];
            }
            // Second-order reconstruction: phi_face = 3/2*phi_{j+1} - 1/2*phi_{j+2}
            three_half * phi.data[(i, j + 1)] - half * phi.data[(i, j + 2)]
        }
    }

    fn order(&self) -> usize {
        2
    }
}

// Keep backward compatibility with old interface (deprecated)
impl<T: RealField + Copy + FromPrimitive + Copy> SpatialDiscretization<T> for SecondOrderUpwind<T> {
    fn compute_derivative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        // This old interface is fundamentally flawed as it doesn't have access to velocity
        // We mark it as deprecated and panic if used
        panic!("The compute_derivative interface is physically incorrect for upwind schemes. Use FaceReconstruction trait instead.");
    }

    fn order(&self) -> usize {
        2
    }

    fn is_conservative(&self) -> bool {
        true
    }
}
