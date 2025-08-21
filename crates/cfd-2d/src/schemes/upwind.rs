//! Upwind discretization schemes

use nalgebra::RealField;
use num_traits::FromPrimitive;
use super::{Grid2D, SpatialDiscretization};

/// First-order upwind scheme
pub struct FirstOrderUpwind<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + FromPrimitive + Copy> FirstOrderUpwind<T> {
    /// Create new first-order upwind scheme
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + FromPrimitive + Copy> SpatialDiscretization<T> for FirstOrderUpwind<T> {
    fn compute_derivative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        let velocity = grid.data[(i, j)];
        
        if velocity > T::zero() {
            // Backward difference for positive velocity
            (grid.data[(i, j)] - grid.data[(i - 1, j)]) / grid.dx
        } else {
            // Forward difference for negative velocity
            (grid.data[(i + 1, j)] - grid.data[(i, j)]) / grid.dx
        }
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

impl<T: RealField + FromPrimitive + Copy> SecondOrderUpwind<T> {
    /// Create new second-order upwind scheme
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + FromPrimitive + Copy> SpatialDiscretization<T> for SecondOrderUpwind<T> {
    fn compute_derivative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        let velocity = grid.data[(i, j)];
        let three = T::from_f64(3.0).unwrap_or_else(T::zero);
        let four = T::from_f64(4.0).unwrap_or_else(T::zero);
        let two = T::from_f64(2.0).unwrap_or_else(T::zero);
        
        if velocity > T::zero() {
            // Backward-biased for positive velocity
            (three * grid.data[(i, j)] - four * grid.data[(i - 1, j)] + grid.data[(i - 2, j)]) 
                / (two * grid.dx)
        } else {
            // Forward-biased for negative velocity
            (-grid.data[(i + 2, j)] + four * grid.data[(i + 1, j)] - three * grid.data[(i, j)]) 
                / (two * grid.dx)
        }
    }
    
    fn order(&self) -> usize {
        2
    }
    
    fn is_conservative(&self) -> bool {
        true
    }
}