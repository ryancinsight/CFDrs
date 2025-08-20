//! Central difference schemes

use nalgebra::RealField;
use num_traits::FromPrimitive;
use super::{Grid2D, SpatialDiscretization};

/// Second-order central difference scheme
pub struct CentralDifference<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> CentralDifference<T> {
    /// Create new central difference scheme
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + FromPrimitive + Copy> SpatialDiscretization<T> for CentralDifference<T> {
    fn compute_derivative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        let two = T::from_f64(2.0).unwrap_or_else(T::zero);
        (grid.data[(i + 1, j)] - grid.data[(i - 1, j)]) / (two * grid.dx)
    }
    
    fn order(&self) -> usize {
        2
    }
    
    fn is_conservative(&self) -> bool {
        true
    }
}

/// Fourth-order central difference scheme
pub struct FourthOrderCentral<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> FourthOrderCentral<T> {
    /// Create new fourth-order central scheme
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + FromPrimitive + Copy> SpatialDiscretization<T> for FourthOrderCentral<T> {
    fn compute_derivative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        let eight = T::from_f64(8.0).unwrap_or_else(T::zero);
        let twelve = T::from_f64(12.0).unwrap_or_else(T::zero);
        
        (-grid.data[(i + 2, j)] + eight * grid.data[(i + 1, j)] 
         - eight * grid.data[(i - 1, j)] + grid.data[(i - 2, j)]) 
         / (twelve * grid.dx)
    }
    
    fn order(&self) -> usize {
        4
    }
    
    fn is_conservative(&self) -> bool {
        true
    }
}