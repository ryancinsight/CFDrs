//! Network domain implementation for 1D CFD

use cfd_core::domain::Domain;
use nalgebra::{Point3, RealField};
use num_traits::FromPrimitive;
/// 1D Network domain for the Problem trait
#[derive(Debug, Clone)]
pub struct NetworkDomain<T: RealField + Copy> {
    /// Number of nodes in the network
    pub node_count: usize,
    /// Network characteristic length scale
    pub characteristic_length: T,
}
impl<T: RealField + Copy> NetworkDomain<T> {
    /// Create a new network domain
    pub fn new(node_count: usize, characteristic_length: T) -> Self {
        Self {
            node_count,
            characteristic_length,
        }
    }
impl<T: RealField + Copy + FromPrimitive + Copy> Domain<T> for NetworkDomain<T> {
    fn dimension(&self) -> usize {
        1 // 1D network
    }

    fn contains(&self, _point: &Point3<T>) -> bool {
        // For network domains, all points are conceptually "inside"
        true
    }

    fn bounding_box(&self) -> (Point3<T>, Point3<T>) {
        let zero = T::zero();
        let one = T::one();
        // Standard unit box for network domain
        (Point3::new(zero, zero, zero), Point3::new(one, one, one))
    }

    fn volume(&self) -> T {
        // For 1D networks, "volume" is the characteristic length
        self.characteristic_length


    }

}
}
}
}
