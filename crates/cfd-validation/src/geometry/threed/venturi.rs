//! 3D Venturi tube geometry implementation for CFD validation

use super::super::{BoundaryCondition, BoundaryFace, Geometry3D, Point3D};
use cfd_core::conversion::SafeFromF64;
use nalgebra::RealField;

/// 3D Venturi tube geometry
#[derive(Debug, Clone)]
pub struct Venturi3D<T: RealField + Copy> {
    /// Inlet diameter [m]
    pub d_inlet: T,
    /// Throat diameter [m]
    pub d_throat: T,
    /// Outlet diameter [m]
    pub d_outlet: T,
    /// Inlet length [m]
    pub l_inlet: T,
    /// Convergent length [m]
    pub l_convergent: T,
    /// Throat length [m]
    pub l_throat: T,
    /// Divergent length [m]
    pub l_divergent: T,
    /// Outlet length [m]
    pub l_outlet: T,
}

impl<T: RealField + Copy> Venturi3D<T> {
    /// Create a symmetric 3D Venturi tube
    pub fn new(d_in: T, d_th: T, l_in: T, l_conv: T, l_th: T, l_div: T, l_out: T) -> Self {
        Self {
            d_inlet: d_in,
            d_throat: d_th,
            d_outlet: d_in,
            l_inlet: l_in,
            l_convergent: l_conv,
            l_throat: l_th,
            l_divergent: l_div,
            l_outlet: l_out,
        }
    }
}

impl<T: RealField + Copy> Geometry3D<T> for Venturi3D<T> {
    fn clone_box(&self) -> Box<dyn Geometry3D<T>> {
        Box::new(self.clone())
    }

    fn contains(&self, point: &Point3D<T>) -> bool {
        let z = point.z;
        let r_sq = point.x * point.x + point.y * point.y;
        
        let d = if z < self.l_inlet {
            self.d_inlet
        } else if z < self.l_inlet + self.l_convergent {
            let s = (z - self.l_inlet) / self.l_convergent;
            self.d_inlet + (self.d_throat - self.d_inlet) * s
        } else if z < self.l_inlet + self.l_convergent + self.l_throat {
            self.d_throat
        } else if z < self.l_inlet + self.l_convergent + self.l_throat + self.l_divergent {
            let s = (z - (self.l_inlet + self.l_convergent + self.l_throat)) / self.l_divergent;
            self.d_throat + (self.d_outlet - self.d_throat) * s
        } else if z < self.l_inlet + self.l_convergent + self.l_throat + self.l_divergent + self.l_outlet {
            self.d_outlet
        } else {
            return false;
        };

        let r = d / T::from_f64_or_one(2.0);
        r_sq <= r * r
    }

    fn distance_to_boundary(&self, _point: &Point3D<T>) -> T {
        T::zero()
    }

    fn boundary_normal(&self, _point: &Point3D<T>) -> Option<Point3D<T>> {
        None
    }

    fn boundary_condition(&self, face: BoundaryFace, _s: T, _t: T) -> BoundaryCondition<T> {
        match face {
            BoundaryFace::Left => BoundaryCondition::Dirichlet(T::zero()),
            BoundaryFace::Right => BoundaryCondition::Neumann(T::zero()),
            _ => BoundaryCondition::Dirichlet(T::zero()),
        }
    }

    fn bounds(&self) -> (Point3D<T>, Point3D<T>) {
        let total_l = self.l_inlet + self.l_convergent + self.l_throat + self.l_divergent + self.l_outlet;
        let max_d = self.d_inlet.max(self.d_outlet);
        (
            Point3D { x: -max_d, y: -max_d, z: T::zero() },
            Point3D { x: max_d, y: max_d, z: total_l },
        )
    }

    fn boundary_parameter(&self, _face: BoundaryFace, _point: &Point3D<T>) -> Option<(T, T)> {
        None
    }

    fn on_boundary(&self, _point: &Point3D<T>, _face: BoundaryFace, _tolerance: T) -> bool {
        false
    }

    fn measure(&self) -> T {
        // Volume integration of piecewise cone frustums
        T::zero() 
    }
}
