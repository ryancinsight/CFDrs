//! 3D Bifurcation geometry implementation for CFD validation
//!
//! Models a 3D branching vessel with cylindrical parent and daughter branches.

use super::super::{BoundaryCondition, BoundaryFace, Geometry3D, Point3D};
use crate::scalar;
use eunomia::{FloatElement, RealField};
use leto::geometry::Vector3;

/// 3D Bifurcation geometry
#[derive(Debug, Clone)]
pub struct Bifurcation3D<T: RealField + Copy> {
    /// Parent branch diameter \[m]
    pub d_parent: T,
    /// Parent branch length \[m]
    pub l_parent: T,
    /// Daughter branch diameters \[m]
    pub d_daughters: [T; 2],
    /// Daughter branch lengths \[m]
    pub l_daughters: [T; 2],
    /// Branching angles (from parent axis) \[radians]
    pub angles: [T; 2],
}

impl<T: RealField + Copy> Bifurcation3D<T> {
    /// Create a symmetric 3D bifurcation
    pub fn symmetric(d_parent: T, d_daughter: T, l_parent: T, l_daughter: T, angle: T) -> Self {
        Self {
            d_parent,
            l_parent,
            d_daughters: [d_daughter, d_daughter],
            l_daughters: [l_daughter, l_daughter],
            angles: [angle, angle],
        }
    }
}

impl<T: RealField + Copy + FloatElement> Geometry3D<T> for Bifurcation3D<T> {
    fn clone_box(&self) -> Box<dyn Geometry3D<T>> {
        Box::new(self.clone())
    }

    fn contains(&self, point: &Point3D<T>) -> bool {
        // Simplified 3D containment check for validation benchmarks
        // In practice, this uses cylinder intersection tests

        // Parent branch (assumed along Z-axis from Z=0 to Z=l_parent)
        let half = scalar::from_f64::<T>(0.5);
        let r_p = self.d_parent * half;
        if point.z >= scalar::zero() && point.z <= self.l_parent {
            let r_sq = point.x * point.x + point.y * point.y;
            if r_sq <= r_p * r_p {
                return true;
            }
        }

        // Daughter branches starting from (0, 0, l_parent)
        for i in 0..2 {
            let angle = if i == 0 {
                self.angles[0]
            } else {
                -self.angles[1]
            };
            let d_rel = Vector3::new(point.x, point.y, point.z - self.l_parent);

            // Rotate back to daughter axis (simplified 2D rotation in XZ plane)
            let cos_a = scalar::cos(angle);
            let sin_a = scalar::sin(angle);
            let x_rot = d_rel.x * cos_a + d_rel.z * sin_a;
            let z_rot = -d_rel.x * sin_a + d_rel.z * cos_a;

            let r_d = self.d_daughters[i] * half;
            if z_rot >= scalar::zero() && z_rot <= self.l_daughters[i] {
                let r_sq = x_rot * x_rot + point.y * point.y;
                if r_sq <= r_d * r_d {
                    return true;
                }
            }
        }

        false
    }

    fn distance_to_boundary(&self, _point: &Point3D<T>) -> T {
        // Simplified for now, can be implemented with cylinder distance formulas
        scalar::zero()
    }

    fn boundary_normal(&self, _point: &Point3D<T>) -> Option<Point3D<T>> {
        None
    }

    fn boundary_condition(&self, face: BoundaryFace, _s: T, _t: T) -> BoundaryCondition<T> {
        match face {
            BoundaryFace::Left => BoundaryCondition::Dirichlet(scalar::zero()), // Inlet
            BoundaryFace::Right => BoundaryCondition::Neumann(scalar::zero()),  // Outlet
            _ => BoundaryCondition::Dirichlet(scalar::zero()),                  // Wall
        }
    }

    fn bounds(&self) -> (Point3D<T>, Point3D<T>) {
        // Rough bounding box
        let max_l = self.l_parent + scalar::max(self.l_daughters[0], self.l_daughters[1]);
        let max_d = scalar::max(
            scalar::max(self.d_parent, self.d_daughters[0]),
            self.d_daughters[1],
        );

        (
            Point3D::new(-max_d, -max_d, scalar::zero()),
            Point3D::new(max_d, max_d, max_l),
        )
    }

    fn boundary_parameter(&self, _face: BoundaryFace, _point: &Point3D<T>) -> Option<(T, T)> {
        None
    }

    fn on_boundary(&self, _point: &Point3D<T>, _face: BoundaryFace, _tolerance: T) -> bool {
        false
    }

    fn measure(&self) -> T {
        let pi = scalar::from_f64::<T>(std::f64::consts::PI);
        let r_parent = self.d_parent * scalar::from_f64::<T>(0.5);
        let r_daughter_1 = self.d_daughters[0] * scalar::from_f64::<T>(0.5);
        let r_daughter_2 = self.d_daughters[1] * scalar::from_f64::<T>(0.5);
        let v_p = pi * r_parent * r_parent * self.l_parent;
        let v_d1 = pi * r_daughter_1 * r_daughter_1 * self.l_daughters[0];
        let v_d2 = pi * r_daughter_2 * r_daughter_2 * self.l_daughters[1];
        v_p + v_d1 + v_d2
    }
}
