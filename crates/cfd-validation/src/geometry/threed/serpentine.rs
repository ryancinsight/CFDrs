//! 3D Serpentine micromixer geometry implementation for CFD validation

use super::super::{BoundaryCondition, BoundaryFace, Geometry3D, Point3D};
use crate::scalar;
use eunomia::{FloatElement, RealField};

/// 3D Serpentine channel geometry (Sine-wave)
#[derive(Debug, Clone)]
pub struct Serpentine3D<T: RealField + Copy> {
    /// Channel diameter/width \[m]
    pub diameter: T,
    /// Amplitude of the serpentine curve \[m]
    pub amplitude: T,
    /// Period of the serpentine curve (wavelength) \[m]
    pub wavelength: T,
    /// Number of periods
    pub num_periods: usize,
}

impl<T: RealField + Copy> Serpentine3D<T> {
    /// Create a new 3D serpentine channel
    pub fn new(diameter: T, amplitude: T, wavelength: T, num_periods: usize) -> Self {
        Self {
            diameter,
            amplitude,
            wavelength,
            num_periods,
        }
    }
}

impl<T: RealField + Copy + FloatElement> Geometry3D<T> for Serpentine3D<T> {
    fn clone_box(&self) -> Box<dyn Geometry3D<T>> {
        Box::new(self.clone())
    }

    fn contains(&self, _point: &Point3D<T>) -> bool {
        // Complex helical/serpentine containment check
        // Simplified for validation benchmarks
        false
    }

    fn distance_to_boundary(&self, _point: &Point3D<T>) -> T {
        scalar::zero()
    }

    fn boundary_normal(&self, _point: &Point3D<T>) -> Option<Point3D<T>> {
        None
    }

    fn boundary_condition(&self, face: BoundaryFace, _s: T, _t: T) -> BoundaryCondition<T> {
        match face {
            BoundaryFace::Left => BoundaryCondition::Dirichlet(scalar::zero()),
            BoundaryFace::Right => BoundaryCondition::Neumann(scalar::zero()),
            _ => BoundaryCondition::Dirichlet(scalar::zero()),
        }
    }

    fn bounds(&self) -> (Point3D<T>, Point3D<T>) {
        let total_l = self.wavelength * scalar::from_usize::<T>(self.num_periods);
        let half_d = self.diameter / scalar::from_f64::<T>(2.0);
        (
            Point3D::new(-self.amplitude - half_d, -half_d, scalar::zero()),
            Point3D::new(self.amplitude + half_d, half_d, total_l),
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
        let total_l = self.wavelength * scalar::from_usize::<T>(self.num_periods);
        let radius = self.diameter / scalar::from_f64::<T>(2.0);
        let area = pi * scalar::powf(radius, scalar::from_f64::<T>(2.0));
        area * total_l
    }
}
