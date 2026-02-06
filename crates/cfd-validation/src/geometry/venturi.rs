//! Venturi geometry for 2D CFD validation
//!
//! This module implements a convergent-divergent nozzle geometry used to
//! validate 2D Venturi flow simulations against Bernoulli and cavitation models.

use super::{BoundaryCondition, BoundaryFace, Geometry2D, Point2D};
use nalgebra::RealField;

/// Venturi geometry with convergent, throat, and divergent sections
#[derive(Debug, Clone)]
pub struct Venturi2D<T: RealField> {
    /// Inlet width (m)
    pub inlet_width: T,
    /// Throat width (m)
    pub throat_width: T,
    /// Outlet width (m)
    pub outlet_width: T,
    /// Length of constant-width inlet section
    pub l_inlet: T,
    /// Length of convergent section
    pub l_converge: T,
    /// Length of constant-width throat section
    pub l_throat: T,
    /// Length of divergent section
    pub l_diverge: T,
    /// Length of constant-width outlet section
    pub l_outlet: T,
}

impl<T: RealField + Copy> Venturi2D<T> {
    /// Create a new standard Venturi geometry
    pub fn new(
        inlet_width: T,
        throat_width: T,
        converge_len: T,
        throat_len: T,
        diverge_len: T,
    ) -> Self {
        Self {
            inlet_width,
            throat_width,
            outlet_width: inlet_width, // Symmetric usually
            l_inlet: inlet_width * T::from_f64(2.0).unwrap(),
            l_converge: converge_len,
            l_throat: throat_len,
            l_diverge: diverge_len,
            l_outlet: inlet_width * T::from_f64(2.0).unwrap(),
        }
    }

    /// Get total length of the Venturi
    pub fn total_length(&self) -> T {
        self.l_inlet + self.l_converge + self.l_throat + self.l_diverge + self.l_outlet
    }

    /// Get local width at a given x-coordinate
    pub fn width_at(&self, x: T) -> T {
        if x < T::zero() {
            return self.inlet_width;
        }

        let mut current_x = T::zero();

        // Inlet
        if x < self.l_inlet {
            return self.inlet_width;
        }
        current_x += self.l_inlet;

        // Convergent
        if x < current_x + self.l_converge {
            let progress = (x - current_x) / self.l_converge;
            return self.inlet_width + progress * (self.throat_width - self.inlet_width);
        }
        current_x += self.l_converge;

        // Throat
        if x < current_x + self.l_throat {
            return self.throat_width;
        }
        current_x += self.l_throat;

        // Divergent
        if x < current_x + self.l_diverge {
            let progress = (x - current_x) / self.l_diverge;
            return self.throat_width + progress * (self.outlet_width - self.throat_width);
        }
        current_x += self.l_diverge;

        // Outlet
        self.outlet_width
    }
}

impl<T: RealField + Copy> Geometry2D<T> for Venturi2D<T> {
    fn clone_box(&self) -> Box<dyn Geometry2D<T>> {
        Box::new(self.clone())
    }

    fn contains(&self, point: &Point2D<T>) -> bool {
        if point.x < T::zero() || point.x > self.total_length() {
            return false;
        }

        let width = self.width_at(point.x);
        let half_width = width / (T::one() + T::one());

        point.y.abs() <= half_width
    }

    fn distance_to_boundary(&self, _point: &Point2D<T>) -> T {
        T::zero()
    }

    fn boundary_normal(&self, _point: &Point2D<T>) -> Option<Point2D<T>> {
        None
    }

    fn boundary_condition(&self, face: BoundaryFace, _s: T) -> BoundaryCondition<T> {
        match face {
            BoundaryFace::Left => BoundaryCondition::Dirichlet(T::one()), // Adjusted in benchmark
            _ => BoundaryCondition::Dirichlet(T::zero()),
        }
    }

    fn bounds(&self) -> (Point2D<T>, Point2D<T>) {
        let y_max = self.inlet_width.max(self.outlet_width) * T::from_f64(0.6).unwrap();
        (
            Point2D { x: T::zero(), y: -y_max },
            Point2D { x: self.total_length(), y: y_max },
        )
    }

    fn boundary_parameter(&self, _face: BoundaryFace, _point: &Point2D<T>) -> Option<T> {
        None
    }

    fn on_boundary(&self, _point: &Point2D<T>, _face: BoundaryFace, _tolerance: T) -> bool {
        false
    }

    fn measure(&self) -> T {
        // Integrate width over length
        let mut area = self.inlet_width * self.l_inlet;
        area += (self.inlet_width + self.throat_width) / (T::one() + T::one()) * self.l_converge;
        area += self.throat_width * self.l_throat;
        area += (self.throat_width + self.outlet_width) / (T::one() + T::one()) * self.l_diverge;
        area += self.outlet_width * self.l_outlet;
        area
    }
}
