//! Symmetric and asymmetric trifurcation geometry for 2D CFD validation
//!
//! This module implements a branching channel geometry with three outlets
//! used to validate 2D trifurcation flow simulations against literature.

use super::{BoundaryCondition, BoundaryFace, Geometry2D, Point2D};
use nalgebra::RealField;

/// Trifurcation geometry consisting of one parent and three daughter branches
#[derive(Debug, Clone)]
pub struct Trifurcation2D<T: RealField> {
    /// Parent channel width
    pub parent_width: T,
    /// Parent channel length
    pub parent_length: T,
    /// Daughter 1 width (Center)
    pub daughter1_width: T,
    /// Daughter 1 length
    pub daughter1_length: T,
    /// Daughter 1 angle (radians from x-axis)
    pub daughter1_angle: T,
    /// Daughter 2 width (Top)
    pub daughter2_width: T,
    /// Daughter 2 length
    pub daughter2_length: T,
    /// Daughter 2 angle (radians from x-axis)
    pub daughter2_angle: T,
    /// Daughter 3 width (Bottom)
    pub daughter3_width: T,
    /// Daughter 3 length
    pub daughter3_length: T,
    /// Daughter 3 angle (radians from x-axis)
    pub daughter3_angle: T,
    /// Junction center point (where parent ends and daughters begin)
    pub junction_center: Point2D<T>,
}

impl<T: RealField + Copy> Trifurcation2D<T> {
    /// Create a new symmetric trifurcation
    pub fn new_symmetric(
        width: T,
        length: T,
        daughter_width: T,
        daughter_length: T,
        angle: T,
    ) -> Self {
        Self {
            parent_width: width,
            parent_length: length,
            daughter1_width: daughter_width,
            daughter1_length: daughter_length,
            daughter1_angle: T::zero(), // Center daughter is horizontal
            daughter2_width: daughter_width,
            daughter2_length: daughter_length,
            daughter2_angle: angle, // Top daughter
            daughter3_width: daughter_width,
            daughter3_length: daughter_length,
            daughter3_angle: -angle, // Bottom daughter
            junction_center: Point2D { x: length, y: T::zero() },
        }
    }

    /// Check if a point is inside a specific branch segment
    fn in_segment(
        point: &Point2D<T>,
        start: &Point2D<T>,
        angle: T,
        length: T,
        width: T,
    ) -> bool {
        // Translate point to local coordinates
        let dx = point.x - start.x;
        let dy = point.y - start.y;

        // Rotate to align with segment
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let local_x = dx * cos_a + dy * sin_a;
        let local_y = -dx * sin_a + dy * cos_a;

        let half_width = width / (T::one() + T::one());

        local_x >= T::zero() && local_x <= length && local_y.abs() <= half_width
    }
}

impl<T: RealField + Copy> Geometry2D<T> for Trifurcation2D<T> {
    fn clone_box(&self) -> Box<dyn Geometry2D<T>> {
        Box::new(self.clone())
    }

    fn contains(&self, point: &Point2D<T>) -> bool {
        // Check parent branch (horizontal from 0 to parent_length)
        let parent_start = Point2D { x: T::zero(), y: T::zero() };
        if Self::in_segment(point, &parent_start, T::zero(), self.parent_length, self.parent_width) {
            return true;
        }

        // Check daughter 1 (Center)
        if Self::in_segment(point, &self.junction_center, self.daughter1_angle, self.daughter1_length, self.daughter1_width) {
            return true;
        }

        // Check daughter 2 (Top)
        if Self::in_segment(point, &self.junction_center, self.daughter2_angle, self.daughter2_length, self.daughter2_width) {
            return true;
        }

        // Check daughter 3 (Bottom)
        if Self::in_segment(point, &self.junction_center, self.daughter3_angle, self.daughter3_length, self.daughter3_width) {
            return true;
        }

        false
    }

    fn distance_to_boundary(&self, _point: &Point2D<T>) -> T {
        T::zero()
    }

    fn boundary_normal(&self, _point: &Point2D<T>) -> Option<Point2D<T>> {
        None
    }

    fn boundary_condition(&self, face: BoundaryFace, _s: T) -> BoundaryCondition<T> {
        match face {
            BoundaryFace::Left => BoundaryCondition::Dirichlet(T::one()), // Inlet
            _ => BoundaryCondition::Dirichlet(T::zero()), // Walls/Outlets
        }
    }

    fn bounds(&self) -> (Point2D<T>, Point2D<T>) {
        // Calculate max bounds
        let x_max = self.parent_length + self.daughter1_length.max(
            self.daughter2_length * self.daughter2_angle.cos().max(self.daughter3_angle.cos())
        );
        let y_max = (self.daughter2_length * self.daughter2_angle.sin()).abs()
            .max((self.daughter3_length * self.daughter3_angle.sin()).abs())
            + self.parent_width;
        
        (
            Point2D { x: T::zero(), y: -y_max },
            Point2D { x: x_max, y: y_max },
        )
    }

    fn boundary_parameter(&self, _face: BoundaryFace, _point: &Point2D<T>) -> Option<T> {
        None
    }

    fn on_boundary(&self, _point: &Point2D<T>, _face: BoundaryFace, _tolerance: T) -> bool {
        false
    }

    fn measure(&self) -> T {
        // Area of the branches
        self.parent_width * self.parent_length 
            + self.daughter1_width * self.daughter1_length 
            + self.daughter2_width * self.daughter2_length
            + self.daughter3_width * self.daughter3_length
    }
}
