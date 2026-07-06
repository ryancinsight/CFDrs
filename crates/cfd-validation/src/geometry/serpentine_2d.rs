//! Serpentine channel geometry for 2D CFD validation
//!
//! This module implements a sinusoidal or multi-bend serpentine channel
//! used to validate 2D/3D Dean vortex flows and mixing efficiency.

use super::{BoundaryCondition, BoundaryFace, Geometry2D, Point2D};
use crate::scalar;
use eunomia::{FloatElement, RealField};

/// Serpentine geometry consisting of multiple curved sections
#[derive(Debug, Clone)]
pub struct Serpentine2D<T: RealField> {
    /// Channel width (m)
    pub width: T,
    /// Amplitude of serpentine (m)
    pub amplitude: T,
    /// Number of periods
    pub periods: usize,
    /// Length of one period (m)
    pub period_length: T,
    /// Inlet length before serpentine (m)
    pub l_inlet: T,
    /// Outlet length after serpentine (m)
    pub l_outlet: T,
}

impl<T: RealField + Copy + FloatElement> Serpentine2D<T> {
    /// Create a new sinusoidal serpentine geometry
    pub fn new(width: T, amplitude: T, period_length: T, periods: usize) -> Self {
        Self {
            width,
            amplitude,
            periods,
            period_length,
            l_inlet: width * scalar::from_f64::<T>(5.0),
            l_outlet: width * scalar::from_f64::<T>(10.0),
        }
    }

    /// Get channel width
    #[must_use]
    pub fn width(&self) -> T {
        self.width
    }

    /// Get amplitude
    #[must_use]
    pub fn amplitude(&self) -> T {
        self.amplitude
    }

    /// Get period length
    #[must_use]
    pub fn period_length(&self) -> T {
        self.period_length
    }

    /// Get total length along x-axis
    pub fn total_length(&self) -> T {
        self.l_inlet + self.period_length * scalar::from_usize::<T>(self.periods) + self.l_outlet
    }

    /// Get centerline y-coordinate at given x
    pub fn centerline_y(&self, x: T) -> T {
        if x < self.l_inlet {
            return scalar::zero();
        }

        let serpentine_start = self.l_inlet;
        let serpentine_end =
            self.l_inlet + self.period_length * scalar::from_usize::<T>(self.periods);

        if x > serpentine_end {
            return scalar::zero();
        }

        let local_x = x - serpentine_start;
        let two_pi = scalar::from_f64::<T>(2.0 * std::f64::consts::PI);
        self.amplitude * scalar::sin(two_pi * local_x / self.period_length)
    }
}

impl<T: RealField + Copy + FloatElement> Geometry2D<T> for Serpentine2D<T> {
    fn clone_box(&self) -> Box<dyn Geometry2D<T>> {
        Box::new(self.clone())
    }

    fn contains(&self, point: &Point2D<T>) -> bool {
        if point.x < scalar::zero() || point.x > self.total_length() {
            return false;
        }

        let cy = self.centerline_y(point.x);
        let dy = point.y - cy;
        let half_width = self.width / (scalar::one::<T>() + scalar::one::<T>());

        scalar::abs(dy) <= half_width
    }

    fn distance_to_boundary(&self, _point: &Point2D<T>) -> T {
        scalar::zero()
    }

    fn boundary_normal(&self, _point: &Point2D<T>) -> Option<Point2D<T>> {
        None
    }

    fn boundary_condition(&self, face: BoundaryFace, _s: T) -> BoundaryCondition<T> {
        match face {
            BoundaryFace::Left => BoundaryCondition::Dirichlet(scalar::one()),
            _ => BoundaryCondition::Dirichlet(scalar::zero()),
        }
    }

    fn bounds(&self) -> (Point2D<T>, Point2D<T>) {
        let y_max = self.amplitude + self.width;
        (
            Point2D::new(scalar::zero(), -y_max),
            Point2D::new(self.total_length(), y_max),
        )
    }

    fn boundary_parameter(&self, _face: BoundaryFace, _point: &Point2D<T>) -> Option<T> {
        None
    }

    fn on_boundary(&self, _point: &Point2D<T>, _face: BoundaryFace, _tolerance: T) -> bool {
        false
    }

    fn measure(&self) -> T {
        // The swept area of a constant-width channel is width × arc_length.
        //
        // For a sinusoidal centerline  y(x) = A sin(2π x / λ),
        //   dy/dx = A (2π/λ) cos(2π x / λ)
        //
        // Arc length of one period:
        //   L_period = ∫₀^λ √(1 + (dy/dx)²) dx
        //
        // We compute this numerically with the trapezoidal rule (512 points
        // per period gives <1e-8 relative error for typical A/λ < 1).
        //
        // Reference: Weisstein, Eric W. "Arc Length." MathWorld.

        let two_pi = scalar::from_f64::<T>(2.0 * std::f64::consts::PI);
        let k = two_pi / self.period_length; // wave number
        let ak = self.amplitude * k; // A·k = A·2π/λ

        // Compute arc length for one period via composite trapezoidal rule
        let n_quad: usize = 512;
        let dx = self.period_length / scalar::from_usize::<T>(n_quad);
        let mut arc_length_period: T = scalar::zero();
        for i in 0..=n_quad {
            let x_loc = dx * scalar::from_usize::<T>(i);
            let dydx = ak * scalar::cos(k * x_loc);
            let integrand = scalar::sqrt(scalar::one::<T>() + dydx * dydx);
            let weight = if i == 0 || i == n_quad {
                scalar::from_f64::<T>(0.5)
            } else {
                scalar::one()
            };
            arc_length_period += weight * integrand * dx;
        }

        // Total arc length = inlet + serpentine periods + outlet
        // Inlet and outlet are straight, so their arc length equals their length.
        let arc_length_total = self.l_inlet
            + arc_length_period * scalar::from_usize::<T>(self.periods)
            + self.l_outlet;

        self.width * arc_length_total
    }
}
