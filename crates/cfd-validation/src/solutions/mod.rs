//! Analytical solutions for CFD verification
//!
//! This module provides exact solutions to fundamental fluid dynamics problems
//! for verification of numerical methods.

use nalgebra::{RealField, Vector3};

/// Trait for analytical solutions
pub trait Solution<T: RealField + Copy> {
    /// Evaluate the solution at given coordinates and time
    fn evaluate(&self, x: T, y: T, z: T, t: T) -> Vector3<T>;

    /// Get the velocity field at given coordinates and time
    fn velocity(&self, x: T, y: T, z: T, t: T) -> Vector3<T> {
        self.evaluate(x, y, z, t)
    }

    /// Get the pressure field at given coordinates and time
    fn pressure(&self, x: T, y: T, z: T, t: T) -> T;

    /// Get the name of the solution
    fn name(&self) -> &str;

    /// Get the domain bounds [x_min, x_max, y_min, y_max, z_min, z_max]
    fn domain_bounds(&self) -> [T; 6];
}

// Import from analytical module
pub use crate::analytical::{CouetteFlow, PoiseuilleFlow, StokesFlow, TaylorGreenVortex};

// Re-export the trait
pub use self::Solution as AnalyticalSolution;
