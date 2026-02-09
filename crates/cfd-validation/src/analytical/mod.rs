//! Analytical solutions for CFD validation
//!
//! This module provides exact analytical solutions for various fluid flow problems
//! that can be used to validate numerical CFD solvers.

pub mod blasius;
pub mod couette;
pub mod poiseuille;
pub mod poiseuille_2d;
pub mod stokes;
pub mod taylor_green;
pub mod utils;
pub mod womersley;

pub use blasius::BlasiusBoundaryLayer;
pub use couette::CouetteFlow;
pub use poiseuille::{PoiseuilleFlow, PoiseuilleGeometry};
pub use poiseuille_2d::{
    PowerLawPoiseuille, CassonPoiseuille, PowerLawModel, RheologicalModel,
};
pub use stokes::StokesFlow;
pub use taylor_green::TaylorGreenVortex;
pub use utils::{AnalyticalUtils, FlowGeometry};
pub use womersley::WomersleyFlow;

use nalgebra::{RealField, Vector3};

/// Trait for analytical solutions
pub trait AnalyticalSolution<T: RealField + Copy> {
    /// Evaluate the solution at given coordinates and time
    fn evaluate(&self, x: T, y: T, z: T, t: T) -> Vector3<T>;

    /// Get the velocity field at given coordinates and time
    fn velocity(&self, x: T, y: T, z: T, t: T) -> Vector3<T> {
        self.evaluate(x, y, z, t)
    }

    /// Get the pressure field at given coordinates and time
    fn pressure(&self, x: T, y: T, z: T, t: T) -> T;

    /// Get the name of the analytical solution
    fn name(&self) -> &str;

    /// Get the domain bounds [`x_min`, `x_max`, `y_min`, `y_max`, `z_min`, `z_max`]
    fn domain_bounds(&self) -> [T; 6];

    /// Check if the solution is valid at given coordinates
    fn is_valid_point(&self, x: T, y: T, z: T) -> bool {
        let bounds = self.domain_bounds();
        x >= bounds[0]
            && x <= bounds[1]
            && y >= bounds[2]
            && y <= bounds[3]
            && z >= bounds[4]
            && z <= bounds[5]
    }

    /// Get the characteristic length scale
    fn length_scale(&self) -> T;

    /// Get the characteristic velocity scale
    fn velocity_scale(&self) -> T;

    /// Get the characteristic time scale
    fn time_scale(&self) -> T {
        self.length_scale() / self.velocity_scale()
    }
}
