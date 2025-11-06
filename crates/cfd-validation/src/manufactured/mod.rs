//! Manufactured solutions for verification of numerical methods
//!
//! This module provides manufactured solutions with known analytical forms
//! that can be used to verify the correctness of numerical implementations.
//!
//! References:
//! - Roache, P.J. (2002) "Code Verification by the Method of Manufactured Solutions"
//! - Salari, K. & Knupp, P. (2000) "Code Verification by the Method of Manufactured Solutions"

use nalgebra::{RealField, ComplexField};

pub mod advanced_physics;
pub mod advection;
pub mod advection_diffusion;
pub mod burgers;
pub mod diffusion;
pub mod multi_physics;
pub mod navier_stokes;
pub mod reynolds_stress_mms;
pub mod richardson_integration;
pub mod turbulent;

pub use advanced_physics::{ManufacturedCompressibleEuler, ManufacturedHypersonic, ManufacturedShockCapturing, ManufacturedTCI};
pub use advection::ManufacturedAdvection;
pub use advection_diffusion::ManufacturedAdvectionDiffusion;
pub use burgers::ManufacturedBurgers;
pub use diffusion::ManufacturedDiffusion;
pub use multi_physics::{ManufacturedConjugateHeatTransfer, ManufacturedMHD, ManufacturedMultiphase, ManufacturedSpeciesTransport};
pub use navier_stokes::{TaylorGreenManufactured, NavierStokesManufacturedSolution, PolynomialNavierStokesMMS};
pub use reynolds_stress_mms::{ManufacturedReynoldsStressMMS, PressureStrainModelMMS, ReynoldsStressConvergenceStudy};
pub use richardson_integration::{MmsRichardsonStudy, RichardsonMmsResult};
pub use turbulent::{ManufacturedKEpsilon, ManufacturedKOmega, ManufacturedReynoldsStress, ManufacturedSpalartAllmaras};

/// Trait for manufactured solutions
pub trait ManufacturedSolution<T: RealField + Copy> {
    /// Evaluate the exact solution at given point and time
    fn exact_solution(&self, x: T, y: T, z: T, t: T) -> T;

    /// Evaluate the source term required to satisfy the governing equation
    fn source_term(&self, x: T, y: T, z: T, t: T) -> T;

    /// Get the initial condition at a given point
    fn initial_condition(&self, x: T, y: T, z: T) -> T {
        self.exact_solution(x, y, z, T::zero())
    }

    /// Get the boundary condition type and value
    fn boundary_condition(&self, x: T, y: T, z: T, t: T) -> T {
        self.exact_solution(x, y, z, t)
    }

    /// Calculate the L2 error norm between numerical and exact solutions
    fn calculate_error_norm(&self, numerical: &[T], exact: &[T]) -> T {
        let mut sum = T::zero();
        let n = T::from_f64(numerical.len() as f64).unwrap();

        for (num, ex) in numerical.iter().zip(exact.iter()) {
            let diff = *num - *ex;
            sum += diff * diff;
        }

        ComplexField::sqrt(sum / n)
    }

    /// Calculate the maximum error between numerical and exact solutions
    fn calculate_max_error(&self, numerical: &[T], exact: &[T]) -> T {
        numerical
            .iter()
            .zip(exact.iter())
            .map(|(num, ex)| ComplexField::abs(*num - *ex))
            .fold(T::zero(), |max, val| if val > max { val } else { max })
}
}

/// Common manufactured solution functions
pub struct ManufacturedFunctions;

impl ManufacturedFunctions {
    /// Sinusoidal solution: sin(kx * x) * sin(ky * y) * exp(-t)
    pub fn sinusoidal<T: RealField + Float>(x: T, y: T, t: T, kx: T, ky: T) -> T {
        Float::sin(kx * x) * Float::sin(ky * y) * Float::exp(-t)
    }

    /// Polynomial solution: x^2 + y^2 + t
    pub fn polynomial<T: RealField + Float>(x: T, y: T, t: T) -> T {
        x * x + y * y + t
    }

    /// Exponential solution: exp(x + y - t)
    pub fn exponential<T: RealField + Float>(x: T, y: T, t: T) -> T {
        Float::exp(x + y - t)
    }

    /// Trigonometric-exponential: cos(x) * sin(y) * exp(-2t)
    pub fn trig_exp<T: RealField + Float>(x: T, y: T, t: T) -> T {
        Float::cos(x) * Float::sin(y) * Float::exp(T::from(-2.0).unwrap() * t)
    }
}
