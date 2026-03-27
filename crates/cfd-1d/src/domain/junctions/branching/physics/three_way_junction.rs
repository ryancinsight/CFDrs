//! Three-way branch junction model with conservation equations.
//!
//! Generalizes the two-way branch formulation to a three-way split
//! (one parent to three daughters).

use crate::domain::channel::Channel;
use cfd_core::conversion::SafeFromF64;
use cfd_core::error::Error;
use cfd_core::physics::fluid::traits::{Fluid as FluidTrait, NonNewtonianFluid};
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

use super::two_way_junction::TwoWayBranchJunction;

/// Three-way branch junction (one parent to three daughters)
///
/// Generalizes the two-way branch formulation to a three-way split.
#[derive(Debug, Clone)]
pub struct ThreeWayBranchJunction<T: RealField + Copy> {
    /// Parent channel (incoming flow)
    pub parent: Channel<T>,
    /// First daughter channel
    pub daughter1: Channel<T>,
    /// Second daughter channel
    pub daughter2: Channel<T>,
    /// Third daughter channel
    pub daughter3: Channel<T>,
    /// Flow distribution ratios: (Q_1/Q_0, Q_2/Q_0, Q_3/Q_0)
    /// Sum should equal 1.0
    pub flow_split_ratios: (T, T, T),
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> ThreeWayBranchJunction<T> {
    /// Create a new three-way branch junction
    pub fn new(
        parent: Channel<T>,
        daughter1: Channel<T>,
        daughter2: Channel<T>,
        daughter3: Channel<T>,
        flow_split_ratios: (T, T, T),
    ) -> Self {
        Self {
            parent,
            daughter1,
            daughter2,
            daughter3,
            flow_split_ratios,
        }
    }

    /// Validate the three-way branch against Murray's law extension:
    ///
    /// ```text
    /// D_0^3 = D_1^3 + D_2^3 + D_3^3
    /// ```
    pub fn murray_law_deviation(&self) -> T {
        let three = T::from_f64_or_one(3.0);
        let d0 = TwoWayBranchJunction::hydraulic_diameter(&self.parent);
        let d1 = TwoWayBranchJunction::hydraulic_diameter(&self.daughter1);
        let d2 = TwoWayBranchJunction::hydraulic_diameter(&self.daughter2);
        let d3 = TwoWayBranchJunction::hydraulic_diameter(&self.daughter3);
        let d0_cubed = d0.powf(three);
        let daughters_sum = d1.powf(three) + d2.powf(three) + d3.powf(three);
        (d0_cubed - daughters_sum).abs() / d0_cubed.max(T::from_f64_or_one(1e-10))
    }

    /// Solve a three-way branch with given parent pressure, flow rate, and thermodynamic state.
    pub fn solve<F: FluidTrait<T> + NonNewtonianFluid<T> + Copy>(
        &self,
        fluid: F,
        q_parent: T,
        p_parent: T,
        temperature: T,
        pressure: T,
    ) -> Result<ThreeWayBranchSolution<T>, Error> {
        let q_1 = self.flow_split_ratios.0 * q_parent;
        let q_2 = self.flow_split_ratios.1 * q_parent;
        let q_3 = self.flow_split_ratios.2 * q_parent;

        let q_sum = q_1 + q_2 + q_3;
        let mass_error = (q_sum - q_parent).abs() / q_parent.max(T::from_f64_or_one(1e-15));

        if mass_error > T::from_f64_or_one(1e-10) {
            use cfd_core::error::ConvergenceErrorKind;
            return Err(Error::Convergence(ConvergenceErrorKind::Diverged {
                norm: mass_error.to_f64().unwrap_or(f64::NAN),
            }));
        }

        let dp_1 = TwoWayBranchJunction::pressure_drop(
            &fluid, q_1, &self.daughter1, temperature, pressure,
        );
        let dp_2 = TwoWayBranchJunction::pressure_drop(
            &fluid, q_2, &self.daughter2, temperature, pressure,
        );
        let dp_3 = TwoWayBranchJunction::pressure_drop(
            &fluid, q_3, &self.daughter3, temperature, pressure,
        );

        let p_1 = p_parent - dp_1;
        let p_2 = p_parent - dp_2;
        let p_3 = p_parent - dp_3;

        let p_max = p_1.max(p_2).max(p_3);
        let p_min = p_1.min(p_2).min(p_3);
        let junction_pressure_error =
            (p_max - p_min).abs() / (p_parent.abs() + T::from_f64_or_one(1.0));

        let gamma_1 = TwoWayBranchJunction::shear_rate(q_1, &self.daughter1);
        let gamma_2 = TwoWayBranchJunction::shear_rate(q_2, &self.daughter2);
        let gamma_3 = TwoWayBranchJunction::shear_rate(q_3, &self.daughter3);

        let mu_1 = TwoWayBranchJunction::apparent_viscosity(
            &fluid, q_1, &self.daughter1, temperature, pressure,
        );
        let mu_2 = TwoWayBranchJunction::apparent_viscosity(
            &fluid, q_2, &self.daughter2, temperature, pressure,
        );
        let mu_3 = TwoWayBranchJunction::apparent_viscosity(
            &fluid, q_3, &self.daughter3, temperature, pressure,
        );

        Ok(ThreeWayBranchSolution {
            q_parent,
            q_1,
            q_2,
            q_3,
            p_parent,
            p_1,
            p_2,
            p_3,
            dp_1,
            dp_2,
            dp_3,
            gamma_1,
            gamma_2,
            gamma_3,
            mu_1,
            mu_2,
            mu_3,
            junction_pressure_error,
            mass_conservation_error: mass_error,
        })
    }
}

/// Solution to the three-way branch problem
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ThreeWayBranchSolution<T: RealField + Copy> {
    /// Parent flow rate [m³/s]
    pub q_parent: T,
    /// Daughter-1 flow rate [m³/s]
    pub q_1: T,
    /// Daughter-2 flow rate [m³/s]
    pub q_2: T,
    /// Daughter-3 flow rate [m³/s]
    pub q_3: T,

    /// Parent pressure [Pa]
    pub p_parent: T,
    /// Daughter-1 pressure [Pa]
    pub p_1: T,
    /// Daughter-2 pressure [Pa]
    pub p_2: T,
    /// Daughter-3 pressure [Pa]
    pub p_3: T,

    /// Parent-to-daughter-1 pressure drop [Pa]
    pub dp_1: T,
    /// Parent-to-daughter-2 pressure drop [Pa]
    pub dp_2: T,
    /// Parent-to-daughter-3 pressure drop [Pa]
    pub dp_3: T,

    /// Daughter-1 wall shear rate [1/s]
    pub gamma_1: T,
    /// Daughter-2 wall shear rate [1/s]
    pub gamma_2: T,
    /// Daughter-3 wall shear rate [1/s]
    pub gamma_3: T,

    /// Daughter-1 apparent viscosity [Pa·s]
    pub mu_1: T,
    /// Daughter-2 apparent viscosity [Pa·s]
    pub mu_2: T,
    /// Daughter-3 apparent viscosity [Pa·s]
    pub mu_3: T,

    /// Junction pressure error
    pub junction_pressure_error: T,
    /// Mass conservation error
    pub mass_conservation_error: T,
}

impl<T: RealField + Copy> ThreeWayBranchSolution<T> {
    /// Check if solution is valid
    pub fn is_valid(&self, tolerance: T) -> bool {
        self.mass_conservation_error < tolerance && self.junction_pressure_error < tolerance
    }
}
