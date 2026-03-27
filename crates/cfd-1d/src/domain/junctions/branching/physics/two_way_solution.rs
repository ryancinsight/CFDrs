//! Solution type for two-way branch junction problems.

use cfd_core::conversion::SafeFromF64;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};
use std::fmt;

/// Solution to the two-way branch problem
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct TwoWayBranchSolution<T: RealField + Copy> {
    /// Parent branch volumetric flow rate [m³/s]
    pub q_parent: T,
    /// Daughter 1 volumetric flow rate [m³/s]
    pub q_1: T,
    /// Daughter 2 volumetric flow rate [m³/s]
    pub q_2: T,

    /// Parent inlet pressure [Pa]
    pub p_parent: T,
    /// Junction pressure [Pa]
    pub p_junction: T,
    /// Daughter 1 outlet pressure [Pa]
    pub p_1: T,
    /// Daughter 2 outlet pressure [Pa]
    pub p_2: T,

    /// Pressure drop in parent [Pa]
    pub dp_parent: T,
    /// Pressure drop in daughter 1 [Pa]
    pub dp_1: T,
    /// Pressure drop in daughter 2 [Pa]
    pub dp_2: T,

    /// Wall shear rate in daughter 1 [1/s]
    pub gamma_1: T,
    /// Wall shear rate in daughter 2 [1/s]
    pub gamma_2: T,

    /// Apparent viscosity in daughter 1 [Pa·s]
    pub mu_1: T,
    /// Apparent viscosity in daughter 2 [Pa·s]
    pub mu_2: T,

    /// Pressure continuity error at junction |P_1 - P_2| / P_parent
    pub junction_pressure_error: T,
    /// Mass conservation error |Q_1 + Q_2 - Q_0| / Q_0
    pub mass_conservation_error: T,
}

impl<T: RealField + Copy> TwoWayBranchSolution<T> {
    /// Check if solution satisfies conservation laws within tolerance
    pub fn is_valid(&self, tolerance: T) -> bool {
        self.mass_conservation_error < tolerance && self.junction_pressure_error < tolerance
    }

    /// Get the flow split ratio Q_1 / Q_parent
    pub fn flow_ratio(&self) -> T {
        self.q_1 / (self.q_parent + T::from_f64_or_one(1e-15))
    }
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> fmt::Display
    for TwoWayBranchSolution<T>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let q_parent_f = self.q_parent.to_f64().unwrap_or(f64::NAN);
        let q_1_f = self.q_1.to_f64().unwrap_or(f64::NAN);
        let q_2_f = self.q_2.to_f64().unwrap_or(f64::NAN);
        let p_parent_f = self.p_parent.to_f64().unwrap_or(f64::NAN);
        let p_1_f = self.p_1.to_f64().unwrap_or(f64::NAN);
        let p_2_f = self.p_2.to_f64().unwrap_or(f64::NAN);
        let junction_error_f = self.junction_pressure_error.to_f64().unwrap_or(f64::NAN);
        let mass_error_f = self.mass_conservation_error.to_f64().unwrap_or(f64::NAN);

        write!(
            f,
            "TwoWayBranchSolution {{\n  Parent: Q={q_parent_f:.2e} m³/s, P={p_parent_f} Pa\n  \
             Daughter1: Q={q_1_f:.2e} m³/s, P={p_1_f} Pa\n  \
             Daughter2: Q={q_2_f:.2e} m³/s, P={p_2_f} Pa\n  \
             Junction P error: {junction_error_f:.2e}, Mass error: {mass_error_f:.2e}\n}}"
        )
    }
}
