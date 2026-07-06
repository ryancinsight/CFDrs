//! Solution type for two-way branch junction problems.

use crate::scalar::Cfd1dScalar;
use cfd_core::conversion::SafeFromF64;
use eunomia::NumericElement;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Solution to the two-way branch problem
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct TwoWayBranchSolution<T: Cfd1dScalar + Copy> {
    /// Parent branch volumetric flow rate [m³/s]
    pub q_parent: T,
    /// Daughter 1 volumetric flow rate [m³/s]
    pub q_1: T,
    /// Daughter 2 volumetric flow rate [m³/s]
    pub q_2: T,

    /// Parent inlet pressure \[Pa]
    pub p_parent: T,
    /// Junction pressure \[Pa]
    pub p_junction: T,
    /// Daughter 1 outlet pressure \[Pa]
    pub p_1: T,
    /// Daughter 2 outlet pressure \[Pa]
    pub p_2: T,

    /// Pressure drop in parent \[Pa]
    pub dp_parent: T,
    /// Pressure drop in daughter 1 \[Pa]
    pub dp_1: T,
    /// Pressure drop in daughter 2 \[Pa]
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

impl<T: Cfd1dScalar + Copy + SafeFromF64> TwoWayBranchSolution<T> {
    /// Check if solution satisfies conservation laws within tolerance
    pub fn is_valid(&self, tolerance: T) -> bool {
        self.mass_conservation_error < tolerance && self.junction_pressure_error < tolerance
    }

    /// Get the flow split ratio Q_1 / Q_parent
    pub fn flow_ratio(&self) -> T {
        self.q_1 / (self.q_parent + T::from_f64_or_one(1e-15))
    }
}

impl<T: Cfd1dScalar + Copy + SafeFromF64> fmt::Display for TwoWayBranchSolution<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let q_parent_f = <T as NumericElement>::to_f64(self.q_parent);
        let q_1_f = <T as NumericElement>::to_f64(self.q_1);
        let q_2_f = <T as NumericElement>::to_f64(self.q_2);
        let p_parent_f = <T as NumericElement>::to_f64(self.p_parent);
        let p_1_f = <T as NumericElement>::to_f64(self.p_1);
        let p_2_f = <T as NumericElement>::to_f64(self.p_2);
        let junction_error_f = <T as NumericElement>::to_f64(self.junction_pressure_error);
        let mass_error_f = <T as NumericElement>::to_f64(self.mass_conservation_error);

        write!(
            f,
            "TwoWayBranchSolution {{\n  Parent: Q={q_parent_f:.2e} m³/s, P={p_parent_f} Pa\n  \
             Daughter1: Q={q_1_f:.2e} m³/s, P={p_1_f} Pa\n  \
             Daughter2: Q={q_2_f:.2e} m³/s, P={p_2_f} Pa\n  \
             Junction P error: {junction_error_f:.2e}, Mass error: {mass_error_f:.2e}\n}}"
        )
    }
}
