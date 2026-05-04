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

use super::pressure_balance::{bisect_monotone_target, ScalarSolveTolerances};
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
    ///
    /// `solve()` computes a pressure-compatible split. Use
    /// `solve_with_prescribed_split()` when these stored ratios are the desired
    /// physical split from an external model or control law.
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

    fn validate_split_ratios(&self) -> Result<(), Error> {
        let (r1, r2, r3) = self.flow_split_ratios;
        if r1 < T::zero() || r2 < T::zero() || r3 < T::zero() {
            return Err(Error::InvalidConfiguration(
                "three-way flow split ratios must be nonnegative".to_string(),
            ));
        }

        let ratio_sum = r1 + r2 + r3;
        let ratio_error = (ratio_sum - T::one()).abs();
        if ratio_error > T::from_f64_or_one(1e-12) {
            return Err(Error::InvalidConfiguration(format!(
                "three-way flow split ratios must sum to 1.0, got {}",
                ratio_sum.to_f64().unwrap_or(f64::NAN)
            )));
        }

        Ok(())
    }

    fn daughter_pressure_drop<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        daughter_index: usize,
        flow_rate: T,
        temperature: T,
        pressure: T,
    ) -> T {
        match daughter_index {
            0 => TwoWayBranchJunction::pressure_drop(
                fluid,
                flow_rate,
                &self.daughter1,
                temperature,
                pressure,
            ),
            1 => TwoWayBranchJunction::pressure_drop(
                fluid,
                flow_rate,
                &self.daughter2,
                temperature,
                pressure,
            ),
            2 => TwoWayBranchJunction::pressure_drop(
                fluid,
                flow_rate,
                &self.daughter3,
                temperature,
                pressure,
            ),
            _ => T::zero(),
        }
    }

    fn flow_for_common_pressure_drop<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        daughter_index: usize,
        target_dp: T,
        q_parent: T,
        temperature: T,
        pressure: T,
    ) -> T {
        if target_dp <= T::zero() || q_parent <= T::zero() {
            return T::zero();
        }

        let lower = T::zero();
        let upper = q_parent;
        let tolerances = ScalarSolveTolerances::for_target(target_dp, q_parent);
        bisect_monotone_target(lower, upper, target_dp, tolerances, |flow_rate| {
            self.daughter_pressure_drop(fluid, daughter_index, flow_rate, temperature, pressure)
        })
    }

    fn total_flow_for_common_pressure_drop<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        target_dp: T,
        q_parent: T,
        temperature: T,
        pressure: T,
    ) -> T {
        self.flow_for_common_pressure_drop(fluid, 0, target_dp, q_parent, temperature, pressure)
            + self.flow_for_common_pressure_drop(
                fluid,
                1,
                target_dp,
                q_parent,
                temperature,
                pressure,
            )
            + self.flow_for_common_pressure_drop(
                fluid,
                2,
                target_dp,
                q_parent,
                temperature,
                pressure,
            )
    }

    fn solve_pressure_balanced_flows<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        q_parent: T,
        temperature: T,
        pressure: T,
    ) -> Result<(T, T, T), Error> {
        let tiny_flow = T::from_f64_or_one(1e-18);
        let q_parent_magnitude = q_parent.abs();
        if q_parent_magnitude <= tiny_flow {
            return Ok((T::zero(), T::zero(), T::zero()));
        }
        let orientation = if q_parent >= T::zero() {
            T::one()
        } else {
            -T::one()
        };

        let lower_dp = T::zero();
        let mut upper_dp = self
            .daughter_pressure_drop(fluid, 0, q_parent_magnitude, temperature, pressure)
            .max(self.daughter_pressure_drop(fluid, 1, q_parent_magnitude, temperature, pressure))
            .max(self.daughter_pressure_drop(fluid, 2, q_parent_magnitude, temperature, pressure));

        let target_sum = q_parent_magnitude;
        let mut upper_sum = self.total_flow_for_common_pressure_drop(
            fluid,
            upper_dp,
            q_parent_magnitude,
            temperature,
            pressure,
        );
        for _ in 0..24 {
            if upper_sum >= target_sum {
                break;
            }
            upper_dp *= T::from_f64_or_one(2.0);
            upper_sum = self.total_flow_for_common_pressure_drop(
                fluid,
                upper_dp,
                q_parent_magnitude,
                temperature,
                pressure,
            );
        }

        if upper_sum < target_sum {
            return Err(Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::StagnatedResidual {
                    residual: (target_sum - upper_sum).to_f64().unwrap_or(f64::NAN),
                },
            ));
        }

        let tolerances = ScalarSolveTolerances::for_target(q_parent_magnitude, upper_dp);
        let final_dp = bisect_monotone_target(
            lower_dp,
            upper_dp,
            q_parent_magnitude,
            tolerances,
            |target_dp| {
                self.total_flow_for_common_pressure_drop(
                    fluid,
                    target_dp,
                    q_parent_magnitude,
                    temperature,
                    pressure,
                )
            },
        );

        let q_1 = self.flow_for_common_pressure_drop(
            fluid,
            0,
            final_dp,
            q_parent_magnitude,
            temperature,
            pressure,
        );
        let q_2 = self.flow_for_common_pressure_drop(
            fluid,
            1,
            final_dp,
            q_parent_magnitude,
            temperature,
            pressure,
        );
        let q_3 = (q_parent_magnitude - q_1 - q_2).max(T::zero());
        Ok((orientation * q_1, orientation * q_2, orientation * q_3))
    }

    fn solve_from_flows<F: FluidTrait<T> + NonNewtonianFluid<T> + Copy>(
        &self,
        fluid: F,
        q_parent: T,
        p_parent: T,
        q_1: T,
        q_2: T,
        q_3: T,
        temperature: T,
        pressure: T,
    ) -> Result<ThreeWayBranchSolution<T>, Error> {
        let q_sum = q_1 + q_2 + q_3;
        let mass_error = (q_sum - q_parent).abs() / q_parent.abs().max(T::from_f64_or_one(1e-15));

        if mass_error > T::from_f64_or_one(1e-10) {
            use cfd_core::error::ConvergenceErrorKind;
            return Err(Error::Convergence(ConvergenceErrorKind::Diverged {
                norm: mass_error.to_f64().unwrap_or(f64::NAN),
            }));
        }

        let dp_parent = TwoWayBranchJunction::pressure_drop(
            &fluid,
            q_parent,
            &self.parent,
            temperature,
            pressure,
        );
        let dp_1 = TwoWayBranchJunction::pressure_drop(
            &fluid,
            q_1,
            &self.daughter1,
            temperature,
            pressure,
        );
        let dp_2 = TwoWayBranchJunction::pressure_drop(
            &fluid,
            q_2,
            &self.daughter2,
            temperature,
            pressure,
        );
        let dp_3 = TwoWayBranchJunction::pressure_drop(
            &fluid,
            q_3,
            &self.daughter3,
            temperature,
            pressure,
        );

        let p_junction = p_parent - dp_parent;
        let p_1 = p_junction - dp_1;
        let p_2 = p_junction - dp_2;
        let p_3 = p_junction - dp_3;

        let p_max = p_1.max(p_2).max(p_3);
        let p_min = p_1.min(p_2).min(p_3);
        let junction_pressure_error =
            (p_max - p_min).abs() / (p_junction.abs() + T::from_f64_or_one(1.0));

        let gamma_1 = TwoWayBranchJunction::shear_rate(q_1, &self.daughter1);
        let gamma_2 = TwoWayBranchJunction::shear_rate(q_2, &self.daughter2);
        let gamma_3 = TwoWayBranchJunction::shear_rate(q_3, &self.daughter3);

        let mu_1 = TwoWayBranchJunction::apparent_viscosity(
            &fluid,
            q_1,
            &self.daughter1,
            temperature,
            pressure,
        );
        let mu_2 = TwoWayBranchJunction::apparent_viscosity(
            &fluid,
            q_2,
            &self.daughter2,
            temperature,
            pressure,
        );
        let mu_3 = TwoWayBranchJunction::apparent_viscosity(
            &fluid,
            q_3,
            &self.daughter3,
            temperature,
            pressure,
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

    /// Solve a three-way branch with given parent pressure, flow rate, and thermodynamic state.
    pub fn solve<F: FluidTrait<T> + NonNewtonianFluid<T> + Copy>(
        &self,
        fluid: F,
        q_parent: T,
        p_parent: T,
        temperature: T,
        pressure: T,
    ) -> Result<ThreeWayBranchSolution<T>, Error> {
        let (q_1, q_2, q_3) =
            self.solve_pressure_balanced_flows(&fluid, q_parent, temperature, pressure)?;

        self.solve_from_flows(
            fluid,
            q_parent,
            p_parent,
            q_1,
            q_2,
            q_3,
            temperature,
            pressure,
        )
    }

    /// Solve a three-way branch using the stored daughter split ratios exactly.
    pub fn solve_with_prescribed_split<F: FluidTrait<T> + NonNewtonianFluid<T> + Copy>(
        &self,
        fluid: F,
        q_parent: T,
        p_parent: T,
        temperature: T,
        pressure: T,
    ) -> Result<ThreeWayBranchSolution<T>, Error> {
        self.validate_split_ratios()?;

        let q_1 = self.flow_split_ratios.0 * q_parent;
        let q_2 = self.flow_split_ratios.1 * q_parent;
        let q_3 = self.flow_split_ratios.2 * q_parent;
        self.solve_from_flows(
            fluid,
            q_parent,
            p_parent,
            q_1,
            q_2,
            q_3,
            temperature,
            pressure,
        )
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
