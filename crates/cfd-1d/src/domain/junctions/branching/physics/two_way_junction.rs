//! Two-way branch junction model with full conservation equations.
//!
//! Implements the junction flow model governing pressure and flow distribution
//! at a bifurcation point.

use crate::domain::channel::{Channel, CrossSection};
use cfd_core::conversion::SafeFromF64;
use cfd_core::error::Error;
use cfd_core::physics::fluid::traits::Fluid as FluidTrait;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};

use super::pressure_balance::{bisect_root, ScalarSolveTolerances};
use super::two_way_solution::TwoWayBranchSolution;

/// Two-way branch junction connecting parent channel to two daughter channels
///
/// # Physics
///
/// A two-way branch junction conserves mass and imposes pressure continuity:
///
/// **Mass conservation:**
/// ```text
/// Q_0 = Q_1 + Q_2
/// ```
///
/// **Pressure relationship:**
/// ```text
/// ΔP_i = P_0 - P_i = (128 μ_i Q_i L_i) / (π D_i^4)  [Poiseuille]
/// ```
///
/// For non-Newtonian fluids, viscosity depends on shear rate in each branch:
/// ```text
/// γ̇_i = (32 Q_i) / (π D_i^3)
/// ```
///
/// # References
///
/// - Huo & Kassab (2012): Validated branching diameter and area ratios
/// - Fung (1993): Branch-junction pressure losses in biological networks
/// - Murray's Law: D_0^3 = D_1^3 + D_2^3
#[derive(Debug, Clone)]
pub struct TwoWayBranchJunction<T: RealField + Copy> {
    /// Parent channel (incoming flow)
    pub parent: Channel<T>,
    /// First daughter channel (outgoing)
    pub daughter1: Channel<T>,
    /// Second daughter channel (outgoing)
    pub daughter2: Channel<T>,
    /// Flow-split seed or prescribed ratio: Q_1 / (Q_1 + Q_2)
    ///
    /// `solve()` uses this value as the initial bracket refinement seed for the
    /// pressure-balanced bifurcation solve. Use `solve_with_prescribed_split()`
    /// to force this exact ratio.
    pub flow_split_ratio: T,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> TwoWayBranchJunction<T> {
    /// Create a new two-way branch junction
    pub fn new(
        parent: Channel<T>,
        daughter1: Channel<T>,
        daughter2: Channel<T>,
        flow_split_ratio: T,
    ) -> Self {
        Self {
            parent,
            daughter1,
            daughter2,
            flow_split_ratio,
        }
    }

    /// Get hydraulic diameter from a channel's cross section
    pub(crate) fn hydraulic_diameter(channel: &Channel<T>) -> T {
        match &channel.geometry.cross_section {
            CrossSection::Circular { diameter } => *diameter,
            CrossSection::Rectangular { width, height } => {
                let four = T::from_f64_or_one(4.0);
                let two = T::from_f64_or_one(2.0);
                let area = *width * *height;
                let perimeter = two * (*width + *height);
                four * area / perimeter
            }
            CrossSection::Elliptical {
                major_axis,
                minor_axis,
            } => {
                let pi = T::from_f64_or_one(std::f64::consts::PI);
                let two = T::from_f64_or_one(2.0);
                let four = T::from_f64_or_one(4.0);

                let a_val = *major_axis / two;
                let b_val = *minor_axis / two;

                let (a, b) = if a_val > b_val {
                    (a_val, b_val)
                } else {
                    (b_val, a_val)
                };

                if a == b || b == T::from_f64_or_one(0.0) {
                    return two * a;
                }

                let m = T::from_f64_or_one(1.0) - (b * b) / (a * a);

                let mut a_n = T::from_f64_or_one(1.0);
                let mut b_n = (T::from_f64_or_one(1.0) - m).sqrt();
                let mut c_n = m.sqrt();

                let mut sum = c_n * c_n / two;
                let mut power = T::from_f64_or_one(1.0);
                let tolerance = T::from_f64_or_one(1e-14);

                for _ in 0..20 {
                    let a_next = (a_n + b_n) / two;
                    let b_next = (a_n * b_n).sqrt();
                    let c_next = (a_n - b_n) / two;

                    a_n = a_next;
                    b_n = b_next;
                    c_n = c_next;

                    sum += power * c_n * c_n;
                    power *= two;

                    if c_n < tolerance || c_n == T::from_f64_or_one(0.0) {
                        break;
                    }
                }

                let e_m = (pi / (two * a_n)) * (T::from_f64_or_one(1.0) - sum);
                let perimeter = four * a * e_m;
                let area = pi * a * b;
                four * area / perimeter
            }
            CrossSection::Trapezoidal {
                top_width,
                bottom_width,
                height,
            } => {
                let two = T::from_f64_or_one(2.0);
                let four = T::from_f64_or_one(4.0);
                let area = (*top_width + *bottom_width) * *height / two;
                let side_length =
                    ((*top_width - *bottom_width).powi(2) / four + height.powi(2)).sqrt();
                let perimeter = *top_width + *bottom_width + two * side_length;
                four * area / perimeter
            }
            CrossSection::Custom {
                hydraulic_diameter, ..
            } => *hydraulic_diameter,
        }
    }

    /// Validate the two-way branch satisfies Murray's law (D_0^3 = D_1^3 + D_2^3)
    ///
    /// Returns deviation from Murray's law as a fraction:
    /// ```text
    /// deviation = |D_0^3 - (D_1^3 + D_2^3)| / D_0^3
    /// ```
    pub fn murray_law_deviation(&self) -> T {
        let three = T::from_f64_or_one(3.0);
        let d0 = Self::hydraulic_diameter(&self.parent);
        let d1 = Self::hydraulic_diameter(&self.daughter1);
        let d2 = Self::hydraulic_diameter(&self.daughter2);
        let d0_cubed = d0.powf(three);
        let d1_cubed = d1.powf(three);
        let d2_cubed = d2.powf(three);
        (d0_cubed - (d1_cubed + d2_cubed)).abs() / d0_cubed.max(T::from_f64_or_one(1e-10))
    }

    /// Calculate shear rate in a channel for given volumetric flow rate
    ///
    /// ```text
    /// γ̇_wall = (32 Q) / (π D^3)
    /// ```
    pub(crate) fn shear_rate(q: T, channel: &Channel<T>) -> T {
        let d = Self::hydraulic_diameter(channel);
        let pi = T::from_f64_or_one(std::f64::consts::PI);
        let thirty_two = T::from_f64_or_one(32.0);
        (thirty_two * q) / (pi * d * d * d)
    }

    /// Calculate apparent viscosity in a channel for given flow rate.
    ///
    /// For non-Newtonian fluids, viscosity depends on wall shear rate via
    /// the constitutive model (Casson, Carreau-Yasuda, etc.).
    pub(crate) fn apparent_viscosity<F: FluidTrait<T>>(
        fluid: &F,
        q: T,
        channel: &Channel<T>,
        temperature: T,
        pressure: T,
    ) -> T {
        let gamma = Self::shear_rate(q, channel);
        fluid
            .viscosity_at_shear(gamma, temperature, pressure)
            .unwrap_or_else(|_| {
                fluid.properties_at(temperature, pressure).map_or_else(
                    |_| T::from_f64_or_one(0.0035),
                    |state| state.dynamic_viscosity,
                )
            })
    }

    /// Calculate pressure drop across a channel using Hagen-Poiseuille equation.
    ///
    /// ```text
    /// ΔP = (128 μ Q L) / (π D^4)
    /// ```
    pub fn pressure_drop<F: cfd_core::physics::fluid::FluidTrait<T>>(
        fluid: &F,
        q: T,
        channel: &Channel<T>,
        temperature: T,
        pressure: T,
    ) -> T {
        let one_two_eight = T::from_f64_or_one(128.0);
        let pi = T::from_f64_or_one(std::f64::consts::PI);
        let mu = Self::apparent_viscosity(fluid, q, channel, temperature, pressure);
        let d = Self::hydraulic_diameter(channel);
        let l = channel.geometry.length;
        (one_two_eight * mu * q * l) / (pi * d * d * d * d)
    }

    fn validate_split_ratio(&self) -> Result<(), Error> {
        if self.flow_split_ratio < T::zero() || self.flow_split_ratio > T::one() {
            return Err(Error::InvalidConfiguration(format!(
                "flow split ratio must be within [0, 1], got {}",
                self.flow_split_ratio.to_f64().unwrap_or(f64::NAN)
            )));
        }
        Ok(())
    }

    fn daughter_pressure_residual<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        q_1: T,
        q_parent: T,
        temperature: T,
        pressure: T,
    ) -> T {
        let q_2 = q_parent - q_1;
        let dp_1 = Self::pressure_drop(fluid, q_1, &self.daughter1, temperature, pressure);
        let dp_2 = Self::pressure_drop(fluid, q_2, &self.daughter2, temperature, pressure);
        dp_1 - dp_2
    }

    fn solve_pressure_balanced_flow<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        q_parent: T,
        temperature: T,
        pressure: T,
    ) -> Result<T, Error> {
        self.validate_split_ratio()?;

        if q_parent < T::zero() {
            return Err(Error::InvalidInput(
                "two-way branch solve requires nonnegative parent flow".to_string(),
            ));
        }

        let tiny_flow = T::from_f64_or_one(1e-18);
        if q_parent.abs() <= tiny_flow {
            return Ok(T::zero());
        }

        let tolerances = ScalarSolveTolerances::for_flow_interval(q_parent);

        let lower_residual =
            self.daughter_pressure_residual(fluid, T::zero(), q_parent, temperature, pressure);
        let upper_residual =
            self.daughter_pressure_residual(fluid, q_parent, q_parent, temperature, pressure);

        if lower_residual > T::zero() || upper_residual < T::zero() {
            return Err(Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::StagnatedResidual {
                    residual: (lower_residual.abs() + upper_residual.abs())
                        .to_f64()
                        .unwrap_or(f64::NAN),
                },
            ));
        }

        Ok(bisect_root(
            T::zero(),
            q_parent,
            Some(self.flow_split_ratio * q_parent),
            tolerances,
            |q_1| self.daughter_pressure_residual(fluid, q_1, q_parent, temperature, pressure),
        ))
    }

    fn solve_from_split<F: FluidTrait<T>>(
        &self,
        fluid: F,
        q_parent: T,
        p_parent: T,
        q_1: T,
        temperature: T,
        pressure: T,
    ) -> Result<TwoWayBranchSolution<T>, Error> {
        let q_2 = q_parent - q_1;

        let q_sum = q_1 + q_2;
        let mass_error = (q_sum - q_parent).abs() / q_parent.max(T::from_f64_or_one(1e-15));
        if mass_error > T::from_f64_or_one(1e-10) {
            use cfd_core::error::ConvergenceErrorKind;
            return Err(Error::Convergence(ConvergenceErrorKind::Diverged {
                norm: mass_error.to_f64().unwrap_or(f64::NAN),
            }));
        }

        let dp_parent = Self::pressure_drop(&fluid, q_parent, &self.parent, temperature, pressure);
        let dp_1 = Self::pressure_drop(&fluid, q_1, &self.daughter1, temperature, pressure);
        let dp_2 = Self::pressure_drop(&fluid, q_2, &self.daughter2, temperature, pressure);

        let p_junction = p_parent - dp_parent;
        let p_1 = p_junction - dp_1;
        let p_2 = p_junction - dp_2;

        let junction_pressure_error =
            (p_1 - p_2).abs() / (p_junction.abs() + T::from_f64_or_one(1.0));

        let gamma_1 = Self::shear_rate(q_1, &self.daughter1);
        let gamma_2 = Self::shear_rate(q_2, &self.daughter2);
        let mu_1 = Self::apparent_viscosity(&fluid, q_1, &self.daughter1, temperature, pressure);
        let mu_2 = Self::apparent_viscosity(&fluid, q_2, &self.daughter2, temperature, pressure);

        Ok(TwoWayBranchSolution {
            q_parent,
            q_1,
            q_2,
            p_parent,
            p_junction,
            p_1,
            p_2,
            dp_parent,
            dp_1,
            dp_2,
            gamma_1,
            gamma_2,
            mu_1,
            mu_2,
            junction_pressure_error,
            mass_conservation_error: mass_error,
        })
    }

    /// Solve a two-way branch by enforcing daughter pressure compatibility.
    ///
    /// The solver finds the unique split `Q_1` in `[0, Q_parent]` such that the
    /// daughter pressure drops match: `ΔP_1(Q_1) = ΔP_2(Q_parent - Q_1)`.
    pub fn solve<F: FluidTrait<T> + Copy>(
        &self,
        fluid: F,
        q_parent: T,
        p_parent: T,
        temperature: T,
        pressure: T,
    ) -> Result<TwoWayBranchSolution<T>, Error> {
        let q_1 = self.solve_pressure_balanced_flow(&fluid, q_parent, temperature, pressure)?;
        self.solve_from_split(fluid, q_parent, p_parent, q_1, temperature, pressure)
    }

    /// Solve a two-way branch using the stored split ratio exactly.
    ///
    /// This is intended for callers that already own an external phase-separation
    /// or control-law split and want the resulting pressure diagnostics without
    /// re-solving the bifurcation compatibility condition.
    pub fn solve_with_prescribed_split<F: FluidTrait<T> + Copy>(
        &self,
        fluid: F,
        q_parent: T,
        p_parent: T,
        temperature: T,
        pressure: T,
    ) -> Result<TwoWayBranchSolution<T>, Error> {
        self.validate_split_ratio()?;
        if q_parent < T::zero() {
            return Err(Error::InvalidInput(
                "two-way branch solve requires nonnegative parent flow".to_string(),
            ));
        }
        let q_1 = self.flow_split_ratio * q_parent;
        self.solve_from_split(fluid, q_parent, p_parent, q_1, temperature, pressure)
    }
}
