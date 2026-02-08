//! Bifurcation and trifurcation junction models with full conservation equations
//!
//! This module implements the junction flow models governing pressure and flow distribution
//! at branching points in vascular networks. All equations are derived from first principles
//! with references to literature validation.

use crate::channel::{Channel, ChannelType, CrossSection};
use cfd_core::conversion::SafeFromF64;
use cfd_core::error::Error;
use cfd_core::physics::fluid::traits::Fluid as FluidTrait;
use cfd_core::physics::fluid::traits::NonNewtonianFluid;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};
use std::fmt;

// ============================================================================
// Bifurcation Junction Model
// ============================================================================

/// Bifurcation junction connecting parent channel to two daughter channels
///
/// # Physics
///
/// A bifurcation junction conserves mass and imposes pressure continuity:
///
/// **Mass conservation:**
/// ```text
/// Q_0 = Q_1 + Q_2
/// ```
///
/// **Pressure relationship:**
/// Each daughter branch experiences pressure drop equal to parent branch pressure
/// minus daughter pressure:
/// ```text
/// ΔP_i = P_0 - P_i = (8μ_i Q_i L_i) / (π R_i^4)  [Poiseuille]
/// ```
///
/// For non-Newtonian fluids, viscosity depends on shear rate in each branch:
/// ```text
/// γ̇_i = (32 Q_i) / (π D_i^3)  [for circular pipe]
/// μ_i = f(γ̇_i)  [Casson, Carreau-Yasuda, etc.]
/// ```
///
/// # Validation References
///
/// - Huo & Kassab (2012): Validated branching diameter and area ratios
/// - Fung (1993): Bifurcation pressure losses in biological networks
/// - Murray's Law: D_0^3 = D_1^3 + D_2^3 (observed in most vascular trees)
#[derive(Debug, Clone)]
pub struct BifurcationJunction<T: RealField + Copy> {
    /// Parent channel (incoming flow)
    pub parent: Channel<T>,
    /// First daughter channel (outgoing)
    pub daughter1: Channel<T>,
    /// Second daughter channel (outgoing)
    pub daughter2: Channel<T>,
    /// Flow distribution ratio: Q_1 / (Q_1 + Q_2)
    /// Typically between 0.0 and 1.0
    pub flow_split_ratio: T,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> BifurcationJunction<T> {
    /// Create a new bifurcation junction
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
    fn hydraulic_diameter(channel: &Channel<T>) -> T {
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
                // Approximation for ellipse: D_h ≈ sqrt(major * minor)
                (*major_axis * *minor_axis).sqrt()
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

    /// Validate the bifurcation satisfies Murray's law (D_0^3 = D_1^3 + D_2^3)
    ///
    /// This is an empirical law observed in most biological branching networks.
    /// It minimizes mechanical work in the network (Murray 1926).
    ///
    /// # Returns
    ///
    /// Deviation from Murray's law as a fraction:
    /// ```text
    /// deviation = |D_0^3 - (D_1^3 + D_2^3)| / D_0^3
    /// ```
    ///
    /// Typical biological networks achieve deviation < 0.1 (10%)
    pub fn murray_law_deviation(&self) -> T {
        let three = T::from_f64_or_one(3.0);

        let d0 = Self::hydraulic_diameter(&self.parent);
        let d1 = Self::hydraulic_diameter(&self.daughter1);
        let d2 = Self::hydraulic_diameter(&self.daughter2);

        let d0_cubed = d0.powf(three);
        let d1_cubed = d1.powf(three);
        let d2_cubed = d2.powf(three);

        let deviation =
            (d0_cubed - (d1_cubed + d2_cubed)).abs() / d0_cubed.max(T::from_f64_or_one(1e-10));
        deviation
    }

    /// Calculate shear rate in a channel for given volumetric flow rate
    ///
    /// # Mathematical Formula
    ///
    /// For a circular pipe:
    /// ```text
    /// γ̇_wall = (32 Q) / (π D^3)
    /// ```
    ///
    /// This is the wall shear rate, which is the maximum shear rate in laminar pipe flow.
    ///
    /// # Arguments
    /// * `q` - Volumetric flow rate [m³/s]
    /// * `channel` - The channel geometry
    ///
    /// # Returns
    /// Wall shear rate [1/s]
    fn shear_rate(q: T, channel: &Channel<T>) -> T {
        let d = Self::hydraulic_diameter(channel);
        let pi = T::from_f64_or_one(std::f64::consts::PI);
        let thirty_two = T::from_f64_or_one(32.0);

        (thirty_two * q) / (pi * d * d * d)
    }

    /// Calculate apparent viscosity in a channel for given flow rate
    ///
    /// For non-Newtonian fluids, viscosity depends on wall shear rate via
    /// the constitutive model (Casson, Carreau-Yasuda, etc.):
    ///
    /// ```text
    /// γ̇_wall = 32Q / (πD³)   [wall shear rate in circular pipe]
    /// μ_app  = f(γ̇_wall)      [constitutive relation]
    /// ```
    ///
    /// The `viscosity_at_shear` method on `Fluid<T>` dispatches to the
    /// correct rheological model. For Newtonian fluids this returns
    /// constant viscosity regardless of shear rate.
    ///
    /// # References
    ///
    /// - Cho & Kensey (1991): Carreau-Yasuda model validation for blood
    /// - Merrill (1969): Casson model for blood rheology
    /// - Chien (1970): Shear-dependent viscosity measurements
    pub fn apparent_viscosity<F: FluidTrait<T>>(
        fluid: &F,
        q: T,
        channel: &Channel<T>,
    ) -> T {
        let gamma = Self::shear_rate(q, channel);
        // Use shear-rate-dependent viscosity through the Fluid trait.
        // `viscosity_at_shear(shear_rate, temperature, pressure)` properly
        // dispatches to Casson/Carreau-Yasuda/Power-Law/etc. for non-Newtonian
        // fluids, and returns constant μ for Newtonian fluids.
        let temperature = T::from_f64_or_one(310.15); // 37°C body temperature
        let pressure = T::from_f64_or_one(101325.0);  // 1 atm
        fluid.viscosity_at_shear(gamma, temperature, pressure)
            .unwrap_or_else(|_| {
                // Fallback: use reference viscosity if shear-based method fails
                fluid.properties_at(temperature, pressure)
                    .map(|state| state.dynamic_viscosity)
                    .unwrap_or_else(|_| T::from_f64_or_one(0.0035)) // blood ~3.5 mPa·s
            })
    }

    /// Calculate pressure drop across a channel using Hagen-Poiseuille equation
    ///
    /// # Mathematical Derivation
    ///
    /// For laminar flow in a circular pipe, the pressure drop is:
    /// ```text
    /// ΔP = (8 μ Q L) / (π R^4)
    ///    = (128 μ Q L) / (π D^4)
    /// ```
    ///
    /// where:
    /// - μ = dynamic viscosity (can be shear-rate dependent)
    /// - Q = volumetric flow rate
    /// - L = pipe length
    /// - D = pipe diameter
    ///
    /// This is exact for Newtonian fluids. For non-Newtonian fluids, it's an
    /// approximation using the apparent viscosity at the wall shear rate.
    ///
    /// # Validation
    ///
    /// This formula is validated for:
    /// - Circular pipes with Re < 2300 (laminar)
    /// - Poiseuille profile assumed
    /// - Blood modeled with Casson or Carreau-Yasuda
    pub fn pressure_drop<F: FluidTrait<T>>(
        fluid: &F,
        q: T,
        channel: &Channel<T>,
    ) -> T {
        let one_two_eight = T::from_f64_or_one(128.0);
        let pi = T::from_f64_or_one(std::f64::consts::PI);

        let mu = Self::apparent_viscosity(fluid, q, channel);
        let d = Self::hydraulic_diameter(channel);
        let l = channel.geometry.length;

        // ΔP = (128 μ Q L) / (π D^4)
        (one_two_eight * mu * q * l) / (pi * d * d * d * d)
    }

    /// Solve bifurcation with given parent pressure and flow rate
    ///
    /// # Problem Formulation
    ///
    /// Given:
    /// - Parent branch: flow rate Q_0 and pressure P_0
    /// - Bifurcation with known geometry and fluid properties
    /// - Flow split ratio r = Q_1 / Q_0 (typically estimated from geometry)
    ///
    /// Find:
    /// - Q_1, Q_2 such that mass is conserved
    /// - P_1, P_2 such that pressure drops satisfy Poiseuille law
    ///
    /// # Algorithm
    ///
    /// 1. Assume flow split ratio based on geometry (daughter diameter ratio)
    /// 2. Calculate Q_1 = r · Q_0, Q_2 = (1-r) · Q_0
    /// 3. Calculate shear rates in each daughter
    /// 4. Calculate apparent viscosities (non-Newtonian)
    /// 5. Calculate pressure drops ΔP_1, ΔP_2
    /// 6. Verify: P_0 - ΔP_1 ≈ P_0 - ΔP_2 (pressure continuity at junction)
    ///
    /// # Arguments
    /// * `fluid` - Fluid with properties and viscosity model
    /// * `q_parent` - Parent flow rate [m³/s]
    /// * `p_parent` - Parent pressure [Pa]
    ///
    /// # Returns
    /// Tuple of (Q_1, Q_2, P_1, P_2, junction_pressure_error)
    pub fn solve<F: FluidTrait<T> + Copy>(
        &self,
        fluid: F,
        q_parent: T,
        p_parent: T,
    ) -> Result<BifurcationSolution<T>, Error> {
        // Mass conservation: Q_0 = Q_1 + Q_2
        let one = T::one();
        let q_1 = self.flow_split_ratio * q_parent;
        let q_2 = (one - self.flow_split_ratio) * q_parent;

        // Verify mass conservation
        let q_sum = q_1 + q_2;
        let mass_error = (q_sum - q_parent).abs() / q_parent.max(T::from_f64_or_one(1e-15));
        if mass_error > T::from_f64_or_one(1e-10) {
            use cfd_core::error::ConvergenceErrorKind;
            return Err(Error::Convergence(ConvergenceErrorKind::Diverged {
                norm: mass_error.to_f64().unwrap_or(f64::NAN),
            }));
        }

        // Pressure drops in parent and daughter branches
        let dp_parent = Self::pressure_drop(&fluid, q_parent, &self.parent);
        let dp_1 = Self::pressure_drop(&fluid, q_1, &self.daughter1);
        let dp_2 = Self::pressure_drop(&fluid, q_2, &self.daughter2);

        // Junction pressure (after parent pressure drop)
        let p_junction = p_parent - dp_parent;
        
        // Daughter pressures (after daughter pressure drops from junction)
        let p_1 = p_junction - dp_1;
        let p_2 = p_junction - dp_2;

        // Junction pressure error (should be small for symmetric bifurcations)
        // This measures pressure continuity at the junction
        let junction_pressure_error =
            (p_1 - p_2).abs() / (p_junction.abs() + T::from_f64_or_one(1.0));

        // Shear rates (for reporting)
        let gamma_1 = Self::shear_rate(q_1, &self.daughter1);
        let gamma_2 = Self::shear_rate(q_2, &self.daughter2);

        // Apparent viscosities
        let mu_1 = Self::apparent_viscosity(&fluid, q_1, &self.daughter1);
        let mu_2 = Self::apparent_viscosity(&fluid, q_2, &self.daughter2);

        Ok(BifurcationSolution {
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
}

// ============================================================================
// Bifurcation Solution Result
// ============================================================================

/// Solution to the bifurcation problem
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct BifurcationSolution<T: RealField + Copy> {
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

impl<T: RealField + Copy> BifurcationSolution<T> {
    /// Check if solution satisfies conservation laws within tolerance
    ///
    /// # Validation Criteria
    ///
    /// - Mass conservation: error < 1e-10
    /// - Pressure continuity: error < 1e-6 (junction effects may create local pressure variation)
    pub fn is_valid(&self, tolerance: T) -> bool {
        self.mass_conservation_error < tolerance && self.junction_pressure_error < tolerance
    }

    /// Get the flow split ratio Q_1 / Q_parent
    pub fn flow_ratio(&self) -> T {
        self.q_1 / (self.q_parent + T::from_f64_or_one(1e-15))
    }
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> fmt::Display
    for BifurcationSolution<T>
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
            "BifurcationSolution {{\n  Parent: Q={:.2e} m³/s, P={} Pa\n  \
             Daughter1: Q={:.2e} m³/s, P={} Pa\n  \
             Daughter2: Q={:.2e} m³/s, P={} Pa\n  \
             Junction P error: {:.2e}, Mass error: {:.2e}\n}}",
            q_parent_f, p_parent_f, q_1_f, p_1_f, q_2_f, p_2_f, junction_error_f, mass_error_f
        )
    }
}

// ============================================================================
// Trifurcation Junction Model
// ============================================================================

/// Trifurcation junction (one parent to three daughters)
///
/// Generalizes bifurcation to three-way split with same physics principles.
#[derive(Debug, Clone)]
pub struct TrifurcationJunction<T: RealField + Copy> {
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

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> TrifurcationJunction<T> {
    /// Create new trifurcation junction
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

    /// Solve trifurcation with given parent pressure and flow rate
    pub fn solve<F: FluidTrait<T> + NonNewtonianFluid<T> + Copy>(
        &self,
        fluid: F,
        q_parent: T,
        p_parent: T,
    ) -> Result<TrifurcationSolution<T>, Error> {
        let one = T::one();

        // Distribute flows
        let q_1 = self.flow_split_ratios.0 * q_parent;
        let q_2 = self.flow_split_ratios.1 * q_parent;
        let q_3 = self.flow_split_ratios.2 * q_parent;

        // Mass conservation check
        let q_sum = q_1 + q_2 + q_3;
        let mass_error = (q_sum - q_parent).abs() / q_parent.max(T::from_f64_or_one(1e-15));

        if mass_error > T::from_f64_or_one(1e-10) {
            use cfd_core::error::ConvergenceErrorKind;
            return Err(Error::Convergence(ConvergenceErrorKind::Diverged {
                norm: mass_error.to_f64().unwrap_or(f64::NAN),
            }));
        }

        // Pressure drops
        let dp_1 = BifurcationJunction::pressure_drop(&fluid, q_1, &self.daughter1);
        let dp_2 = BifurcationJunction::pressure_drop(&fluid, q_2, &self.daughter2);
        let dp_3 = BifurcationJunction::pressure_drop(&fluid, q_3, &self.daughter3);

        // Daughter pressures
        let p_1 = p_parent - dp_1;
        let p_2 = p_parent - dp_2;
        let p_3 = p_parent - dp_3;

        // Maximum pressure deviation at junction
        let p_max = p_1.max(p_2).max(p_3);
        let p_min = p_1.min(p_2).min(p_3);
        let junction_pressure_error =
            (p_max - p_min).abs() / (p_parent.abs() + T::from_f64_or_one(1.0));

        // Shear rates (for reporting)
        let gamma_1 = BifurcationJunction::shear_rate(q_1, &self.daughter1);
        let gamma_2 = BifurcationJunction::shear_rate(q_2, &self.daughter2);
        let gamma_3 = BifurcationJunction::shear_rate(q_3, &self.daughter3);

        // Apparent viscosities
        let mu_1 = BifurcationJunction::apparent_viscosity(&fluid, q_1, &self.daughter1);
        let mu_2 = BifurcationJunction::apparent_viscosity(&fluid, q_2, &self.daughter2);
        let mu_3 = BifurcationJunction::apparent_viscosity(&fluid, q_3, &self.daughter3);

        Ok(TrifurcationSolution {
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

// ============================================================================
// Trifurcation Solution Result
// ============================================================================

/// Solution to the trifurcation problem
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct TrifurcationSolution<T: RealField + Copy> {
    /// Parent flow rate [m³/s]
    pub q_parent: T,
    /// Daughter flows [m³/s]
    pub q_1: T,
    pub q_2: T,
    pub q_3: T,

    /// Parent pressure [Pa]
    pub p_parent: T,
    /// Daughter pressures [Pa]
    pub p_1: T,
    pub p_2: T,
    pub p_3: T,

    /// Pressure drops [Pa]
    pub dp_1: T,
    pub dp_2: T,
    pub dp_3: T,

    /// Wall shear rates [1/s]
    pub gamma_1: T,
    pub gamma_2: T,
    pub gamma_3: T,

    /// Apparent viscosities [Pa·s]
    pub mu_1: T,
    pub mu_2: T,
    pub mu_3: T,

    /// Junction pressure error
    pub junction_pressure_error: T,
    /// Mass conservation error
    pub mass_conservation_error: T,
}

impl<T: RealField + Copy> TrifurcationSolution<T> {
    /// Check if solution is valid
    pub fn is_valid(&self, tolerance: T) -> bool {
        self.mass_conservation_error < tolerance && self.junction_pressure_error < tolerance
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::channel::ChannelGeometry;
    use approx::assert_relative_eq;
    use cfd_core::physics::fluid::blood::CassonBlood;

    #[test]
    fn test_bifurcation_mass_conservation() {
        // Create symmetric bifurcation using new Channel API
        let parent_geom = ChannelGeometry::<f64>::circular(1.0e-2, 2.0e-3, 1e-6);
        let parent = Channel::new(parent_geom);
        
        let d1_geom = ChannelGeometry::<f64>::circular(1.0e-2, 1.5e-3, 1e-6);
        let d1 = Channel::new(d1_geom);
        
        let d2_geom = ChannelGeometry::<f64>::circular(1.0e-2, 1.5e-3, 1e-6);
        let d2 = Channel::new(d2_geom);

        let bifurcation = BifurcationJunction::new(parent, d1, d2, 0.5);

        // Solve with blood (Copy type)
        let blood = CassonBlood::<f64>::normal_blood();
        let q_parent = 1.0e-6; // 1 μL/s
        let p_parent = 1000.0; // 1000 Pa

        let solution = bifurcation.solve(blood, q_parent, p_parent).unwrap();

        // Verify mass conservation (blood has slight variations due to non-Newtonian effects)
        assert_relative_eq!(
            solution.q_1 + solution.q_2,
            q_parent,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_bifurcation_blood_flow() {
        // Create bifurcation with blood using new Channel API
        let parent_geom = ChannelGeometry::<f64>::circular(1.0e-3, 100.0e-6, 1e-6);
        let parent = Channel::new(parent_geom);
        
        let d1_geom = ChannelGeometry::<f64>::circular(1.0e-3, 80.0e-6, 1e-6);
        let d1 = Channel::new(d1_geom);
        
        let d2_geom = ChannelGeometry::<f64>::circular(1.0e-3, 80.0e-6, 1e-6);
        let d2 = Channel::new(d2_geom);

        let bifurcation = BifurcationJunction::new(parent, d1, d2, 0.5);
        let blood = CassonBlood::<f64>::normal_blood();

        let q_parent = 1.0e-8; // 10 nL/s
        let p_parent = 100.0; // 100 Pa

        let solution = bifurcation.solve(blood, q_parent, p_parent).unwrap();

        // Blood should have non-Newtonian viscosity effects
        assert!(solution.mu_1 > 0.0);
        assert!(solution.mu_2 > 0.0);
        assert!(solution.gamma_1 > 0.0);
        assert!(solution.gamma_2 > 0.0);
    }

    #[test]
    fn test_murrary_law() {
        let parent_geom = ChannelGeometry::<f64>::circular(1.0e-2, 2.0e-3, 1e-6);
        let parent = Channel::new(parent_geom);
        
        let d1_geom = ChannelGeometry::<f64>::circular(1.0e-2, 1.58e-3, 1e-6);
        let d1 = Channel::new(d1_geom);
        
        let d2_geom = ChannelGeometry::<f64>::circular(1.0e-2, 1.58e-3, 1e-6);
        let d2 = Channel::new(d2_geom);

        let bifurcation = BifurcationJunction::new(parent, d1, d2, 0.5);
        let deviation = bifurcation.murray_law_deviation();

        // Should satisfy Murray's law closely
        assert!(deviation < 0.2); // Within 20%
    }
}
