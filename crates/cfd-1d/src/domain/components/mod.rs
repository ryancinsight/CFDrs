//! Microfluidic and millifluidic components for 1D network modeling.
//!
//! This module provides models for various microfluidic components including
//! channels, pumps, valves, mixers, and sensors with their characteristic
//! hydraulic properties and behaviors.

use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::ConstantPropertyFluid;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use std::collections::HashMap;

// Re-export submodules
pub mod channels;
pub mod constants;
pub mod factory;
pub mod membranes;
pub mod mixers;
pub mod pumps;
pub mod sensors;
pub mod valves;

// Re-export commonly used types
pub use channels::{CircularChannel, RectangularChannel};
pub use factory::ComponentFactory;
pub use membranes::{OrganCompartment, PorousMembrane};
pub use mixers::{Micromixer, MixerType};
pub use pumps::{Micropump, PumpType};
pub use sensors::{FlowSensor, SensorType};
pub use valves::{Microvalve, ValveType};

/// Trait for all microfluidic components
pub trait Component<T: RealField + Copy> {
    /// Get the hydraulic resistance of the component (linear part R)
    fn resistance(&self, fluid: &ConstantPropertyFluid<T>) -> T;

    /// Get the linear (R) and quadratic (k) resistance coefficients
    /// such that ΔP = R·Q + k·Q|Q|
    fn coefficients(&self, fluid: &ConstantPropertyFluid<T>) -> (T, T) {
        (self.resistance(fluid), T::zero())
    }

    /// Get the pressure drop across the component for a given flow rate
    fn pressure_drop(&self, flow_rate: T, fluid: &ConstantPropertyFluid<T>) -> T {
        let (r, k) = self.coefficients(fluid);
        r * flow_rate + k * flow_rate * flow_rate.abs()
    }

    /// Get the component type identifier
    fn component_type(&self) -> &str;

    /// Get component parameters
    fn parameters(&self) -> &HashMap<String, T>;

    /// Update component parameters
    fn set_parameter(&mut self, key: &str, value: T) -> Result<()>;

    /// Check if the component is active (e.g., pump, valve)
    fn is_active(&self) -> bool {
        false
    }

    /// Get the volume of the component (for transient analysis)
    fn volume(&self) -> Option<T> {
        None
    }
}

/// Convert a compile-time `f64` constant into a real scalar.
///
/// This is used for infallible component defaults where a panic would be a
/// worse failure mode than a saturated fallback for exotic scalar types.
#[inline]
pub(crate) fn real_from_f64<T>(value: f64) -> T
where
    T: RealField + Copy + FromPrimitive,
{
    T::from_f64(value).unwrap_or_else(|| {
        if value.is_sign_negative() {
            -(T::one() / T::zero())
        } else {
            T::one() / T::zero()
        }
    })
}

/// Convert a compile-time `f64` constant into a real scalar or return an error.
#[inline]
pub(crate) fn try_real_from_f64<T>(value: f64, context: &str) -> Result<T>
where
    T: RealField + Copy + FromPrimitive,
{
    T::from_f64(value).ok_or_else(|| {
        Error::InvalidConfiguration(format!(
            "{context}: constant {value} cannot be represented in the target scalar type"
        ))
    })
}
