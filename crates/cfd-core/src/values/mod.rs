//! Value objects for CFD domain modeling.
//!
//! This module provides immutable value objects that represent important
//! domain concepts with built-in validation and behavior.

pub mod dimensionless;
pub mod flow;
pub mod pressure;
pub mod temperature;
pub mod velocity;

pub use dimensionless::{DimensionlessNumber, DimensionlessType};
pub use flow::{FlowGeometry, ReynoldsNumber};
pub use pressure::{Pressure, PressureUnit};
pub use temperature::{Temperature, TemperatureUnit};
pub use velocity::Velocity;
