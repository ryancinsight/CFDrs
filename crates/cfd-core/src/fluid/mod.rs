//! Fluid properties and models module.

pub mod properties;
pub mod temperature;
pub mod viscosity;

pub use properties::{Fluid, FluidBuilder};
pub use temperature::TemperatureDependence;
pub use viscosity::{ViscosityModel, ViscosityCalculator};