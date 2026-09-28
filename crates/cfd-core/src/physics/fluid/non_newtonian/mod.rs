//! Non-Newtonian fluid models.
//!
//! References:
//! - Bird, R.B., Armstrong, R.C., Hassager, O. (1987) "Dynamics of Polymeric Liquids"
//! - Chhabra, R.P., Richardson, J.F. (2008) "Non-Newtonian Flow and Applied Rheology"

/// Bingham plastic fluid model.
mod bingham;
/// Carreau–Yasuda fluid model for blood and polymer solutions.
mod carreau_yasuda;
/// Casson fluid model for blood.
mod casson;
/// Herschel–Bulkley generalised Bingham plastic model.
mod herschel_bulkley;
/// Power-law (Ostwald–de Waele) fluid model.
mod power_law;

pub use bingham::BinghamPlastic;
pub use carreau_yasuda::CarreauYasuda;
pub use casson::Casson;
pub use herschel_bulkley::HerschelBulkley;
pub use power_law::PowerLawFluid;

#[cfg(test)]
mod tests;
