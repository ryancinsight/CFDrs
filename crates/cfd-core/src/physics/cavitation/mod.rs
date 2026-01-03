//! Cavitation physics models for CFD simulations
//!
//! This module implements hydrodynamic cavitation models based on established
//! literature including Brennen (1995), Franc & Michel (2004), and Rayleigh-Plesset
//! dynamics for bubble growth and collapse.

pub mod constants;
pub mod damage;
pub mod models;
pub mod number;
pub mod rayleigh_plesset;
pub mod venturi;

// Re-export main types
pub use damage::CavitationDamage;
pub use models::CavitationModel;
pub use number::CavitationNumber;
pub use rayleigh_plesset::RayleighPlesset;
pub use venturi::VenturiCavitation;
