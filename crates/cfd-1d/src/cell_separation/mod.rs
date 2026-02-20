//! Inertial microfluidic cell separation for millifluidic device design.
//!
//! This module provides physics-based models for predicting the lateral
//! equilibrium positions of cell populations in rectangular microchannels,
//! enabling design of devices that separate cancer cells from healthy blood
//! cells using inertial focusing and Dean flow.
//!
//! # Physical basis
//!
//! In a straight rectangular channel at finite Reynolds number, particles
//! experience inertial lift forces that drive them to stable equilibrium
//! positions.  Larger, stiffer particles (cancer cells) focus near the
//! channel center; smaller, more deformable particles (RBCs) focus near
//! the walls.  In curved channels, Dean flow secondary circulation enhances
//! this separation.
//!
//! # Module structure
//!
//! | Module | Contents |
//! |--------|----------|
//! | [`properties`] | Cell physical properties (size, deformability, density) |
//! | [`margination`] | Inertial lift and Dean drag force models |
//! | [`separation_model`] | High-level separation analysis and efficiency metrics |
//!
//! # Quick start
//!
//! ```rust,no_run
//! use cfd_1d::cell_separation::{CellProperties, CellSeparationModel};
//!
//! let cancer = CellProperties::mcf7_breast_cancer();
//! let healthy = CellProperties::red_blood_cell();
//!
//! let model = CellSeparationModel::new(500e-6, 200e-6, None);
//! let analysis = model.analyze(&cancer, &healthy, 1060.0, 3.5e-3, 0.05)
//!     .expect("cells must focus (κ > 0.07)");
//!
//! println!("Separation efficiency: {:.2}", analysis.separation_efficiency);
//! println!("Cancer center fraction: {:.2}", analysis.target_center_fraction);
//! ```
//!
//! # References
//! - Di Carlo, D. (2009). Inertial microfluidics. *Lab Chip*, 9, 3038–3046.
//! - Gossett, D. R. & Di Carlo, D. (2009). *Anal. Chem.*, 81, 8459–8465.
//! - Hur, S. C. et al. (2011). *Lab Chip*, 11, 912–920.

pub mod margination;
pub mod properties;
pub mod separation_model;

pub use margination::{
    dean_drag_force_n, dean_number, inertial_lift_force_n, lateral_equilibrium, EquilibriumResult,
};
pub use properties::CellProperties;
pub use separation_model::{CellSeparationAnalysis, CellSeparationModel};
