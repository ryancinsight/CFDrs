//! Smagorinsky Large Eddy Simulation (LES) model
//!
//! The Smagorinsky model is one of the most widely used SGS (Sub-Grid Scale)
//! models for LES. It assumes that the SGS stress tensor is proportional to
//! the magnitude of the resolved strain rate tensor.
//!
//! ## Mathematical Formulation
//!
//! The SGS stress tensor is modeled as:
//!
//! τ_ij = -2ν_SGS S_ij
//!
//! where ν_SGS = (C_S Δ)² |S| is the SGS viscosity,
//! C_S is the Smagorinsky constant (typically 0.1-0.2),
//! Δ is the filter width, and |S| is the magnitude of the strain rate.
//!
//! ## Dynamic Procedure
//!
//! The dynamic Smagorinsky model computes C_S locally using:
//!
//! C_S² = <L_ij M_ij> / <M_ij M_ij>
//!
//! where L_ij and M_ij are computed from resolved scales.
//!
//! ## References
//!
//! - Smagorinsky, J. (1963). General circulation experiments with the primitive equations.
//! - Germano, M., et al. (1991). A dynamic subgrid-scale eddy viscosity model.
//! - Lilly, D. K. (1992). A proposed modification of the Germano subgrid-scale closure method.
//!
//! ## WALE Model
//!
//! The Wall-Adapting Local Eddy-viscosity (WALE) model provides better near-wall
//! behavior than the Smagorinsky model by using a different invariant of the
//! velocity gradient tensor.
//!
//! **SGS Viscosity (WALE):**
//! ```math
//! \nu_{sgs} = (C_w \Delta)^2 \frac{(S_{ij}^d S_{ij}^d)^{3/2}}{(\bar{S}_{ij} \bar{S}_{ij})^{5/2} + (S_{ij}^d S_{ij}^d)^{5/4}}
//! ```
//!
//! where $S_{ij}^d$ is the traceless symmetric part of the square of the velocity gradient tensor.
//!
//! **References:**
//! - Nicoud, F. & Ducros, F. (1999). Subgrid-scale stress modelling based on the square of the velocity gradient tensor.
//! - Ducros, F., et al. (1998). Wall-adapting local eddy-viscosity model.

pub mod config;
pub mod dynamic;
pub mod gpu;
pub mod miles;
pub mod model;
pub mod sigma;
pub mod strain;
pub mod viscosity;
pub mod vreman;
pub mod wale;

// Re-export main types
pub use config::SmagorinskyConfig;
pub use miles::MilesLES;
pub use model::SmagorinskyLES;
pub use sigma::SigmaModel;
pub use vreman::VremanModel;
pub use wale::WaleModel;
