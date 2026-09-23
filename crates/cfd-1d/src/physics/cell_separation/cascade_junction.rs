#![cfg_attr(test, expect(clippy::print_stderr, reason = "test diagnostic output"))]
//! Zweifach-Fung junction routing model for selective split-sequence trees.
//!
//! Used by primitive selective branching topologies to predict the
//! cell-type-specific distribution between the designated treatment arm
//! (cancer/WBC-enriched, routed to the therapy zone) and the peripheral
//! bypass arms (RBC-enriched, low shear).
//!
//! # Physical basis
//!
//! At a trifurcation junction, the **Zweifach-Fung law** (1969) predicts that large,
//! stiff particles (cancer cells, WBCs) preferentially enter the arm carrying the
//! highest volumetric flow.  The probability follows a power law:
//!
//! ```text
//! P_center(cell) = r_c^beta / (r_c^beta + 2*r_p^beta)
//! ```
//!
//! where `r_c = Q_center / Q_total`, `r_p = (1 - r_c) / 2` (each peripheral
//! arm), and `beta` is a stiffness exponent:
//!
//! | Cell type | beta | Basis |
//! |-----------|------|-------|
//! | Cancer (MCF-7, ~17.5 um, stiff)   | 1.85 | Size-enhanced stiff sphere (Hou 2012, Karabacak 2014) |
//! | WBC (~10 um, semi-rigid)           | 1.40 | Intermediate deformability |
//! | RBC (~7 um, highly deformable)     | 1.00 | Deformable cell - flow-weighted |
//!
//! The cascade iterates this routing over `n_levels` trifurcation levels, each time
//! narrowing the center arm by `center_frac` and carrying `q_center_frac` of the
//! current level's flow.  After all levels the center-arm fraction of each cell type
//! accumulates multiplicatively.
//!
//! # Module structure
//!
//! This module root is a passthrough facade: it declares the routing leaves and
//! the shared descriptor types, then re-exports their public surface so every
//! path that existed before the split still resolves.
//!
//! | Module | Contents |
//! |--------|----------|
//! | [`types`] | Shared cascade descriptors (`CascadeStage`, `PeripheralRecovery`, results) |
//! | [`routing_probability`] | Core Zweifach-Fung probability functions and cell constants |
//! | [`cascade_routing`] | Cascade trifurcation, cross-junction, and mixed Bi/Tri routing |
//! | [`incremental_filtration`] | CIF staged selective-routing with pre-trifurcation skimming |
//!
//! # References
//! - Fung, Y. C. (1969). Biorheology of soft tissues. *Biorheology*, 6, 409-419.
//! - Yang et al. (2017). Red blood cell phase separation in symmetric and asymmetric microchannel
//!   networks: effect of capillary dilation and inflow velocity. *PMC5114676*.
//! - Doyeux, V. et al. (2011). Spheres in the vicinity of a bifurcation: elucidating
//!   the Zweifach-Fung effect. *J. Fluid Mech.*, 674, 359-388.
//! - Di Carlo, D. (2009). Inertial microfluidics. *Lab Chip*, 9, 3038-3046.

mod types;
pub mod cascade_routing;
pub mod incremental_filtration;
pub mod routing_probability;

pub use types::{
    CascadeJunctionResult, CascadeStage, IncrementalFiltrationResult, PeripheralRecovery,
};

pub use cascade_routing::{
    cascade_junction_separation, cascade_junction_separation_cross_junction,
    cascade_junction_separation_from_qfracs, checked_cascade_junction_separation,
    checked_cascade_junction_separation_cross_junction,
    checked_cascade_junction_separation_from_qfracs, checked_mixed_cascade_separation,
    checked_mixed_cascade_separation_kappa_aware, checked_treatment_bifurcation_separation,
    checked_tri_asymmetric_q_fracs, checked_tri_center_q_frac,
    checked_tri_center_q_frac_cross_junction, mixed_cascade_separation,
    mixed_cascade_separation_kappa_aware, treatment_bifurcation_separation, tri_asymmetric_q_fracs,
    tri_center_q_frac, tri_center_q_frac_cross_junction,
};
pub use incremental_filtration::{
    checked_cif_pretri_stage_center_fracs, checked_cif_pretri_stage_q_fracs,
    checked_cif_pretri_stage_q_fracs_cross_junction,
    checked_incremental_filtration_separation_cross_junction,
    checked_incremental_filtration_separation_from_qfracs,
    checked_incremental_filtration_separation_staged, cif_pretri_stage_center_fracs,
    cif_pretri_stage_q_fracs, cif_pretri_stage_q_fracs_cross_junction,
    incremental_filtration_separation_cross_junction,
    incremental_filtration_separation_from_qfracs, incremental_filtration_separation_staged,
};

#[cfg(test)]
mod tests;
