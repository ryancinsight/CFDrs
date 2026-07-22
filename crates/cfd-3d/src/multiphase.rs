//! `cfd-3d` multiphase coupling.  See
//! `repos/CFDrs/docs/book/turbulence_multiphase.md` § "Multiphase Coupling"
//! for the integration contract.
//!
//! Stub module — volume-of-fluid style coupling for liquid+gas mixtures.
//! Selects between co-located and staggered momentum schemes.  Full
//! implementation is filled in by the cfd-3d FEM team; this file is
//! intentionally minimal so the link target resolves while the FEM
//! team whitelists `pub mod multiphase;` in `lib.rs`.

/// Momentum-coupling scheme selector (co-located vs staggered mesh).
///
/// Mirrors the `MomentumCoupling` enum spelled out in
/// `turbulence_multiphase.md` § "Multiphase Coupling".
#[derive(Default, Clone, Copy, Debug, PartialEq, Eq)]
pub enum MomentumCoupling {
    /// Co-located velocity-pressure grid.
    #[default]
    Colocated,
    /// Staggered (Marker-and-Cell) grid.
    Staggered,
}

/// Volume-weighted mixture properties from a phase volume fraction.
///
/// `alpha`: gas-phase volume fraction ∈ [0, 1].
/// `rho_l`, `rho_g`: liquid / gas densities.
/// `mu_l`, `mu_g`: liquid / gas dynamic viscosities.
///
/// Returns `(rho_mix, mu_mix)`.
///
/// Mirrors the chapter's pseudocode verbatim:
/// ```text
/// let rho_mix = rho_l * (1.0 - alpha) + rho_g * alpha;
/// let mu_mix  = mu_l  * (1.0 - alpha) + mu_g  * alpha;
/// ```
pub fn exchange(alpha: f64, rho_l: f64, rho_g: f64, mu_l: f64, mu_g: f64) -> (f64, f64) {
    let rho_mix = rho_l * (1.0 - alpha) + rho_g * alpha;
    let mu_mix = mu_l * (1.0 - alpha) + mu_g * alpha;
    (rho_mix, mu_mix)
}
