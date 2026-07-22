//! `cfd-3d` cavitation physics.  See
//! `repos/CFDrs/docs/book/turbulence_multiphase.md` § "Cavitation Physics"
//! for the integration contract.
//!
//! Stub module — phase-transition model whose rate depends on local
//! pressure vs. vapor pressure.  Default is Rayleigh-Plesset
//! (single-bubble); Eulerian-Eulerian is a cloud-cavitation variant.
//! Full implementation is filled in by the cfd-3d FEM team; this
//! file is intentionally minimal so the link target resolves while
//! the FEM team whitelists `pub mod cavitation;` in `lib.rs`.

/// Cavitation-model trait.  Inception pressure + collapse rate.
///
/// Mirrors the `Cavitation` trait spelled out in
/// `turbulence_multiphase.md` § "Cavitation Physics".
pub trait Cavitation {
    /// Local cavity inception pressure (Pa).
    fn inception_pressure(&self) -> f64;
    /// Collapse rate (kg/m²/s²) at local pressure `p` (Pa).
    fn collapse_rate(&self, p: f64) -> f64;
}

/// Rayleigh-Plesset closure (default for single-bubble dynamics).
#[derive(Default, Clone, Debug)]
pub struct RayleighPlesset {
    /// Vapor pressure of the working fluid (Pa).  Defaults to water at ~20 °C.
    pub vapor_pressure: f64,
}

/// Eulerian-Eulerian closure for cloud cavitation.
#[derive(Default, Clone, Debug)]
pub struct EulerianEulerian {
    /// Threshold nucleation density (bubbles/m³).
    pub nucleation_density: f64,
}

/// Accumulate cavitation damage for one timestep.
///
/// Reports a delta in `kg/m²/s` so downstream fatigue models can
/// post-process without re-scaling.  Mirrors the chapter:
///
/// ```text
/// damage += (inception_indicator) · (collapse_rate · dt)
/// ```
pub fn damage_step(rate: f64, dt: f64) -> f64 {
    rate * dt
}

impl Cavitation for RayleighPlesset {
    fn inception_pressure(&self) -> f64 {
        self.vapor_pressure
    }
    fn collapse_rate(&self, p: f64) -> f64 {
        // TODO(cfd-3d, atlas-team): Rayleigh-Plesset bubble dynamics.
        let _ = p;
        0.0
    }
}

impl Cavitation for EulerianEulerian {
    fn inception_pressure(&self) -> f64 {
        // Cloud-cavitation inception ≠ single-bubble vapor pressure; scale by density.
        0.0_f64.max(self.nucleation_density.sqrt() * 1.0e3)
    }
    fn collapse_rate(&self, p: f64) -> f64 {
        // TODO(cfd-3d, atlas-team): Eulerian-Eulerian collapse rate.
        let _ = p;
        0.0
    }
}
