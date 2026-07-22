//! `cfd-3d` turbulence closures: k-ε, k-ω SST, Smagorinsky.
//!
//! Stub module — chapter-derived public API.  See
//! `repos/CFDrs/docs/book/turbulence_multiphase.md` for the integration
//! contract.  Full physics implementations are filled in by the cfd-3d
//! FEM team per the chapter's equations; this file is intentionally
//! minimal so the link target resolves while the FEM team whitelists
//! `pub mod turbulence;` in `lib.rs` once the equations are merged.

/// Turbulence closure trait.  Each step returns turbulent viscosity `ν_t`
/// per cell so the momentum solver can apply it as an additional
/// diffusion term.
///
/// Mirrors the `TurbulenceModel` trait spelled out in
/// `turbulence_multiphase.md` § "Turbulence Closures".
pub trait TurbulenceModel {
    /// Scalar field element type (typically `f64`).
    type Field;

    /// Advance the closure by one timestep.  Mutates `nu_t` in place.
    ///
    /// `u`, `v`, `w`: cell-centered velocity components (read-only).
    /// `nu_t`: turbulent-viscosity buffer (mutated).
    /// `dt`: timestep size.
    fn step(
        &mut self,
        u: &[f64],
        v: &[f64],
        w: &[f64],
        nu_t: &mut [f64],
        dt: f64,
    );

    /// Human-readable closure name.
    fn label(&self) -> &'static str;
}

/// Standard k-ε closure.
#[derive(Default, Clone, Debug)]
pub struct KEpsilon;

/// k-ω SST closure (Menter 1994).
#[derive(Default, Clone, Debug)]
pub struct KOmegaSST;

/// Smagorinsky LES closure.
#[derive(Default, Clone, Debug)]
pub struct Smagorinsky;

impl TurbulenceModel for KEpsilon {
    type Field = f64;
    fn step(
        &mut self,
        _u: &[f64],
        _v: &[f64],
        _w: &[f64],
        _nu_t: &mut [f64],
        _dt: f64,
    ) {
        // TODO(cfd-3d, atlas-team): implement k-ε production + dissipation.
    }
    fn label(&self) -> &'static str {
        "k-ε"
    }
}

impl TurbulenceModel for KOmegaSST {
    type Field = f64;
    fn step(
        &mut self,
        _u: &[f64],
        _v: &[f64],
        _w: &[f64],
        _nu_t: &mut [f64],
        _dt: f64,
    ) {
        // TODO(cfd-3d, atlas-team): implement k-ω SST blend.
    }
    fn label(&self) -> &'static str {
        "k-ω SST"
    }
}

impl TurbulenceModel for Smagorinsky {
    type Field = f64;
    fn step(
        &mut self,
        _u: &[f64],
        _v: &[f64],
        _w: &[f64],
        _nu_t: &mut [f64],
        _dt: f64,
    ) {
        // TODO(cfd-3d, atlas-team): implement Smagorinsky LES model.
    }
    fn label(&self) -> &'static str {
        "Smagorinsky"
    }
}
