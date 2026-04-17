//! Hydraulic resistance models for 1D CFD networks.
//!
//! This module provides comprehensive resistance modeling for various
//! microfluidic components and flow conditions, including analytical
//! solutions and empirical correlations.
//!
//! ## Theorem — Hagen-Poiseuille (1839/1840)
//!
//! **Theorem**: For incompressible, steady, laminar flow of a Newtonian fluid in
//! a straight circular tube of radius R and length L, the exact solution of the
//! Navier-Stokes equations gives:
//!
//! ```text
//! Q = π R⁴ ΔP / (8 μ L) = ΔP / R_HP
//! ```
//!
//! where R_HP = 128 μ L / (π D⁴) is the hydraulic resistance.
//!
//! **Proof sketch**: The NS equations reduce to ∇²u = (1/μ) dP/dz in cylindrical
//! coordinates. Applying no-slip BC u(R) = 0 and symmetry du/dr(0) = 0 gives the
//! parabolic profile u(r) = (ΔP/4μL)(R²-r²). Integration over the cross-section
//! yields the volumetric flow rate Q = π R⁴ ΔP / (8μL).
//!
//! **Validity Conditions**: Re < 2300 (laminar), L/D > 0.06 Re (fully developed),
//! Kn < 0.01 (no-slip continuum).
//!
//! ## Theorem — Rectangular Duct Resistance (Shah & London 1978)
//!
//! **Theorem**: For a rectangular channel of width w and height h (aspect ratio
//! α = h/w ≤ 1), the hydraulic resistance via the series expansion in α is:
//!
//! ```text
//! R_rect = 12 μ L / (w h³) · 1/φ(α)     where
//! φ(α) = 1 - (192α/π⁵) Σₙ (tanh(n π/(2α)) / n⁵)  (n = 1, 3, 5, ...)
//! ```
//!
//! The series converges exponentially; truncating at N=5 odd terms gives relative
//! error < 10⁻⁶ for all α ∈ (0,1].
//!
//! ## Theorem — Dean Number Correction (Ito 1959)
//!
//! **Theorem**: For curved channels with centerline radius of curvature R_c and
//! hydraulic diameter D_h, the resistance is augmented by secondary Dean vortices:
//!
//! ```text
//! R_curved = R_HP · f(De)     where De = Re · √(D_h / 2R_c)
//! ```
//!
//! For De < 11.6 (weak curvature): f(De) = 1 + 0.033 (log₁₀ De)⁴.
//! For De > 11.6 (strong curvature): f(De) ≈ 0.1033 De^{0.5}.
//!
//! ## Theorem — Kirchhoff's Laws for Hydraulic Networks
//!
//! **Theorem (Mass Conservation at Junction)**:
//!
//! ```text
//! Σ Q_in = Σ Q_out    (conservation of volume for incompressible flow)
//! ```
//!
//! **Theorem (Pressure Continuity at Junction)**:
//! The pressure at a junction node is unique — all channels meeting at a node
//! share the same pressure value (no pressure jump in the absence of a membrane).
//!
//! ## Theorem — Knudsen Number Slip-Flow Regime
//!
//! **Theorem**: The Knudsen number Kn = λ/L determines flow regime:
//! - Kn < 0.001: Continuum, Hagen-Poiseuille applies.
//! - 0.001 < Kn < 0.1: Slip-flow; velocity slip at wall corrects Q by factor (1 + 4Kn).
//! - Kn > 0.1: Transition/molecular regime; continuum models fail.
//!
//! **Implementation**: The `NumericalParameters::Kn` field enforces this check and
//! selects the appropriate resistance model variant.

pub mod calculator;
pub mod factory;
pub mod geometry;
pub mod models;
pub mod traits;

// Re-export main types
pub use calculator::ResistanceCalculator;
pub use factory::ResistanceModelFactory;
pub use geometry::ChannelGeometry;
pub use models::CombinationMethod;
pub use models::{
    bayat_rezai_enhancement, cascade_treatment_flow_fractions, durst_entrance_k,
    durst_entrance_length, durst_resistance_multiplier, parallel_channel_flow_fractions, BendType,
    DarcyWeisbachModel, EntranceEffectsModel, ExpansionType, FlowConditions, HagenPoiseuilleModel,
    JunctionFlowDirection, JunctionLossModel, JunctionType, MembranePoreModel,
    RectangularChannelModel, ResistanceModel, SerpentineAnalysis, SerpentineCrossSection,
    SerpentineModel, VenturiAnalysis, VenturiGeometry, VenturiModel,
};

// Convenience re-export of the traits facade
pub use traits::{FlowConditions as TraitFlowConditions, ResistanceModel as ResistanceModelTrait};
