use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationRow {
    pub track: String,
    pub id: String,
    pub topology: String,
    pub dp_1d_bernoulli_pa: f64,
    pub dp_2d_fvm_pa: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dp_3d_fem_pa: Option<f64>,
    pub agreement_1d_2d_pct: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub agreement_2d_3d_pct: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mass_error_3d_pct: Option<f64>,
    pub sigma_1d: f64,
    pub sigma_2d: f64,
    pub score: f64,
    /// Whether the 2D FVM SIMPLE solver converged and the velocity field is physical.
    #[serde(default)]
    pub two_d_converged: bool,
    /// Whether the throat Reynolds number exceeds the Stokes-flow validity limit (Re > 1).
    /// When true the 3D Stokes FEM ΔP will be << Bernoulli (inertia-dominated flow).
    #[serde(default)]
    pub high_re_stokes_mismatch: bool,
}
