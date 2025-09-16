//! Turbulence model constants

/// von Kármán constant
pub const KAPPA: f64 = cfd_core::constants::fluid::VON_KARMAN;

/// Roughness parameter for smooth walls
pub const E_WALL_FUNCTION: f64 = cfd_core::constants::fluid::WALL_FUNCTION_E;

/// k-ε model constants
pub const C_MU: f64 = 0.09;
/// First closure coefficient for turbulent kinetic energy production (Launder & Spalding 1974)
pub const C1_EPSILON: f64 = 1.44;
/// Second closure coefficient for turbulent kinetic energy dissipation (Launder & Spalding 1974)
pub const C2_EPSILON: f64 = 1.92;
/// Turbulent kinetic energy Schmidt number (Launder & Spalding 1974)
pub const SIGMA_K: f64 = 1.0;
/// Turbulent dissipation rate Schmidt number (Launder & Spalding 1974)
pub const SIGMA_EPSILON: f64 = 1.3;

/// SST model constants (Menter 1994)
pub const SST_ALPHA_1: f64 = 0.31;
/// First closure coefficient for turbulent kinetic energy dissipation (Menter 1994)
pub const SST_BETA_1: f64 = 0.075;
/// Second closure coefficient for turbulent kinetic energy dissipation (Menter 1994)
pub const SST_BETA_2: f64 = 0.0828;
/// Universal closure coefficient for turbulent kinetic energy (Menter 1994)
pub const SST_BETA_STAR: f64 = 0.09;
/// γ₁ = β₁/β* - σω₁κ²/√β* (Menter 1994 Eq. 11)
pub const SST_GAMMA_1: f64 = 0.5532; // = 0.075/0.09 - 0.5*0.41²/√0.09
/// γ₂ = β₂/β* - σω₂κ²/√β* (Menter 1994 Eq. 11)
pub const SST_GAMMA_2: f64 = 0.4403; // = 0.0828/0.09 - 0.856*0.41²/√0.09
/// Turbulent kinetic energy diffusion coefficient for inner region (Menter 1994)
pub const SST_SIGMA_K1: f64 = 0.85;
/// Turbulent kinetic energy diffusion coefficient for outer region (Menter 1994)
pub const SST_SIGMA_K2: f64 = 1.0;
/// Specific turbulent dissipation rate diffusion coefficient for inner region (Menter 1994)
pub const SST_SIGMA_OMEGA1: f64 = 0.5;
/// Specific turbulent dissipation rate diffusion coefficient for outer region (Menter 1994)
pub const SST_SIGMA_OMEGA2: f64 = 0.856;

/// Numerical stability
pub const EPSILON_MIN: f64 = 1e-10;
/// Minimum allowable specific turbulent dissipation rate for numerical stability
pub const OMEGA_MIN: f64 = 1e-10;

/// Wall treatment thresholds
pub const Y_PLUS_VISCOUS_SUBLAYER: f64 = 5.0;
/// Dimensionless wall distance threshold for log-law region onset
pub const Y_PLUS_LOG_LAW: f64 = cfd_core::constants::fluid::Y_PLUS_LAMINAR;
/// Dimensionless wall distance for buffer layer start (Spalding 1961)
pub const Y_PLUS_BUFFER_START: f64 = 5.0;
/// Dimensionless wall distance for buffer layer end (Pope 2000)
pub const Y_PLUS_BUFFER_END: f64 = 30.0;

/// Wall function coefficients
pub const K_VISC_COEFFICIENT: f64 = 11.0;
/// Wall boundary condition coefficient for specific dissipation rate (Menter 1994)
pub const OMEGA_WALL_COEFFICIENT: f64 = 60.0;
/// Blending factor for SST model near-wall treatment (Menter 1994)
pub const BLENDING_FACTOR: f64 = 0.01;
