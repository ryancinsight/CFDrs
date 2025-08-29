//! Turbulence model constants

/// von Kármán constant
pub const KAPPA: f64 = cfd_core::constants::fluid::VON_KARMAN;

/// Roughness parameter for smooth walls
pub const E_WALL_FUNCTION: f64 = cfd_core::constants::fluid::WALL_FUNCTION_E;

/// k-ε model constants
pub const C_MU: f64 = 0.09;
pub const C1_EPSILON: f64 = 1.44;
pub const C2_EPSILON: f64 = 1.92;
pub const SIGMA_K: f64 = 1.0;
pub const SIGMA_EPSILON: f64 = 1.3;

/// SST model constants (Menter 1994)
pub const SST_ALPHA_1: f64 = 0.31;
pub const SST_BETA_1: f64 = 0.075;
pub const SST_BETA_2: f64 = 0.0828;
pub const SST_BETA_STAR: f64 = 0.09;
/// γ₁ = β₁/β* - σω₁κ²/√β* (Menter 1994 Eq. 11)
pub const SST_GAMMA_1: f64 = 0.5532; // = 0.075/0.09 - 0.5*0.41²/√0.09
pub const SST_GAMMA_2: f64 = 0.4403; // = 0.0828/0.09 - 0.856*0.41²/√0.09
pub const SST_SIGMA_K1: f64 = 0.85;
pub const SST_SIGMA_K2: f64 = 1.0;
pub const SST_SIGMA_OMEGA1: f64 = 0.5;
pub const SST_SIGMA_OMEGA2: f64 = 0.856;

/// Numerical stability
pub const EPSILON_MIN: f64 = 1e-10;
pub const OMEGA_MIN: f64 = 1e-10;

/// Wall treatment thresholds
pub const Y_PLUS_VISCOUS_SUBLAYER: f64 = 5.0;
pub const Y_PLUS_LOG_LAW: f64 = cfd_core::constants::fluid::Y_PLUS_LAMINAR;
pub const Y_PLUS_BUFFER_START: f64 = 5.0;
pub const Y_PLUS_BUFFER_END: f64 = 30.0;

/// Wall function coefficients
pub const K_VISC_COEFFICIENT: f64 = 11.0;
pub const OMEGA_WALL_COEFFICIENT: f64 = 60.0;
pub const BLENDING_FACTOR: f64 = 0.01;

/// Common numerical constants
pub const HALF: f64 = 0.5;
pub const TWO: f64 = 2.0;
pub const THREE: f64 = 3.0;
pub const FOUR: f64 = 4.0;
pub const FIVE: f64 = 5.0;
