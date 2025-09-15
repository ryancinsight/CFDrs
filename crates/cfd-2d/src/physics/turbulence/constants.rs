//! Turbulence model constants

/// von Kármán constant
pub const KAPPA: f64 = cfd_core::constants::fluid::VON_KARMAN;

/// Roughness parameter for smooth walls
pub const E_WALL_FUNCTION: f64 = cfd_core::constants::fluid::WALL_FUNCTION_E;

/// k-ε model constants

/// C_μ constant in k-ε model (Launder & Spalding 1974)
/// Used in eddy viscosity calculation: μ_t = ρ C_μ k²/ε
pub const C_MU: f64 = 0.09;

/// C₁ε constant in ε equation source term (Launder & Spalding 1974)
/// Appears in production term coefficient
pub const C1_EPSILON: f64 = 1.44;

/// C₂ε constant in ε equation destruction term (Launder & Spalding 1974)
/// Appears in dissipation term coefficient
pub const C2_EPSILON: f64 = 1.92;

/// σ_k diffusion coefficient for turbulent kinetic energy equation
/// Controls diffusion of k: ∇·((μ + μ_t/σ_k)∇k)
pub const SIGMA_K: f64 = 1.0;

/// σ_ε diffusion coefficient for dissipation rate equation
/// Controls diffusion of ε: ∇·((μ + μ_t/σ_ε)∇ε)
pub const SIGMA_EPSILON: f64 = 1.3;

/// SST model constants (Menter 1994)

/// α₁ constant for SST blending function calculation
/// Used in F₁ blending function for near-wall treatment
pub const SST_ALPHA_1: f64 = 0.31;

/// β₁ coefficient for inner k-ω model (Menter 1994)
/// Destruction term in ω equation for near-wall region
pub const SST_BETA_1: f64 = 0.075;

/// β₂ coefficient for outer k-ε model (Menter 1994)
/// Destruction term in ω equation for outer region
pub const SST_BETA_2: f64 = 0.0828;

/// β* coefficient in SST model (Menter 1994)
/// Appears in k equation destruction term: β* ρ ω k
pub const SST_BETA_STAR: f64 = 0.09;
/// γ₁ = β₁/β* - σω₁κ²/√β* (Menter 1994 Eq. 11)
pub const SST_GAMMA_1: f64 = 0.5532; // = 0.075/0.09 - 0.5*0.41²/√0.09
/// γ₂ coefficient for outer k-ε model (Menter 1994 Eq. 11)
/// γ₂ = β₂/β* - σω₂κ²/√β* = 0.0828/0.09 - 0.856*0.41²/√0.09
pub const SST_GAMMA_2: f64 = 0.4403; // = 0.0828/0.09 - 0.856*0.41²/√0.09

/// σ_k1 diffusion coefficient for k equation in inner region
/// Controls turbulent kinetic energy diffusion near walls
pub const SST_SIGMA_K1: f64 = 0.85;

/// σ_k2 diffusion coefficient for k equation in outer region
/// Controls turbulent kinetic energy diffusion in freestream
pub const SST_SIGMA_K2: f64 = 1.0;

/// σ_ω1 diffusion coefficient for ω equation in inner region
/// Controls specific dissipation rate diffusion near walls
pub const SST_SIGMA_OMEGA1: f64 = 0.5;

/// σ_ω2 diffusion coefficient for ω equation in outer region
/// Controls specific dissipation rate diffusion in freestream
pub const SST_SIGMA_OMEGA2: f64 = 0.856;

/// Numerical stability

/// Minimum value for dissipation rate ε to prevent division by zero
/// Prevents numerical instability in eddy viscosity calculation
pub const EPSILON_MIN: f64 = 1e-10;

/// Minimum value for specific dissipation rate ω to prevent division by zero
/// Prevents numerical instability in SST model calculations
pub const OMEGA_MIN: f64 = 1e-10;

/// Wall treatment thresholds

/// y⁺ threshold for viscous sublayer (y⁺ < 5)
/// Below this value, linear velocity profile applies
pub const Y_PLUS_VISCOUS_SUBLAYER: f64 = 5.0;

/// y⁺ threshold for log-law region transition
/// References core fluid mechanics constants
pub const Y_PLUS_LOG_LAW: f64 = cfd_core::constants::fluid::Y_PLUS_LAMINAR;

/// y⁺ value where buffer layer begins
/// Transition between viscous sublayer and log-law region
pub const Y_PLUS_BUFFER_START: f64 = 5.0;

/// y⁺ value where buffer layer ends
/// Full log-law region begins above this value
pub const Y_PLUS_BUFFER_END: f64 = 30.0;

/// Wall function coefficients

/// Coefficient for near-wall k boundary condition
/// Used in viscous sublayer treatment for turbulent kinetic energy
pub const K_VISC_COEFFICIENT: f64 = 11.0;

/// Coefficient for near-wall ω boundary condition (Menter 1994)
/// ω_wall = 60 μ/(ρ β₁ (Δy₁)²) where Δy₁ is first cell height
pub const OMEGA_WALL_COEFFICIENT: f64 = 60.0;

/// Blending factor for smooth transition between wall treatments
/// Used to blend between different near-wall formulations
pub const BLENDING_FACTOR: f64 = 0.01;
