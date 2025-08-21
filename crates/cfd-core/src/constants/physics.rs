//! Physical constants for CFD calculations
//! 
//! All values are in SI units and referenced from standard literature:
//! - Fluid properties: White, F.M. (2011) "Fluid Mechanics", 7th Edition
//! - Numerical parameters: Ferziger & Perić (2002) "Computational Methods for Fluid Dynamics"

// Reynolds number regime boundaries
pub const REYNOLDS_LAMINAR_LIMIT: f64 = crate::constants::LAMINAR_THRESHOLD;  // Below this: laminar flow
pub const REYNOLDS_TRANSITION_START: f64 = crate::constants::LAMINAR_THRESHOLD;  // Start of transition
pub const REYNOLDS_TRANSITION_END: f64 = crate::constants::TURBULENT_THRESHOLD;  // End of transition
pub const REYNOLDS_TURBULENT_THRESHOLD: f64 = crate::constants::TURBULENT_THRESHOLD;  // Above this: turbulent

// Fluid properties at 20°C and 1 atm
pub const WATER_DENSITY_20C: f64 = 998.2;  // kg/m³
pub const WATER_VISCOSITY_20C: f64 = 1.002e-3;  // Pa·s (dynamic)
pub const WATER_KINEMATIC_VISCOSITY_20C: f64 = 1.004e-6;  // m²/s
pub const WATER_BULK_MODULUS: f64 = 2.2e9;  // Pa
pub const WATER_SPECIFIC_HEAT: f64 = 4182.0;  // J/(kg·K)
pub const WATER_THERMAL_CONDUCTIVITY: f64 = 0.598;  // W/(m·K)
pub const WATER_PRANDTL_NUMBER: f64 = 7.0;  // Pr = μ·Cp/k

pub const AIR_DENSITY_20C: f64 = 1.204;  // kg/m³
pub const AIR_VISCOSITY_20C: f64 = 1.825e-5;  // Pa·s (dynamic)
pub const AIR_KINEMATIC_VISCOSITY_20C: f64 = 1.516e-5;  // m²/s
pub const AIR_SPECIFIC_HEAT_CP: f64 = 1005.0;  // J/(kg·K)
pub const AIR_SPECIFIC_HEAT_CV: f64 = 718.0;  // J/(kg·K)
pub const AIR_GAMMA: f64 = 1.4;  // Cp/Cv
pub const AIR_GAS_CONSTANT: f64 = 287.0;  // J/(kg·K)
pub const AIR_THERMAL_CONDUCTIVITY: f64 = 0.0257;  // W/(m·K)
pub const AIR_PRANDTL_NUMBER: f64 = 0.71;  // Pr

// Universal physical constants
pub const GRAVITY_STANDARD: f64 = crate::constants::E_WALL_FUNCTION0665;  // m/s²
pub const ATMOSPHERIC_PRESSURE: f64 = 101325.0;  // Pa
pub const STEFAN_BOLTZMANN: f64 = 5.67e-8;  // W/(m²·K⁴)
pub const UNIVERSAL_GAS_CONSTANT: f64 = 8314.46;  // J/(kmol·K)

// Numerical solver parameters (Patankar 1980, Ferziger & Perić 2002)
pub const RELAXATION_FACTOR_VELOCITY: f64 = 0.7;  // Under-relaxation for velocity
pub const RELAXATION_FACTOR_PRESSURE: f64 = 0.3;  // Under-relaxation for pressure
pub const RELAXATION_FACTOR_TEMPERATURE: f64 = 0.8;  // Under-relaxation for temperature
pub const RELAXATION_FACTOR_TURBULENCE: f64 = 0.8;  // Under-relaxation for k-ε

pub const CONVERGENCE_TOLERANCE_VELOCITY: f64 = 1e-6;  // m/s
pub const CONVERGENCE_TOLERANCE_PRESSURE: f64 = 1e-6;  // Pa
pub const CONVERGENCE_TOLERANCE_CONTINUITY: f64 = 1e-5;  // kg/s
pub const CONVERGENCE_TOLERANCE_ENERGY: f64 = 1e-8;  // J

pub const MAX_ITERATIONS_OUTER: usize = 1000;  // Outer iterations
pub const MAX_ITERATIONS_INNER: usize = 100;   // Inner iterations (pressure correction)
pub const MAX_ITERATIONS_LINEAR_SOLVER: usize = 1000;  // Linear solver iterations

// CFL conditions
pub const CFL_NUMBER_EXPLICIT: f64 = 0.5;  // For explicit schemes
pub const CFL_NUMBER_IMPLICIT: f64 = 1.0;  // For implicit schemes
pub const CFL_NUMBER_MAXIMUM: f64 = 10.0;  // Maximum allowable

// Discretization scheme parameters
pub const PECLET_THRESHOLD: f64 = 2.0;  // Pe > 2: use upwind
pub const ASPECT_RATIO_LIMIT: f64 = 10.0;  // Maximum cell aspect ratio
pub const SKEWNESS_LIMIT: f64 = 0.95;  // Maximum cell skewness
pub const ORTHOGONALITY_LIMIT: f64 = 0.2;  // Minimum orthogonality

// Turbulence model constants (k-ε model, Launder & Spalding 1974)
pub const TURBULENCE_CMU: f64 = 0.09;
pub const TURBULENCE_C1E: f64 = 1.44;
pub const TURBULENCE_C2E: f64 = 1.92;
pub const TURBULENCE_SIGMA_K: f64 = 1.0;
pub const TURBULENCE_SIGMA_E: f64 = 1.3;

// Wall function parameters
pub const WALL_YPLUS_LAMINAR: f64 = crate::constants::Y_PLUS_LAMINAR;  // y+ < crate::constants::Y_PLUS_LAMINAR: viscous sublayer
pub const WALL_YPLUS_BUFFER_START: f64 = 5.0;  // Start of buffer layer
pub const WALL_YPLUS_BUFFER_END: f64 = 30.0;  // End of buffer layer
pub const WALL_KAPPA: f64 = crate::constants::VON_KARMAN;  // von Kármán constant
pub const WALL_E: f64 = crate::constants::E_WALL_FUNCTION;  // Wall roughness parameter

// Lattice Boltzmann D2Q9 weights (Sukop & Thorne 2007)
pub const D2Q9_WEIGHT_CENTER: f64 = 4.0 / 9.0;
pub const D2Q9_WEIGHT_CARDINAL: f64 = 1.0 / 9.0;
pub const D2Q9_WEIGHT_DIAGONAL: f64 = 1.0 / 36.0;

// Numerical precision
pub const EPSILON_MACHINE: f64 = 2.220446049250313e-16;  // Machine epsilon for f64
pub const EPSILON_TOLERANCE: f64 = 1e-10;  // General tolerance
pub const EPSILON_ZERO: f64 = 1e-14;  // Near-zero threshold