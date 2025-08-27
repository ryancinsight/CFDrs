//! Physical constants for microfluidic components

/// Default surface roughness for smooth channels [m]
pub const DEFAULT_ROUGHNESS: f64 = 1e-6;

/// Minimum Reynolds number for laminar flow
pub const RE_LAMINAR_MIN: f64 = 0.1;

/// Maximum Reynolds number for laminar flow
pub const RE_LAMINAR_MAX: f64 = cfd_core::constants::dimensionless::reynolds::PIPE_CRITICAL_LOWER;

/// Transition Reynolds number
pub const RE_TRANSITION: f64 = cfd_core::constants::dimensionless::reynolds::PIPE_CRITICAL_LOWER;

/// Minimum Reynolds number for turbulent flow
pub const RE_TURBULENT_MIN: f64 = cfd_core::constants::dimensionless::reynolds::PIPE_CRITICAL_UPPER;

/// Default pump efficiency
pub const DEFAULT_PUMP_EFFICIENCY: f64 = 0.7;

/// Default valve flow coefficient
pub const DEFAULT_VALVE_CV: f64 = 0.1;

/// Default mixing efficiency
pub const DEFAULT_MIXING_EFFICIENCY: f64 = 0.95;

/// Friction factor constants for rectangular channels
pub const RECT_CHANNEL_C1: f64 = 24.0;
pub const RECT_CHANNEL_C2: f64 = 64.0;
pub const RECT_CHANNEL_C3: f64 = 96.0;

/// Friction factor for circular channels (laminar)
pub const CIRCULAR_FRICTION_LAMINAR: f64 = 64.0;

/// Colebrook-White equation tolerance
pub const COLEBROOK_TOLERANCE: f64 = 1e-6;

/// Maximum iterations for iterative solvers
pub const MAX_ITERATIONS: usize = 100;

/// Laminar flow friction factor coefficient for circular pipes (f = 64/Re)
pub const LAMINAR_FRICTION_COEFFICIENT: f64 = 64.0;

/// Hydraulic diameter factor for rectangular channels
pub const HYDRAULIC_DIAMETER_FACTOR: f64 = 4.0;

/// Diameter exponent for Hagen-Poiseuille law
pub const HAGEN_POISEUILLE_EXPONENT: f64 = 4.0;

/// Factor of 2 for perimeter calculation
pub const PERIMETER_FACTOR: f64 = 2.0;
