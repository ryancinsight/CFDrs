//! Physical constants for microfluidic components

/// Default surface roughness for smooth channels [m]
pub const DEFAULT_ROUGHNESS: f64 = 1e-6;

/// Minimum Reynolds number for laminar flow
pub const RE_LAMINAR_MIN: f64 = 0.1;

/// Maximum Reynolds number for laminar flow
pub const RE_LAMINAR_MAX: f64 = 2300.0;

/// Transition Reynolds number
pub const RE_TRANSITION: f64 = 2300.0;

/// Minimum Reynolds number for turbulent flow
pub const RE_TURBULENT_MIN: f64 = 4000.0;

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