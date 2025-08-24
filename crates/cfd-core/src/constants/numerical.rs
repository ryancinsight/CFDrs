//! Numerical method constants for CFD computations
//!
//! Constants for numerical solvers, discretization schemes, and convergence criteria.

/// Solver convergence and iteration parameters
pub mod solver {
    /// Convergence tolerance for iterative solvers
    pub const CONVERGENCE_TOLERANCE: f64 = 1e-6;
    
    /// Machine epsilon tolerance
    pub const EPSILON_TOLERANCE: f64 = 1e-10;
    
    /// Zero threshold for numerical comparisons
    pub const ZERO_THRESHOLD: f64 = 1e-14;
    
    /// Maximum iterations for outer loops
    pub const MAX_ITERATIONS_OUTER: usize = 1000;
    
    /// Maximum iterations for inner loops
    pub const MAX_ITERATIONS_INNER: usize = 100;
    
    /// Maximum iterations for linear solvers
    pub const MAX_ITERATIONS_LINEAR: usize = 1000;
    
    /// Residual reduction factor for convergence
    pub const RESIDUAL_REDUCTION_FACTOR: f64 = 1e-3;
}

/// Relaxation factors for iterative methods
pub mod relaxation {
    /// Velocity under-relaxation factor (SIMPLE algorithm)
    pub const VELOCITY: f64 = 0.7;
    
    /// Pressure under-relaxation factor (SIMPLE algorithm)
    pub const PRESSURE: f64 = 0.3;
    
    /// Temperature under-relaxation factor
    pub const TEMPERATURE: f64 = 0.8;
    
    /// Turbulence under-relaxation factor
    pub const TURBULENCE: f64 = 0.8;
    
    /// Density under-relaxation factor
    pub const DENSITY: f64 = 0.5;
    
    /// Default relaxation factor
    pub const DEFAULT: f64 = 0.7;
}

/// Discretization scheme parameters
pub mod discretization {
    /// CFL number for explicit schemes
    pub const CFL_EXPLICIT: f64 = 0.5;
    
    /// CFL number for implicit schemes
    pub const CFL_IMPLICIT: f64 = 1.0;
    
    /// Maximum allowable CFL number
    pub const CFL_MAX: f64 = 10.0;
    
    /// Peclet number threshold for upwind switching
    pub const PECLET_THRESHOLD: f64 = 2.0;
    
    /// Maximum cell aspect ratio
    pub const ASPECT_RATIO_MAX: f64 = 10.0;
    
    /// Maximum cell skewness
    pub const SKEWNESS_MAX: f64 = 0.95;
    
    /// Minimum orthogonality quality
    pub const ORTHOGONALITY_MIN: f64 = 0.2;
    
    /// Courant number for VOF advection
    pub const COURANT_VOF: f64 = 0.25;
}

/// Time integration parameters
pub mod time {
    /// Default time step safety factor
    pub const SAFETY_FACTOR: f64 = 0.9;
    
    /// Minimum time step ratio
    pub const MIN_STEP_RATIO: f64 = 0.1;
    
    /// Maximum time step ratio
    pub const MAX_STEP_RATIO: f64 = 2.0;
    
    /// Runge-Kutta 4 coefficients
    pub const RK4_C2: f64 = 0.5;
    pub const RK4_C3: f64 = 0.5;
    pub const RK4_C4: f64 = 1.0;
    pub const RK4_B1: f64 = 1.0 / 6.0;
    pub const RK4_B2: f64 = 1.0 / 3.0;
    pub const RK4_B3: f64 = 1.0 / 3.0;
    pub const RK4_B4: f64 = 1.0 / 6.0;
}

/// Lattice Boltzmann method parameters
pub mod lbm {
    /// D2Q9 lattice weights
    pub const D2Q9_W0: f64 = 4.0 / 9.0;   // Center
    pub const D2Q9_W1: f64 = 1.0 / 9.0;   // Cardinal directions
    pub const D2Q9_W2: f64 = 1.0 / 36.0;  // Diagonal directions
    
    /// D3Q19 lattice weights
    pub const D3Q19_W0: f64 = 1.0 / 3.0;
    pub const D3Q19_W1: f64 = 1.0 / 18.0;
    pub const D3Q19_W2: f64 = 1.0 / 36.0;
    
    /// BGK relaxation time bounds
    pub const TAU_MIN: f64 = 0.5;
    pub const TAU_MAX: f64 = 2.0;
}

/// Common mathematical constants for numerical computations
pub mod math {
    /// Half value constant
    pub const HALF: f64 = 0.5;
    
    /// One and a half
    pub const ONE_AND_HALF: f64 = 1.5;
    
    /// Two
    pub const TWO: f64 = 2.0;
    
    /// Three
    pub const THREE: f64 = 3.0;
    
    /// Four
    pub const FOUR: f64 = 4.0;
    
    /// Quarter value
    pub const QUARTER: f64 = 0.25;
    
    /// Three quarters
    pub const THREE_QUARTERS: f64 = 0.75;
    
    /// One and a quarter
    pub const ONE_AND_QUARTER: f64 = 1.25;
    
    /// Two and a half
    pub const TWO_AND_HALF: f64 = 2.5;
    
    /// Two thirds
    pub const TWO_THIRDS: f64 = 2.0 / 3.0;
    
    /// Pi squared
    pub const PI_SQUARED: f64 = std::f64::consts::PI * std::f64::consts::PI;
    
    /// Two Pi
    pub const TWO_PI: f64 = 2.0 * std::f64::consts::PI;
}