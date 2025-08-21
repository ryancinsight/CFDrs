//! Constants module for CFD simulations

// Re-export commonly used constants at module level
pub use self::physics::*;

pub mod physics {
    //! Physical constants
    
    /// Critical Reynolds number for transition
    pub const REYNOLDS_CRITICAL: f64 = 2300.0;
    
    /// Von Karman constant
    pub const VON_KARMAN: f64 = 0.41;
    
    /// Default convergence tolerance
    pub const CONVERGENCE_TOLERANCE: f64 = 1e-10;
    
    /// Default maximum iterations
    pub const MAX_ITERATIONS_DEFAULT: usize = 100;
    
    /// Half constant
    pub const HALF: f64 = 0.5;
    
    /// Two constant  
    pub const TWO: f64 = 2.0;
    
    /// One constant
    pub const ONE: f64 = 1.0;
    
    /// Three constant
    pub const THREE: f64 = 3.0;
    
    /// Four constant
    pub const FOUR: f64 = 4.0;
    
    /// One tenth constant
    pub const ONE_TENTH: f64 = 0.1;
    
    /// Laminar flow threshold (Reynolds number)
    pub const LAMINAR_THRESHOLD: f64 = 2300.0;
    
    /// Turbulent flow threshold (Reynolds number)
    pub const TURBULENT_THRESHOLD: f64 = 4000.0;
    
    /// E value for wall functions
    pub const E_WALL_FUNCTION: f64 = 9.8;
    
    /// Default tolerance for iterative solvers
    pub const DEFAULT_TOLERANCE: f64 = 1e-6;
    
    /// Default CFL number for stability
    pub const DEFAULT_CFL_NUMBER: f64 = 0.5;
    
    /// Default time step factor
    pub const DEFAULT_TIME_STEP_FACTOR: f64 = 0.8;
    
    /// Velocity under-relaxation factor
    pub const VELOCITY_UNDER_RELAXATION: f64 = 0.7;
    
    /// Pressure under-relaxation factor
    pub const PRESSURE_UNDER_RELAXATION: f64 = 0.3;
    
    /// Convergence tolerance for velocity
    pub const CONVERGENCE_TOLERANCE_VELOCITY: f64 = 1e-5;
    
    /// Convergence tolerance for pressure
    pub const CONVERGENCE_TOLERANCE_PRESSURE: f64 = 1e-4;
    
    /// Convergence tolerance for continuity equation
    pub const CONVERGENCE_TOLERANCE_CONTINUITY: f64 = 1e-5;
    
    /// Maximum outer iterations
    pub const MAX_ITERATIONS_OUTER: usize = 100;
    
    /// Y+ threshold for laminar sublayer
    pub const Y_PLUS_LAMINAR: f64 = 11.63;
}

pub mod numerical {
    //! Numerical method constants
    
    /// Default convergence tolerance
    pub const DEFAULT_CONVERGENCE_TOLERANCE: f64 = 1e-6;
    
    /// Default maximum iterations
    pub const DEFAULT_MAX_ITERATIONS: usize = 1000;
}