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
}