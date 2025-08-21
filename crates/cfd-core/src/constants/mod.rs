//! Constants module for CFD simulations

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
}