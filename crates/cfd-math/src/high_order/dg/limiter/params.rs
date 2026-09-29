/// Type of slope limiter
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LimiterType {
    /// No limiting (identity operator)
    None,
    /// Minmod limiter (most diffusive)
    Minmod,
    /// TVB (Total Variation Bounded) limiter
    TVB,
    /// Moment limiter (for high-order methods)
    Moment,
    /// WENO (Weighted Essentially Non-Oscillatory) limiter
    WENO,
}

/// Parameters for slope limiters
#[derive(Debug, Clone)]
pub struct LimiterParams {
    /// Type of limiter
    pub limiter_type: LimiterType,
    /// TVB parameter (for TVB limiter)
    pub tvb_m: f64,
    /// WENO parameters (epsilon and p)
    pub weno_epsilon: f64,
    /// Power parameter for WENO weights
    pub weno_p: f64,
    /// Whether to apply the limiter adaptively
    pub adaptive: bool,
    /// Tolerance for detecting troubled cells
    pub tolerance: f64,
}

impl Default for LimiterParams {
    fn default() -> Self {
        Self {
            limiter_type: LimiterType::Minmod,
            tvb_m: 1.0,
            weno_epsilon: 1e-6,
            weno_p: 2.0,
            adaptive: true,
            tolerance: 1e-4,
        }
    }
}

impl LimiterParams {
    /// Create a new set of limiter parameters
    pub fn new(limiter_type: LimiterType) -> Self {
        Self {
            limiter_type,
            ..Default::default()
        }
    }

    /// Set the TVB parameter
    pub fn with_tvb_m(mut self, m: f64) -> Self {
        self.tvb_m = m;
        self
    }

    /// Set the WENO parameters
    pub fn with_weno_params(mut self, epsilon: f64, p: f64) -> Self {
        self.weno_epsilon = epsilon;
        self.weno_p = p;
        self
    }

    /// Set the adaptive flag
    pub fn with_adaptive(mut self, adaptive: bool) -> Self {
        self.adaptive = adaptive;
        self
    }

    /// Set the tolerance for detecting troubled cells
    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance;
        self
    }
}
