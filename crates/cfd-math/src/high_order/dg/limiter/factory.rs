use super::params::LimiterType;
use super::traits::Limiter;
use super::{MinmodLimiter, MomentLimiter, NoLimiter, TVBLimiter, WENOLimiter};

/// Factory for creating limiter instances
pub struct LimiterFactory;

impl LimiterFactory {
    /// Create a new limiter
    pub fn create(limiter_type: LimiterType) -> Box<dyn Limiter> {
        match limiter_type {
            LimiterType::None => Box::new(NoLimiter),
            LimiterType::TVB => Box::new(TVBLimiter),
            LimiterType::Moment => Box::new(MomentLimiter),
            LimiterType::WENO => Box::new(WENOLimiter::new(3)),
            LimiterType::Minmod | LimiterType::MC | LimiterType::Superbee => {
                Box::new(MinmodLimiter)
            }
        }
    }
}
