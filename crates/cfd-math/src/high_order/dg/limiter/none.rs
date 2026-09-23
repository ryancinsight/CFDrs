use super::traits::Limiter;
use super::params::LimiterParams;
use super::super::DGSolution;
use crate::error::Result;

/// No limiting (identity operator)
pub struct NoLimiter;

impl Limiter for NoLimiter {
    fn limit(
        &self,
        _solution: &mut DGSolution,
        _neighbors: &[DGSolution],
        _params: &LimiterParams,
    ) -> Result<()> {
        // Do nothing
        Ok(())
    }

    fn is_troubled_cell(
        &self,
        _solution: &DGSolution,
        _neighbors: &[DGSolution],
        _params: &LimiterParams,
    ) -> bool {
        false
    }
}
