use super::params::LimiterParams;
use super::super::DGSolution;
use crate::error::Result;

/// Trait for slope limiters
pub trait Limiter: Send + Sync {
    /// Apply the limiter to a DG solution
    fn limit(
        &self,
        solution: &mut DGSolution,
        neighbors: &[DGSolution],
        params: &LimiterParams,
    ) -> Result<()>;

    /// Check if a cell is troubled and needs limiting
    fn is_troubled_cell(
        &self,
        solution: &DGSolution,
        neighbors: &[DGSolution],
        params: &LimiterParams,
    ) -> bool;
}
