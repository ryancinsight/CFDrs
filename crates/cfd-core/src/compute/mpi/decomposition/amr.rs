//! Adaptive mesh refinement support for distributed grids.

use super::domain::DomainDecomposition;
use super::load_balancer::LoadBalancer;
use crate::compute::mpi::error::MpiResult;

/// Adaptive mesh refinement support
pub struct AdaptiveMeshRefinement {
    /// Current refinement level
    refinement_level: usize,
    /// Maximum allowed refinement level
    max_refinement_level: usize,
    /// Refinement criteria (error thresholds)
    refinement_criteria: RefinementCriteria,
    /// Load balancer for refined meshes
    load_balancer: Option<LoadBalancer>,
}

impl AdaptiveMeshRefinement {
    /// Create new AMR system
    pub fn new(
        max_refinement_level: usize,
        refinement_criteria: RefinementCriteria,
        load_balancer: Option<LoadBalancer>,
    ) -> Self {
        Self {
            refinement_level: 0,
            max_refinement_level,
            refinement_criteria,
            load_balancer,
        }
    }

    /// Check if cells need refinement based on error estimates
    pub fn needs_refinement(&self, error_estimates: &[f64]) -> Vec<bool> {
        error_estimates
            .iter()
            .map(|&error| {
                error > self.refinement_criteria.error_threshold
                    && self.refinement_level < self.max_refinement_level
            })
            .collect()
    }

    /// Check if cells need coarsening
    pub fn needs_coarsening(&self, error_estimates: &[f64]) -> Vec<bool> {
        error_estimates
            .iter()
            .map(|&error| {
                error < self.refinement_criteria.coarsening_threshold
                    && self.refinement_level > 0
            })
            .collect()
    }

    /// Perform mesh adaptation
    pub fn adapt_mesh(
        &mut self,
        error_estimates: &[f64],
        current_decomp: &DomainDecomposition,
    ) -> MpiResult<DomainDecomposition> {
        let refine_flags = self.needs_refinement(error_estimates);
        let coarsen_flags = self.needs_coarsening(error_estimates);

        // Count cells needing adaptation
        let refine_count = refine_flags.iter().filter(|&&x| x).count();
        let coarsen_count = coarsen_flags.iter().filter(|&&x| x).count();

        // If significant adaptation needed, trigger load rebalancing
        if refine_count > error_estimates.len() / 10
            || coarsen_count > error_estimates.len() / 10
        {
            if let Some(ref mut balancer) = self.load_balancer {
                // Create workload estimate based on refinement
                let local_workload = error_estimates.len() + refine_count - coarsen_count;
                let metrics = balancer.assess_load_balance(local_workload)?;

                if balancer.should_repartition(&metrics) {
                    let new_workloads =
                        vec![local_workload; current_decomp.communicator.size() as usize];
                    return balancer.repartition(&new_workloads);
                }
            }
        }

        // Return current decomposition if no rebalancing needed
        Ok(current_decomp.clone())
    }

    /// Get current refinement level
    pub fn refinement_level(&self) -> usize {
        self.refinement_level
    }

    /// Increment refinement level
    pub fn increment_level(&mut self) {
        if self.refinement_level < self.max_refinement_level {
            self.refinement_level += 1;
        }
    }
}

/// Refinement criteria for adaptive mesh refinement
#[derive(Debug, Clone)]
pub struct RefinementCriteria {
    /// Error threshold above which cells are refined
    pub error_threshold: f64,
    /// Error threshold below which cells are coarsened
    pub coarsening_threshold: f64,
    /// Maximum refinement ratio between adjacent cells
    pub max_refinement_ratio: usize,
}
