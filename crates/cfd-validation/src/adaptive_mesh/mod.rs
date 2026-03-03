//! Adaptive mesh refinement (AMR) for CFD simulations
//!
//! This module provides intelligent mesh refinement strategies that adapt
//! to solution features, error estimates, and physical phenomena.
//!
//! # Submodules
//!
//! | Module              | Purpose                                             |
//! |---------------------|-----------------------------------------------------|
//! | `error_estimators`  | Richardson, adjoint, residual, smoothness estimators |
//! | `indicators`        | Feature-based and physics-based refinement indicators|

mod error_estimators;
mod indicators;

use cfd_core::error::Result;
use nalgebra::DMatrix;

/// Refinement criteria for adaptive mesh refinement
#[derive(Debug, Clone, PartialEq)]
pub enum RefinementCriteria {
    /// Gradient-based refinement (refine where solution gradient is high)
    Gradient {
        /// Gradient threshold for refinement
        threshold: f64,
        /// Minimum refinement level
        min_level: u32,
        /// Maximum refinement level
        max_level: u32,
    },
    /// Error-based refinement (refine where error estimate is high)
    ErrorEstimate {
        /// Error threshold for refinement
        threshold: f64,
        /// Minimum refinement level
        min_level: u32,
        /// Maximum refinement level
        max_level: u32,
    },
    /// Feature-based refinement (refine near specific features)
    FeatureBased {
        /// Feature detection function name (for serialization)
        feature_name: String,
        /// Minimum refinement level
        min_level: u32,
        /// Maximum refinement level
        max_level: u32,
    },
    /// Physics-based refinement (refine based on physical phenomena)
    PhysicsBased {
        /// Physical phenomena to track
        phenomena: Vec<PhysicsPhenomena>,
        /// Minimum refinement level
        min_level: u32,
        /// Maximum refinement level
        max_level: u32,
    },
}

/// Physical phenomena that require refined mesh
#[derive(Debug, Clone, PartialEq)]
pub enum PhysicsPhenomena {
    /// Shock waves and discontinuities
    ShockWaves,
    /// Boundary layers
    BoundaryLayers {
        /// Wall distance threshold
        wall_distance: f64,
    },
    /// Vortices and turbulent structures
    Vortices {
        /// Vorticity threshold
        vorticity_threshold: f64,
    },
    /// Multiphase interfaces
    Interfaces {
        /// Interface thickness
        interface_thickness: f64,
    },
    /// Thermal gradients
    ThermalGradients {
        /// Temperature gradient threshold
        temp_gradient_threshold: f64,
    },
}

/// Adaptive mesh refinement manager
pub struct AdaptiveMeshRefinement {
    /// Current refinement level for each cell
    refinement_levels: DMatrix<u32>,
    /// Refinement criteria
    criteria: RefinementCriteria,
    /// Maximum cells allowed (for memory management)
    max_cells: usize,
    /// Refinement history
    refinement_history: Vec<RefinementStep>,
}

/// Record of a refinement step
#[derive(Debug, Clone)]
pub struct RefinementStep {
    /// Step number
    pub step: u32,
    /// Number of cells refined
    pub cells_refined: u32,
    /// Number of cells coarsened
    pub cells_coarsened: u32,
    /// Total cells after refinement
    pub total_cells: u32,
    /// Refinement criterion used
    pub criterion: String,
    /// Maximum refinement level
    pub max_level: u32,
}

impl AdaptiveMeshRefinement {
    /// Create new AMR manager
    pub fn new(nx: usize, ny: usize, criteria: RefinementCriteria) -> Self {
        Self {
            refinement_levels: DMatrix::zeros(nx, ny),
            criteria,
            max_cells: nx * ny * 8, // Allow up to 8x refinement
            refinement_history: Vec::new(),
        }
    }

    /// Compute refinement indicators based on solution
    pub fn compute_refinement_indicators(&self, solution: &DMatrix<f64>) -> Result<DMatrix<f64>> {
        let nx = solution.nrows();
        let ny = solution.ncols();
        let mut indicators = DMatrix::zeros(nx, ny);

        match &self.criteria {
            RefinementCriteria::Gradient { threshold, .. } => {
                // Compute gradient magnitude using central differences
                for i in 1..nx - 1 {
                    for j in 1..ny - 1 {
                        let dx = (solution[(i + 1, j)] - solution[(i - 1, j)]) / 2.0;
                        let dy = (solution[(i, j + 1)] - solution[(i, j - 1)]) / 2.0;
                        indicators[(i, j)] = (dx * dx + dy * dy).sqrt();

                        // Normalize by threshold
                        if indicators[(i, j)] > *threshold {
                            indicators[(i, j)] /= *threshold;
                        } else {
                            indicators[(i, j)] = 0.0;
                        }
                    }
                }
            }
            RefinementCriteria::ErrorEstimate { threshold, .. } => {
                for i in 0..nx {
                    for j in 0..ny {
                        indicators[(i, j)] =
                            Self::compute_sophisticated_error_estimate(solution, i, j, *threshold);
                    }
                }
            }
            RefinementCriteria::FeatureBased { feature_name, .. } => {
                for i in 0..nx {
                    for j in 0..ny {
                        indicators[(i, j)] =
                            Self::compute_feature_based_indicator(solution, i, j, feature_name);
                    }
                }
            }
            RefinementCriteria::PhysicsBased { phenomena, .. } => {
                for i in 0..nx {
                    for j in 0..ny {
                        indicators[(i, j)] =
                            Self::compute_physics_based_indicator(solution, i, j, phenomena);
                    }
                }
            }
        }

        Ok(indicators)
    }

    /// Apply refinement based on indicators
    pub fn apply_refinement(&mut self, indicators: &DMatrix<f64>) -> Result<()> {
        let (nx, ny) = indicators.shape();
        let mut cells_refined = 0;
        let mut cells_coarsened = 0;

        let (min_level, max_level) = match &self.criteria {
            RefinementCriteria::Gradient {
                min_level,
                max_level,
                ..
            }
            | RefinementCriteria::ErrorEstimate {
                min_level,
                max_level,
                ..
            }
            | RefinementCriteria::FeatureBased {
                min_level,
                max_level,
                ..
            }
            | RefinementCriteria::PhysicsBased {
                min_level,
                max_level,
                ..
            } => (*min_level, *max_level),
        };

        for i in 0..nx {
            for j in 0..ny {
                let current_level = self.refinement_levels[(i, j)];
                let indicator = indicators[(i, j)];

                if indicator > 0.5 && current_level < max_level {
                    // Refine this cell
                    self.refinement_levels[(i, j)] += 1;
                    cells_refined += 1;
                } else if indicator < 0.1 && current_level > min_level {
                    // Coarsen this cell
                    self.refinement_levels[(i, j)] -= 1;
                    cells_coarsened += 1;
                }
            }
        }

        // Check memory constraint
        let mut total_cells = self.count_total_cells();
        if total_cells > self.max_cells {
            let mut candidates: Vec<(f64, u32, usize, usize)> = Vec::new();
            for i in 0..nx {
                for j in 0..ny {
                    let level = self.refinement_levels[(i, j)];
                    if level > min_level {
                        candidates.push((indicators[(i, j)], level, i, j));
                    }
                }
            }

            candidates.sort_by(|a, b| {
                a.0.partial_cmp(&b.0)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| b.1.cmp(&a.1))
            });

            for &(_, level, i, j) in &candidates {
                if total_cells <= self.max_cells {
                    break;
                }
                if level <= min_level {
                    continue;
                }

                let current_cells = 4_usize.pow(level);
                let reduced_cells = 4_usize.pow(level - 1);
                let reduction = current_cells.saturating_sub(reduced_cells);
                if reduction == 0 {
                    continue;
                }

                self.refinement_levels[(i, j)] -= 1;
                total_cells = total_cells.saturating_sub(reduction);
                cells_coarsened += 1;
            }

            if total_cells > self.max_cells {
                return Err(cfd_core::error::Error::InvalidInput(format!(
                    "AMR would exceed memory limit: {} > {}",
                    total_cells, self.max_cells
                )));
            }
        }

        // Record refinement step
        let step = u32::try_from(self.refinement_history.len()).map_err(|_| {
            cfd_core::error::Error::InvalidInput("AMR refinement step exceeds u32::MAX".to_string())
        })?;
        let total_cells_u32 = u32::try_from(total_cells).map_err(|_| {
            cfd_core::error::Error::InvalidInput(
                "AMR total cell count exceeds u32::MAX".to_string(),
            )
        })?;
        self.refinement_history.push(RefinementStep {
            step,
            cells_refined,
            cells_coarsened,
            total_cells: total_cells_u32,
            criterion: format!("{:?}", self.criteria),
            max_level,
        });

        Ok(())
    }

    /// Count total cells considering refinement levels
    fn count_total_cells(&self) -> usize {
        let (nx, ny) = self.refinement_levels.shape();
        let mut total = 0;

        for i in 0..nx {
            for j in 0..ny {
                let level = self.refinement_levels[(i, j)];
                // Each refinement level doubles resolution in each direction
                let cells_per_level = 2_usize.pow(level) * 2_usize.pow(level);
                total += cells_per_level;
            }
        }

        total
    }

    /// Get current refinement levels
    pub fn refinement_levels(&self) -> &DMatrix<u32> {
        &self.refinement_levels
    }

    /// Get refinement history
    pub fn refinement_history(&self) -> &[RefinementStep] {
        &self.refinement_history
    }

    /// Reset refinement levels
    pub fn reset(&mut self) {
        self.refinement_levels.fill(0);
        self.refinement_history.clear();
    }

    /// Estimate memory usage of refined mesh
    pub fn estimate_memory_usage(&self) -> usize {
        let total_cells = self.count_total_cells();
        // Assume 8 bytes per value (f64) and 5 values per cell (typical CFD variables)
        total_cells * 8 * 5
    }
}

/// Factory for common refinement criteria
impl RefinementCriteria {
    /// Create gradient-based refinement criterion
    pub fn gradient(threshold: f64, min_level: u32, max_level: u32) -> Self {
        Self::Gradient {
            threshold,
            min_level,
            max_level,
        }
    }

    /// Create error-based refinement criterion
    pub fn error_estimate(threshold: f64, min_level: u32, max_level: u32) -> Self {
        Self::ErrorEstimate {
            threshold,
            min_level,
            max_level,
        }
    }

    /// Create shock wave refinement criterion
    pub fn shock_waves(min_level: u32, max_level: u32) -> Self {
        Self::PhysicsBased {
            phenomena: vec![PhysicsPhenomena::ShockWaves],
            min_level,
            max_level,
        }
    }

    /// Create boundary layer refinement criterion
    pub fn boundary_layers(wall_distance: f64, min_level: u32, max_level: u32) -> Self {
        Self::PhysicsBased {
            phenomena: vec![PhysicsPhenomena::BoundaryLayers { wall_distance }],
            min_level,
            max_level,
        }
    }

    /// Create vortex refinement criterion
    pub fn vortices(vorticity_threshold: f64, min_level: u32, max_level: u32) -> Self {
        Self::PhysicsBased {
            phenomena: vec![PhysicsPhenomena::Vortices {
                vorticity_threshold,
            }],
            min_level,
            max_level,
        }
    }

    /// Create feature-based refinement criterion
    pub fn feature_based(feature_name: &str, min_level: u32, max_level: u32) -> Self {
        Self::FeatureBased {
            feature_name: feature_name.to_string(),
            min_level,
            max_level,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gradient_refinement() {
        let criteria = RefinementCriteria::gradient(0.1, 0, 3);
        let mut amr = AdaptiveMeshRefinement::new(10, 10, criteria);

        // Create a solution with high gradient region
        let mut solution = DMatrix::zeros(10, 10);
        for i in 0..10 {
            for j in 0..10 {
                solution[(i, j)] = if i > 5 { 1.0 } else { 0.0 };
            }
        }

        let indicators = amr.compute_refinement_indicators(&solution).unwrap();
        amr.apply_refinement(&indicators).unwrap();

        // Check that refinement was applied near the gradient
        // Note: we check column 1 because central differences skip the boundary (j=0)
        assert!(amr.refinement_levels[(5, 1)] > 0 || amr.refinement_levels[(6, 1)] > 0);
    }

    #[test]
    fn test_memory_limit() {
        let criteria = RefinementCriteria::gradient(0.1, 0, 10);
        let mut amr = AdaptiveMeshRefinement::new(100, 100, criteria);
        amr.max_cells = 100; // Very low limit

        let solution = DMatrix::from_fn(100, 100, |i, j| {
            (i as f64 * 0.1).sin() * (j as f64 * 0.1).cos()
        });
        let indicators = amr.compute_refinement_indicators(&solution).unwrap();

        // Should fail due to memory limit
        assert!(amr.apply_refinement(&indicators).is_err());
    }
}
