//! Adaptive mesh refinement (AMR) for CFD simulations
//!
//! This module provides intelligent mesh refinement strategies that adapt
//! to solution features, error estimates, and physical phenomena.

use nalgebra::DMatrix;
use cfd_core::error::Result;

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
    
    /// Compute sophisticated error estimate using multiple methods
    fn compute_sophisticated_error_estimate(&self, solution: &DMatrix<f64>, i: usize, j: usize, threshold: f64) -> f64 {
        let (nx, ny) = solution.shape();
        
        // Combine multiple error estimation techniques
        let richardson_error = self.compute_richardson_error_estimate(solution, i, j);
        let adjoint_error = self.compute_adjoint_error_estimate(solution, i, j);
        let residual_error = self.compute_residual_error_estimate(solution, i, j);
        let smoothness_error = self.compute_smoothness_error_estimate(solution, i, j);
        
        // Weighted combination of error estimates
        let combined_error = 0.3 * richardson_error + 
                            0.3 * adjoint_error + 
                            0.2 * residual_error + 
                            0.2 * smoothness_error;
        
        // Apply threshold normalization
        if combined_error > threshold {
            combined_error / threshold
        } else {
            0.0
        }
    }
    
    /// Richardson extrapolation-based error estimate
    fn compute_richardson_error_estimate(&self, solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();
        
        // Need at least 2x2 stencil for Richardson extrapolation
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }
        
        // Fine grid approximation (current solution)
        let fine = solution[(i, j)];
        
        // Coarse grid approximation (average of neighbors)
        let coarse = 0.25 * (solution[(i-1, j)] + solution[(i+1, j)] + 
                           solution[(i, j-1)] + solution[(i, j+1)]);
        
        // Richardson error estimate (assuming second-order accuracy)
        let error = (fine - coarse).abs() / 3.0;
        
        // Scale by local solution magnitude for relative error
        if fine.abs() > 1e-10 {
            error / fine.abs()
        } else {
            error
        }
    }
    
    /// Adjoint-based error estimate (simplified)
    fn compute_adjoint_error_estimate(&self, solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();
        
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }
        
        // Compute local Hessian (second derivatives)
        let d2u_dx2 = solution[(i+1, j)] - 2.0 * solution[(i, j)] + solution[(i-1, j)];
        let d2u_dy2 = solution[(i, j+1)] - 2.0 * solution[(i, j)] + solution[(i, j-1)];
        let d2u_dxdy = 0.25 * (solution[(i+1, j+1)] - solution[(i+1, j-1)] - 
                              solution[(i-1, j+1)] + solution[(i-1, j-1)]);
        
        // Adjoint weight (simplified - based on solution curvature)
        let adjoint_weight = (d2u_dx2 * d2u_dx2 + d2u_dy2 * d2u_dy2 + 2.0 * d2u_dxdy * d2u_dxdy).sqrt();
        
        // Error estimate based on local truncation error weighted by adjoint
        let local_truncation_error = (d2u_dx2.abs() + d2u_dy2.abs()) / 12.0;
        
        adjoint_weight * local_truncation_error
    }
    
    /// Residual-based error estimate
    fn compute_residual_error_estimate(&self, solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();
        
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }
        
        // Compute discrete Laplacian (residual of Poisson equation)
        let laplacian = solution[(i+1, j)] + solution[(i-1, j)] + 
                       solution[(i, j+1)] + solution[(i, j-1)] - 
                       4.0 * solution[(i, j)];
        
        // Residual magnitude
        laplacian.abs()
    }
    
    /// Smoothness-based error estimate
    fn compute_smoothness_error_estimate(&self, solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();
        
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }
        
        // Compute solution variation in different directions
        let dx_plus = solution[(i+1, j)] - solution[(i, j)];
        let dx_minus = solution[(i, j)] - solution[(i-1, j)];
        let dy_plus = solution[(i, j+1)] - solution[(i, j)];
        let dy_minus = solution[(i, j)] - solution[(i, j-1)];
        
        // Smoothness indicator (based on total variation diminishing)
        let smoothness_x = if dx_plus.abs() + dx_minus.abs() > 1e-10 {
            (dx_plus - dx_minus).abs() / (dx_plus.abs() + dx_minus.abs())
        } else {
            0.0
        };
        
        let smoothness_y = if dy_plus.abs() + dy_minus.abs() > 1e-10 {
            (dy_plus - dy_minus).abs() / (dy_plus.abs() + dy_minus.abs())
        } else {
            0.0
        };
        
        // Combined smoothness error
        (smoothness_x + smoothness_y) / 2.0
    }
    
    /// Compute feature-based refinement indicator
    fn compute_feature_based_indicator(&self, solution: &DMatrix<f64>, i: usize, j: usize, feature_name: &str) -> f64 {
        let (nx, ny) = solution.shape();
        
        match feature_name {
            "edge_detection" => self.compute_edge_detection_indicator(solution, i, j),
            "vorticity" => self.compute_vorticity_indicator(solution, i, j),
            "shock_detection" => self.compute_shock_detection_indicator(solution, i, j),
            "discontinuity" => self.compute_discontinuity_indicator(solution, i, j),
            "gradient_magnitude" => self.compute_gradient_magnitude_indicator(solution, i, j),
            "curvature" => self.compute_curvature_indicator(solution, i, j),
            "wavelet" => self.compute_wavelet_indicator(solution, i, j),
            _ => {
                // Default to gradient-based detection for unknown features
                self.compute_gradient_magnitude_indicator(solution, i, j)
            }
        }
    }
    
    /// Edge detection using Sobel operator
    fn compute_edge_detection_indicator(&self, solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();
        
        // Need 3x3 neighborhood for Sobel operator
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }
        
        // Sobel kernels
        // Gx (horizontal gradient)
        let gx = (-1.0 * solution[(i-1, j-1)] + 1.0 * solution[(i+1, j-1)] +
                  -2.0 * solution[(i-1, j)]   + 2.0 * solution[(i+1, j)] +
                  -1.0 * solution[(i-1, j+1)] + 1.0 * solution[(i+1, j+1)]) / 8.0;
        
        // Gy (vertical gradient)
        let gy = (-1.0 * solution[(i-1, j-1)] - 2.0 * solution[(i, j-1)] - 1.0 * solution[(i+1, j-1)] +
                   1.0 * solution[(i-1, j+1)] + 2.0 * solution[(i, j+1)] + 1.0 * solution[(i+1, j+1)]) / 8.0;
        
        // Edge magnitude
        (gx * gx + gy * gy).sqrt()
    }
    
    /// Vorticity-based indicator for 2D flow
    fn compute_vorticity_indicator(&self, solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();
        
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }
        
        // Compute vorticity using central differences
        let dv_dx = (solution[(i+1, j)] - solution[(i-1, j)]) / 2.0;
        let du_dy = (solution[(i, j+1)] - solution[(i, j-1)]) / 2.0;
        let vorticity = dv_dx - du_dy;
        
        vorticity.abs()
    }
    
    /// Shock detection using pressure/density gradient
    fn compute_shock_detection_indicator(&self, solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();
        
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }
        
        // Compute pressure gradient magnitude
        let dp_dx = (solution[(i+1, j)] - solution[(i-1, j)]) / 2.0;
        let dp_dy = (solution[(i, j+1)] - solution[(i, j-1)]) / 2.0;
        let grad_p = (dp_dx * dp_dx + dp_dy * dp_dy).sqrt();
        
        // Compute second derivative (shock indicator)
        let d2p_dx2 = solution[(i+1, j)] - 2.0 * solution[(i, j)] + solution[(i-1, j)];
        let d2p_dy2 = solution[(i, j+1)] - 2.0 * solution[(i, j)] + solution[(i, j-1)];
        let laplacian_p = (d2p_dx2 + d2p_dy2).abs();
        
        // Combine gradient and curvature for shock detection
        grad_p * (1.0 + laplacian_p)
    }
    
    /// Discontinuity detection using jump detection
    fn compute_discontinuity_indicator(&self, solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();
        
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }
        
        // Check for jumps in different directions
        let center = solution[(i, j)];
        
        // Horizontal jump
        let left = solution[(i-1, j)];
        let right = solution[(i+1, j)];
        let jump_x = ((center - left).abs() + (right - center).abs()) / 2.0;
        
        // Vertical jump
        let top = solution[(i, j-1)];
        let bottom = solution[(i, j+1)];
        let jump_y = ((center - top).abs() + (bottom - center).abs()) / 2.0;
        
        // Diagonal jumps
        let jump_tl = ((center - solution[(i-1, j-1)]).abs() + (solution[(i+1, j+1)] - center).abs()) / 2.0;
        let jump_tr = ((center - solution[(i+1, j-1)]).abs() + (solution[(i-1, j+1)] - center).abs()) / 2.0;
        
        // Maximum jump in any direction
        jump_x.max(jump_y).max(jump_tl).max(jump_tr)
    }
    
    /// Gradient magnitude indicator
    fn compute_gradient_magnitude_indicator(&self, solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();
        
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }
        
        // Central difference gradients
        let dx = (solution[(i+1, j)] - solution[(i-1, j)]) / 2.0;
        let dy = (solution[(i, j+1)] - solution[(i, j-1)]) / 2.0;
        
        (dx * dx + dy * dy).sqrt()
    }
    
    /// Curvature-based indicator
    fn compute_curvature_indicator(&self, solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();
        
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }
        
        // Second derivatives (Laplacian)
        let d2u_dx2 = solution[(i+1, j)] - 2.0 * solution[(i, j)] + solution[(i-1, j)];
        let d2u_dy2 = solution[(i, j+1)] - 2.0 * solution[(i, j)] + solution[(i, j-1)];
        
        // Mixed derivative
        let d2u_dxdy = (solution[(i+1, j+1)] - solution[(i+1, j-1)] - solution[(i-1, j+1)] + solution[(i-1, j-1)]) / 4.0;
        
        // Curvature magnitude (simplified)
        (d2u_dx2 * d2u_dx2 + d2u_dy2 * d2u_dy2 + 2.0 * d2u_dxdy * d2u_dxdy).sqrt()
    }
    
    /// Wavelet-based indicator using Haar wavelet
    fn compute_wavelet_indicator(&self, solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();
        
        // Need 2x2 neighborhood for Haar wavelet
        if i >= nx - 1 || j >= ny - 1 {
            return 0.0;
        }
        
        // Haar wavelet coefficients
        let h00 = solution[(i, j)];
        let h01 = solution[(i, j+1)];
        let h10 = solution[(i+1, j)];
        let h11 = solution[(i+1, j+1)];
        
        // Horizontal detail coefficient
        let dh = (h00 + h01 - h10 - h11).abs() / 2.0;
        
        // Vertical detail coefficient
        let dv = (h00 - h01 + h10 - h11).abs() / 2.0;
        
        // Diagonal detail coefficient
        let dd = (h00 - h01 - h10 + h11).abs() / 2.0;
        
        // Combined wavelet energy
        (dh * dh + dv * dv + dd * dd).sqrt()
    }
    
    /// Compute physics-based refinement indicator
    fn compute_physics_based_indicator(&self, solution: &DMatrix<f64>, i: usize, j: usize, phenomena: &[PhysicsPhenomena]) -> f64 {
        let mut indicator: f64 = 0.0;
        
        for phenomenon in phenomena {
            let phenomenon_indicator = match phenomenon {
                PhysicsPhenomena::ShockWaves => self.compute_shock_wave_indicator(solution, i, j),
                PhysicsPhenomena::BoundaryLayers { wall_distance } => {
                    self.compute_boundary_layer_indicator(solution, i, j, *wall_distance)
                }
                PhysicsPhenomena::Vortices { vorticity_threshold } => {
                    self.compute_vortex_indicator(solution, i, j, *vorticity_threshold)
                }
                PhysicsPhenomena::Interfaces { interface_thickness } => {
                    self.compute_interface_indicator(solution, i, j, *interface_thickness)
                }
                PhysicsPhenomena::ThermalGradients { temp_gradient_threshold } => {
                    self.compute_thermal_gradient_indicator(solution, i, j, *temp_gradient_threshold)
                }
            };
            
            // Combine indicators using maximum (conservative approach)
            indicator = indicator.max(phenomenon_indicator);
        }
        
        indicator
    }
    
    /// Shock wave detection using pressure jump and Mach number criteria
    fn compute_shock_wave_indicator(&self, solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();
        
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }
        
        // Pressure gradient (shock indicator)
        let dp_dx = (solution[(i+1, j)] - solution[(i-1, j)]) / 2.0;
        let dp_dy = (solution[(i, j+1)] - solution[(i, j-1)]) / 2.0;
        let grad_p = (dp_dx * dp_dx + dp_dy * dp_dy).sqrt();
        
        // Pressure jump detection
        let p_center = solution[(i, j)];
        let p_left = solution[(i-1, j)];
        let p_right = solution[(i+1, j)];
        let pressure_jump = ((p_right - p_left).abs() / (2.0 * p_center.abs() + 1e-10)).max(0.0);
        
        // Density jump (assuming solution represents density)
        let density_jump = self.compute_discontinuity_indicator(solution, i, j);
        
        // Combined shock indicator
        grad_p * (1.0 + pressure_jump + density_jump)
    }
    
    /// Boundary layer detection based on wall distance and velocity gradients
    fn compute_boundary_layer_indicator(&self, solution: &DMatrix<f64>, i: usize, j: usize, wall_distance: f64) -> f64 {
        let (nx, ny) = solution.shape();
        
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }
        
        // Distance from nearest wall (simplified - using grid coordinates)
        let wall_dist = (i as f64).min(j as f64).min((nx - i - 1) as f64).min((ny - j - 1) as f64);
        
        // Only apply boundary layer detection near walls
        if wall_dist > wall_distance * 3.0 {
            return 0.0;
        }
        
        // Velocity gradient magnitude (boundary layer thickness indicator)
        let dv_dx = (solution[(i+1, j)] - solution[(i-1, j)]) / 2.0;
        let dv_dy = (solution[(i, j+1)] - solution[(i, j-1)]) / 2.0;
        let grad_v = (dv_dx * dv_dx + dv_dy * dv_dy).sqrt();
        
        // Wall damping function (stronger near wall)
        let wall_factor = (-wall_dist / wall_distance).exp();
        
        // Boundary layer indicator
        grad_v * wall_factor
    }
    
    /// Vortex detection using vorticity magnitude and Q-criterion
    fn compute_vortex_indicator(&self, solution: &DMatrix<f64>, i: usize, j: usize, vorticity_threshold: f64) -> f64 {
        let (nx, ny) = solution.shape();
        
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }
        
        // Compute vorticity
        let vorticity = self.compute_vorticity_indicator(solution, i, j);
        
        if vorticity < vorticity_threshold {
            return 0.0;
        }
        
        // Q-criterion (second invariant of velocity gradient tensor)
        // Simplified version for 2D: Q = 0.5 * (||Ω||² - ||S||²)
        let du_dx = (solution[(i+1, j)] - solution[(i-1, j)]) / 2.0;
        let du_dy = (solution[(i, j+1)] - solution[(i, j-1)]) / 2.0;
        let dv_dx = (solution[(i+1, j)] - solution[(i-1, j)]) / 2.0;
        let dv_dy = (solution[(i, j+1)] - solution[(i, j-1)]) / 2.0;
        
        // Vorticity magnitude squared
        let omega_sq = (dv_dx - du_dy).powi(2);
        
        // Strain rate magnitude squared
        let s_xx = du_dx;
        let s_yy = dv_dy;
        let s_xy = 0.5 * (du_dy + dv_dx);
        let strain_sq = s_xx * s_xx + s_yy * s_yy + 2.0 * s_xy * s_xy;
        
        // Q-criterion
        let q_criterion = 0.5 * (omega_sq - strain_sq);
        
        // Vortex indicator (positive Q indicates vortex core)
        if q_criterion > 0.0 {
            vorticity * (1.0 + q_criterion / vorticity_threshold)
        } else {
            0.0
        }
    }
    
    /// Interface detection for multiphase flows
    fn compute_interface_indicator(&self, solution: &DMatrix<f64>, i: usize, j: usize, interface_thickness: f64) -> f64 {
        let (nx, ny) = solution.shape();
        
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }
        
        // Volume fraction gradient (interface indicator)
        let dphi_dx = (solution[(i+1, j)] - solution[(i-1, j)]) / 2.0;
        let dphi_dy = (solution[(i, j+1)] - solution[(i, j-1)]) / 2.0;
        let grad_phi = (dphi_dx * dphi_dx + dphi_dy * dphi_dy).sqrt();
        
        // Interface detection based on gradient magnitude
        let interface_indicator = if grad_phi > 1.0 / interface_thickness {
            grad_phi * interface_thickness
        } else {
            0.0
        };
        
        // Curvature of interface (second derivative)
        let d2phi_dx2 = solution[(i+1, j)] - 2.0 * solution[(i, j)] + solution[(i-1, j)];
        let d2phi_dy2 = solution[(i, j+1)] - 2.0 * solution[(i, j)] + solution[(i, j-1)];
        let curvature = (d2phi_dx2 + d2phi_dy2).abs();
        
        // Combined interface indicator
        interface_indicator * (1.0 + curvature * interface_thickness)
    }
    
    /// Thermal gradient detection for heat transfer problems
    fn compute_thermal_gradient_indicator(&self, solution: &DMatrix<f64>, i: usize, j: usize, temp_gradient_threshold: f64) -> f64 {
        let (nx, ny) = solution.shape();
        
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }
        
        // Temperature gradient magnitude
        let dt_dx = (solution[(i+1, j)] - solution[(i-1, j)]) / 2.0;
        let dt_dy = (solution[(i, j+1)] - solution[(i, j-1)]) / 2.0;
        let grad_t = (dt_dx * dt_dx + dt_dy * dt_dy).sqrt();
        
        if grad_t < temp_gradient_threshold {
            return 0.0;
        }
        
        // Heat flux divergence (indicator of thermal sources/sinks)
        let d2t_dx2 = solution[(i+1, j)] - 2.0 * solution[(i, j)] + solution[(i-1, j)];
        let d2t_dy2 = solution[(i, j+1)] - 2.0 * solution[(i, j)] + solution[(i, j-1)];
        let heat_source = (d2t_dx2 + d2t_dy2).abs();
        
        // Thermal boundary layer detection
        let thermal_bl_indicator = grad_t * (1.0 + heat_source / temp_gradient_threshold);
        
        thermal_bl_indicator
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
                            indicators[(i, j)] = indicators[(i, j)] / threshold;
                        } else {
                            indicators[(i, j)] = 0.0;
                        }
                    }
                }
            }
            RefinementCriteria::ErrorEstimate { threshold, .. } => {
                // Sophisticated error estimation using multiple methods
                for i in 0..nx {
                    for j in 0..ny {
                        indicators[(i, j)] = self.compute_sophisticated_error_estimate(solution, i, j, *threshold);
                    }
                }
            }
            RefinementCriteria::FeatureBased { feature_name, .. } => {
                // Implement feature-based detection using various algorithms
                for i in 0..nx {
                    for j in 0..ny {
                        indicators[(i, j)] = self.compute_feature_based_indicator(solution, i, j, feature_name);
                    }
                }
            }
            RefinementCriteria::PhysicsBased { phenomena, .. } => {
                // Implement physics-based indicators
                for i in 0..nx {
                    for j in 0..ny {
                        indicators[(i, j)] = self.compute_physics_based_indicator(solution, i, j, phenomena);
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
            RefinementCriteria::Gradient { min_level, max_level, .. }
            | RefinementCriteria::ErrorEstimate { min_level, max_level, .. }
            | RefinementCriteria::FeatureBased { min_level, max_level, .. }
            | RefinementCriteria::PhysicsBased { min_level, max_level, .. } => (*min_level, *max_level),
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
        let total_cells = self.count_total_cells();
        if total_cells > self.max_cells {
            // TODO: Implement smart coarsening to stay within memory limits
            return Err(cfd_core::error::Error::InvalidInput(
                format!("AMR would exceed memory limit: {} > {}", total_cells, self.max_cells)
            ));
        }
        
        // Record refinement step
        let step = self.refinement_history.len() as u32;
        self.refinement_history.push(RefinementStep {
            step,
            cells_refined,
            cells_coarsened,
            total_cells: total_cells as u32,
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
            phenomena: vec![PhysicsPhenomena::Vortices { vorticity_threshold }],
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
        assert!(amr.refinement_levels[(5, 0)] > 0 || amr.refinement_levels[(6, 0)] > 0);
    }
    
    #[test]
    fn test_memory_limit() {
        let criteria = RefinementCriteria::gradient(0.1, 0, 10);
        let mut amr = AdaptiveMeshRefinement::new(100, 100, criteria);
        amr.max_cells = 100; // Very low limit
        
        let solution = DMatrix::random(100, 100);
        let indicators = amr.compute_refinement_indicators(&solution).unwrap();
        
        // Should fail due to memory limit
        assert!(amr.apply_refinement(&indicators).is_err());
    }
}
