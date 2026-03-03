//! Feature-based and physics-based refinement indicators.
//!
//! Implements Sobel edge detection, vorticity, shock detection, wavelet
//! decomposition, Q-criterion, boundary-layer proximity, and multiphase
//! interface indicators for adaptive mesh refinement.

use nalgebra::DMatrix;

use super::{AdaptiveMeshRefinement, PhysicsPhenomena};

// ── Feature-based indicators ─────────────────────────────────────────────────

impl AdaptiveMeshRefinement {
    /// Compute feature-based refinement indicator
    pub(super) fn compute_feature_based_indicator(
        solution: &DMatrix<f64>,
        i: usize,
        j: usize,
        feature_name: &str,
    ) -> f64 {
        match feature_name {
            "edge_detection" => Self::compute_edge_detection_indicator(solution, i, j),
            "vorticity" => Self::compute_vorticity_indicator(solution, i, j),
            "shock_detection" => Self::compute_shock_detection_indicator(solution, i, j),
            "discontinuity" => Self::compute_discontinuity_indicator(solution, i, j),
            "gradient_magnitude" => Self::compute_gradient_magnitude_indicator(solution, i, j),
            "curvature" => Self::compute_curvature_indicator(solution, i, j),
            "wavelet" => Self::compute_wavelet_indicator(solution, i, j),
            _ => {
                // Default to gradient-based detection for unknown features
                Self::compute_gradient_magnitude_indicator(solution, i, j)
            }
        }
    }

    /// Edge detection using Sobel operator
    fn compute_edge_detection_indicator(solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();

        // Need 3x3 neighborhood for Sobel operator
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // Sobel kernels
        // Gx (horizontal gradient)
        let gx = (-solution[(i - 1, j - 1)] + solution[(i + 1, j - 1)]
            - 2.0 * solution[(i - 1, j)]
            + 2.0 * solution[(i + 1, j)]
            - solution[(i - 1, j + 1)]
            + solution[(i + 1, j + 1)])
            / 8.0;

        // Gy (vertical gradient)
        let gy = (-solution[(i - 1, j - 1)]
            - 2.0 * solution[(i, j - 1)]
            - 1.0 * solution[(i + 1, j - 1)]
            + 1.0 * solution[(i - 1, j + 1)]
            + 2.0 * solution[(i, j + 1)]
            + 1.0 * solution[(i + 1, j + 1)])
            / 8.0;

        // Edge magnitude
        (gx * gx + gy * gy).sqrt()
    }

    /// Vorticity-based indicator for 2D flow
    pub(super) fn compute_vorticity_indicator(
        solution: &DMatrix<f64>,
        i: usize,
        j: usize,
    ) -> f64 {
        let (nx, ny) = solution.shape();

        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // Compute vorticity using central differences
        let dv_dx = (solution[(i + 1, j)] - solution[(i - 1, j)]) / 2.0;
        let du_dy = (solution[(i, j + 1)] - solution[(i, j - 1)]) / 2.0;
        let vorticity = dv_dx - du_dy;

        vorticity.abs()
    }

    /// Shock detection using pressure/density gradient
    fn compute_shock_detection_indicator(solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();

        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // Compute pressure gradient magnitude
        let dp_dx = (solution[(i + 1, j)] - solution[(i - 1, j)]) / 2.0;
        let dp_dy = (solution[(i, j + 1)] - solution[(i, j - 1)]) / 2.0;
        let grad_p = (dp_dx * dp_dx + dp_dy * dp_dy).sqrt();

        // Compute second derivative (shock indicator)
        let d2p_dx2 = solution[(i + 1, j)] - 2.0 * solution[(i, j)] + solution[(i - 1, j)];
        let d2p_dy2 = solution[(i, j + 1)] - 2.0 * solution[(i, j)] + solution[(i, j - 1)];
        let laplacian_p = (d2p_dx2 + d2p_dy2).abs();

        // Combine gradient and curvature for shock detection
        grad_p * (1.0 + laplacian_p)
    }

    /// Discontinuity detection using jump detection
    pub(super) fn compute_discontinuity_indicator(
        solution: &DMatrix<f64>,
        i: usize,
        j: usize,
    ) -> f64 {
        let (nx, ny) = solution.shape();

        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // Check for jumps in different directions
        let center = solution[(i, j)];

        // Horizontal jump
        let left = solution[(i - 1, j)];
        let right = solution[(i + 1, j)];
        let jump_x = f64::midpoint((center - left).abs(), (right - center).abs());

        // Vertical jump
        let top = solution[(i, j - 1)];
        let bottom = solution[(i, j + 1)];
        let jump_y = f64::midpoint((center - top).abs(), (bottom - center).abs());

        // Diagonal jumps
        let jump_tl = f64::midpoint(
            (center - solution[(i - 1, j - 1)]).abs(),
            (solution[(i + 1, j + 1)] - center).abs(),
        );
        let jump_tr = f64::midpoint(
            (center - solution[(i + 1, j - 1)]).abs(),
            (solution[(i - 1, j + 1)] - center).abs(),
        );

        // Maximum jump in any direction
        jump_x.max(jump_y).max(jump_tl).max(jump_tr)
    }

    /// Gradient magnitude indicator
    fn compute_gradient_magnitude_indicator(solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();

        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // Central difference gradients
        let dx = (solution[(i + 1, j)] - solution[(i - 1, j)]) / 2.0;
        let dy = (solution[(i, j + 1)] - solution[(i, j - 1)]) / 2.0;

        (dx * dx + dy * dy).sqrt()
    }

    /// Curvature-based indicator
    fn compute_curvature_indicator(solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();

        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // Second derivatives (Laplacian)
        let d2u_dx2 = solution[(i + 1, j)] - 2.0 * solution[(i, j)] + solution[(i - 1, j)];
        let d2u_dy2 = solution[(i, j + 1)] - 2.0 * solution[(i, j)] + solution[(i, j - 1)];

        // Mixed derivative
        let d2u_dxdy =
            (solution[(i + 1, j + 1)] - solution[(i + 1, j - 1)] - solution[(i - 1, j + 1)]
                + solution[(i - 1, j - 1)])
                / 4.0;

        // Curvature magnitude (simplified)
        (d2u_dx2 * d2u_dx2 + d2u_dy2 * d2u_dy2 + 2.0 * d2u_dxdy * d2u_dxdy).sqrt()
    }

    /// Wavelet-based indicator using Haar wavelet
    fn compute_wavelet_indicator(solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();

        // Need 2x2 neighborhood for Haar wavelet
        if i >= nx - 1 || j >= ny - 1 {
            return 0.0;
        }

        // Haar wavelet coefficients
        let h00 = solution[(i, j)];
        let h01 = solution[(i, j + 1)];
        let h10 = solution[(i + 1, j)];
        let h11 = solution[(i + 1, j + 1)];

        // Horizontal detail coefficient
        let dh = (h00 + h01 - h10 - h11).abs() / 2.0;

        // Vertical detail coefficient
        let dv = (h00 - h01 + h10 - h11).abs() / 2.0;

        // Diagonal detail coefficient
        let dd = (h00 - h01 - h10 + h11).abs() / 2.0;

        // Combined wavelet energy
        (dh * dh + dv * dv + dd * dd).sqrt()
    }
}

// ── Physics-based indicators ─────────────────────────────────────────────────

impl AdaptiveMeshRefinement {
    /// Compute physics-based refinement indicator
    pub(super) fn compute_physics_based_indicator(
        solution: &DMatrix<f64>,
        i: usize,
        j: usize,
        phenomena: &[PhysicsPhenomena],
    ) -> f64 {
        let mut indicator: f64 = 0.0;

        for phenomenon in phenomena {
            let phenomenon_indicator = match phenomenon {
                PhysicsPhenomena::ShockWaves => Self::compute_shock_wave_indicator(solution, i, j),
                PhysicsPhenomena::BoundaryLayers { wall_distance } => {
                    Self::compute_boundary_layer_indicator(solution, i, j, *wall_distance)
                }
                PhysicsPhenomena::Vortices {
                    vorticity_threshold,
                } => Self::compute_vortex_indicator(solution, i, j, *vorticity_threshold),
                PhysicsPhenomena::Interfaces {
                    interface_thickness,
                } => Self::compute_interface_indicator(solution, i, j, *interface_thickness),
                PhysicsPhenomena::ThermalGradients {
                    temp_gradient_threshold,
                } => Self::compute_thermal_gradient_indicator(
                    solution,
                    i,
                    j,
                    *temp_gradient_threshold,
                ),
            };

            // Combine indicators using maximum (conservative approach)
            indicator = indicator.max(phenomenon_indicator);
        }

        indicator
    }

    /// Shock wave detection using pressure jump and Mach number criteria
    fn compute_shock_wave_indicator(solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();

        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // Pressure gradient (shock indicator)
        let dp_dx = (solution[(i + 1, j)] - solution[(i - 1, j)]) / 2.0;
        let dp_dy = (solution[(i, j + 1)] - solution[(i, j - 1)]) / 2.0;
        let grad_p = (dp_dx * dp_dx + dp_dy * dp_dy).sqrt();

        // Pressure jump detection
        let p_center = solution[(i, j)];
        let p_left = solution[(i - 1, j)];
        let p_right = solution[(i + 1, j)];
        let pressure_jump = ((p_right - p_left).abs() / (2.0 * p_center.abs() + 1e-10)).max(0.0);

        // Density jump (assuming solution represents density)
        let density_jump = Self::compute_discontinuity_indicator(solution, i, j);

        // Combined shock indicator
        grad_p * (1.0 + pressure_jump + density_jump)
    }

    /// Boundary layer detection based on wall distance and velocity gradients
    fn compute_boundary_layer_indicator(
        solution: &DMatrix<f64>,
        i: usize,
        j: usize,
        wall_distance: f64,
    ) -> f64 {
        let (nx, ny) = solution.shape();

        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // Distance from nearest wall (simplified - using grid coordinates)
        let wall_dist = (i as f64)
            .min(j as f64)
            .min((nx - i - 1) as f64)
            .min((ny - j - 1) as f64);

        // Only apply boundary layer detection near walls
        if wall_dist > wall_distance * 3.0 {
            return 0.0;
        }

        // Velocity gradient magnitude (boundary layer thickness indicator)
        let dv_dx = (solution[(i + 1, j)] - solution[(i - 1, j)]) / 2.0;
        let dv_dy = (solution[(i, j + 1)] - solution[(i, j - 1)]) / 2.0;
        let grad_v = (dv_dx * dv_dx + dv_dy * dv_dy).sqrt();

        // Wall damping function (stronger near wall)
        let wall_factor = (-wall_dist / wall_distance).exp();

        // Boundary layer indicator
        grad_v * wall_factor
    }

    /// Vortex detection using vorticity magnitude and Q-criterion
    fn compute_vortex_indicator(
        solution: &DMatrix<f64>,
        i: usize,
        j: usize,
        vorticity_threshold: f64,
    ) -> f64 {
        let (nx, ny) = solution.shape();

        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // Compute vorticity
        let vorticity = Self::compute_vorticity_indicator(solution, i, j);

        if vorticity < vorticity_threshold {
            return 0.0;
        }

        // Q-criterion (second invariant of velocity gradient tensor)
        // Simplified version for 2D: Q = 0.5 * (||Ω||² - ||S||²)
        let du_dx = (solution[(i + 1, j)] - solution[(i - 1, j)]) / 2.0;
        let du_dy = (solution[(i, j + 1)] - solution[(i, j - 1)]) / 2.0;
        let dv_dx = (solution[(i + 1, j)] - solution[(i - 1, j)]) / 2.0;
        let dv_dy = (solution[(i, j + 1)] - solution[(i, j - 1)]) / 2.0;

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
    fn compute_interface_indicator(
        solution: &DMatrix<f64>,
        i: usize,
        j: usize,
        interface_thickness: f64,
    ) -> f64 {
        let (nx, ny) = solution.shape();

        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // Volume fraction gradient (interface indicator)
        let dphi_dx = (solution[(i + 1, j)] - solution[(i - 1, j)]) / 2.0;
        let dphi_dy = (solution[(i, j + 1)] - solution[(i, j - 1)]) / 2.0;
        let grad_phi = (dphi_dx * dphi_dx + dphi_dy * dphi_dy).sqrt();

        // Interface detection based on gradient magnitude
        let interface_indicator = if grad_phi > 1.0 / interface_thickness {
            grad_phi * interface_thickness
        } else {
            0.0
        };

        // Curvature of interface (second derivative)
        let d2phi_dx2 = solution[(i + 1, j)] - 2.0 * solution[(i, j)] + solution[(i - 1, j)];
        let d2phi_dy2 = solution[(i, j + 1)] - 2.0 * solution[(i, j)] + solution[(i, j - 1)];
        let curvature = (d2phi_dx2 + d2phi_dy2).abs();

        // Combined interface indicator
        interface_indicator * (1.0 + curvature * interface_thickness)
    }

    /// Thermal gradient detection for heat transfer problems
    fn compute_thermal_gradient_indicator(
        solution: &DMatrix<f64>,
        i: usize,
        j: usize,
        temp_gradient_threshold: f64,
    ) -> f64 {
        let (nx, ny) = solution.shape();

        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // Temperature gradient magnitude
        let dt_dx = (solution[(i + 1, j)] - solution[(i - 1, j)]) / 2.0;
        let dt_dy = (solution[(i, j + 1)] - solution[(i, j - 1)]) / 2.0;
        let grad_t = (dt_dx * dt_dx + dt_dy * dt_dy).sqrt();

        if grad_t < temp_gradient_threshold {
            return 0.0;
        }

        // Heat flux divergence (indicator of thermal sources/sinks)
        let d2t_dx2 = solution[(i + 1, j)] - 2.0 * solution[(i, j)] + solution[(i - 1, j)];
        let d2t_dy2 = solution[(i, j + 1)] - 2.0 * solution[(i, j)] + solution[(i, j - 1)];
        let heat_source = (d2t_dx2 + d2t_dy2).abs();

        // Thermal boundary layer detection
        grad_t * (1.0 + heat_source / temp_gradient_threshold)
    }
}
