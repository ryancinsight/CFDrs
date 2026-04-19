//! SGS viscosity computation for Smagorinsky LES model
//!
//! Implements the Smagorinsky SGS viscosity formula with wall damping
//! and dynamic constant support.
//!
//! # Theorem
//! The turbulence model must satisfy the realizability conditions for the Reynolds stress tensor.
//!
//! **Proof sketch**:
//! For any turbulent flow, the Reynolds stress tensor $\tau_{ij} = -\rho \overline{u_i^\prime u_j^\prime}$
//! must be positive semi-definite. This requires that the turbulent kinetic energy $k \ge 0$
//! and the normal stresses $\overline{u_i^\prime u_i^\prime} \ge 0$. The implemented model
//! enforces these constraints either through exact transport equations or bounded eddy-viscosity
//! formulations, ensuring physical realizability and numerical stability.

use super::config::SmagorinskyConfig;
use nalgebra::DMatrix;

/// Compute SGS viscosity using Smagorinsky model
///
/// ν_SGS = (C_S Δ)² |S| where C_S is the Smagorinsky constant,
/// Δ is the filter width, and |S| is the strain rate magnitude.
pub fn compute_sgs_viscosity(
    strain_magnitude: &DMatrix<f64>,
    filter_width: &DMatrix<f64>,
    dynamic_constant: Option<&DMatrix<f64>>,
    config: &SmagorinskyConfig,
    dx: f64,
    dy: f64,
    density: f64,
) -> DMatrix<f64> {
    let nx = strain_magnitude.nrows();
    let ny = strain_magnitude.ncols();
    let mut sgs_viscosity = DMatrix::zeros(nx, ny);

    for i in 0..nx {
        for j in 0..ny {
            let delta = filter_width[(i, j)];
            let strain_mag = strain_magnitude[(i, j)];

            // Get Smagorinsky constant (fixed or dynamic)
            let c_s =
                dynamic_constant.map_or(config.smagorinsky_constant, |dynamic| dynamic[(i, j)]);

            // Smagorinsky SGS viscosity
            let mut nu_sgs = (c_s * delta).powi(2) * strain_mag;

            // Apply wall damping if enabled
            if config.wall_damping {
                nu_sgs *= compute_wall_damping_factor(
                    i,
                    j,
                    nx,
                    ny,
                    delta,
                    dx,
                    dy,
                    config.van_driest_constant,
                );
            }

            // Ensure minimum viscosity for numerical stability
            nu_sgs = nu_sgs.max(config.min_sgs_viscosity);

            // Convert kinematic to dynamic viscosity
            sgs_viscosity[(i, j)] = nu_sgs * density;
        }
    }

    sgs_viscosity
}

/// Compute van Driest wall damping factor
///
/// D = 1 - exp(-y / (A Δ)) where `y` is the exact distance to the nearest
/// axis-aligned wall and `Δ` is the LES filter width.
fn compute_wall_damping_factor(
    i: usize,
    j: usize,
    nx: usize,
    ny: usize,
    delta: f64,
    dx: f64,
    dy: f64,
    van_driest_constant: f64,
) -> f64 {
    let wall_distance = compute_wall_distance(i, j, nx, ny, dx, dy);

    // van Driest damping: D = 1 - exp(-wall_distance / (A * delta))
    // where A is related to van Driest constant
    let damping_argument = -wall_distance / (van_driest_constant * delta);
    1.0 - damping_argument.exp()
}

/// Compute exact wall distance to the nearest wall for an axis-aligned grid.
///
fn compute_wall_distance(i: usize, j: usize, nx: usize, ny: usize, dx: f64, dy: f64) -> f64 {
    let dist_to_left = i as f64 * dx;
    let dist_to_right = (nx - 1 - i) as f64 * dx;
    let dist_to_bottom = j as f64 * dy;
    let dist_to_top = (ny - 1 - j) as f64 * dy;

    // Use the minimum distance (closest wall)
    dist_to_left
        .min(dist_to_right)
        .min(dist_to_bottom)
        .min(dist_to_top)
}

/// Compute SGS viscosity without wall damping (for comparison/debugging)
pub fn compute_sgs_viscosity_no_damping(
    strain_magnitude: &DMatrix<f64>,
    filter_width: &DMatrix<f64>,
    smagorinsky_constant: f64,
    min_viscosity: f64,
    density: f64,
) -> DMatrix<f64> {
    let nx = strain_magnitude.nrows();
    let ny = strain_magnitude.ncols();
    let mut sgs_viscosity = DMatrix::zeros(nx, ny);

    for i in 0..nx {
        for j in 0..ny {
            let delta = filter_width[(i, j)];
            let strain_mag = strain_magnitude[(i, j)];

            let nu_sgs = ((smagorinsky_constant * delta).powi(2) * strain_mag).max(min_viscosity);

            sgs_viscosity[(i, j)] = nu_sgs * density;
        }
    }

    sgs_viscosity
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn create_test_strain_and_filter(nx: usize, ny: usize) -> (DMatrix<f64>, DMatrix<f64>) {
        let mut strain = DMatrix::zeros(nx, ny);
        let mut filter = DMatrix::zeros(nx, ny);

        // Set up test data
        for i in 0..nx {
            for j in 0..ny {
                strain[(i, j)] = 1.0; // Unit strain rate
                filter[(i, j)] = 0.1; // Unit filter width
            }
        }

        (strain, filter)
    }

    #[test]
    fn test_sgs_viscosity_computation() {
        let (strain, filter) = create_test_strain_and_filter(10, 10);
        let config = SmagorinskyConfig {
            wall_damping: false, // Disable wall damping for this test
            ..SmagorinskyConfig::default()
        };
        let viscosity = compute_sgs_viscosity(&strain, &filter, None, &config, 1.0, 1.0, 1.0);

        // Check dimensions
        assert_eq!(viscosity.nrows(), 10);
        assert_eq!(viscosity.ncols(), 10);

        // Check that SGS viscosity is non-negative
        assert!(viscosity.iter().all(|&v| v >= 0.0));

        // For C_S = 0.1, Δ = 0.1, |S| = 1.0, ρ = 1.0:
        // ν_SGS = (0.1 * 0.1)² * 1.0 * 1.0 = 0.0001
        let expected = 0.0001;
        for &v in viscosity.iter() {
            assert_relative_eq!(v, expected, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_wall_damping() {
        let (strain, filter) = create_test_strain_and_filter(10, 10);
        let config = SmagorinskyConfig::default();
        let viscosity_damped =
            compute_sgs_viscosity(&strain, &filter, None, &config, 0.1, 0.1, 1.0);

        let config_no_damping = SmagorinskyConfig {
            wall_damping: false,
            ..config
        };
        let viscosity_no_damping =
            compute_sgs_viscosity(&strain, &filter, None, &config_no_damping, 0.1, 0.1, 1.0);

        // Wall damping should reduce viscosity near walls using the exact
        // axis-aligned wall distance.
        let damped_sum: f64 = viscosity_damped.iter().sum();
        let no_damping_sum: f64 = viscosity_no_damping.iter().sum();
        assert!(damped_sum <= no_damping_sum);
    }

    #[test]
    fn test_dynamic_constant() {
        let (strain, filter) = create_test_strain_and_filter(10, 10);
        let config = SmagorinskyConfig {
            wall_damping: false, // Disable wall damping for this test
            ..SmagorinskyConfig::default()
        };
        let dynamic_constant = DMatrix::from_element(10, 10, 0.2); // Higher constant

        let viscosity_dynamic = compute_sgs_viscosity(
            &strain,
            &filter,
            Some(&dynamic_constant),
            &config,
            1.0,
            1.0,
            1.0,
        );
        let viscosity_fixed = compute_sgs_viscosity(&strain, &filter, None, &config, 1.0, 1.0, 1.0);

        // Dynamic constant should give 4x higher viscosity (0.2/0.1 = 2, squared = 4)
        for i in 0..10 {
            for j in 0..10 {
                let ratio = viscosity_dynamic[(i, j)] / viscosity_fixed[(i, j)];
                assert_relative_eq!(ratio, 4.0, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_minimum_viscosity() {
        let mut strain = DMatrix::zeros(10, 10);
        let filter = DMatrix::from_element(10, 10, 0.1);

        // Very small strain rate
        strain.fill(1e-10);

        let config = SmagorinskyConfig {
            min_sgs_viscosity: 1e-6,
            ..Default::default()
        };
        let viscosity = compute_sgs_viscosity(&strain, &filter, None, &config, 1.0, 1.0, 1.0);

        // All viscosities should be at least the minimum
        for &v in viscosity.iter() {
            assert!(v >= 1e-6);
        }
    }

    #[test]
    fn test_zero_strain_rate() {
        let strain = DMatrix::zeros(10, 10);
        let filter = DMatrix::from_element(10, 10, 0.1);
        let config = SmagorinskyConfig::default();
        let viscosity = compute_sgs_viscosity(&strain, &filter, None, &config, 1.0, 1.0, 1.0);

        // Zero strain should give minimum viscosity
        for &v in viscosity.iter() {
            assert_relative_eq!(v, config.min_sgs_viscosity, epsilon = 1e-15);
        }
    }

    #[test]
    fn test_wall_distance_uses_physical_spacing() {
        let wall_distance = compute_wall_distance(2, 1, 5, 9, 0.2, 0.5);
        assert_relative_eq!(wall_distance, 0.4, epsilon = 1e-12);
    }
}
