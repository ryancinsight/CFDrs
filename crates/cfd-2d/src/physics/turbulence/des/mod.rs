//! Detached Eddy Simulation (DES) model
//!
//! DES is a hybrid RANS-LES approach where RANS models are used near walls
//! and LES is used in the detached regions away from walls. This provides
//! accurate boundary layer prediction with reduced computational cost.
//!
//! ## DES97 Formulation
//!
//! The DES length scale is defined as:
//!
//! L_DES = min(L_RANS, C_DES * Δ)
//!
//! where L_RANS is the RANS length scale, C_DES ≈ 0.65, and Δ is the grid spacing.
//!
//! ## DDES (Delayed DES)
//!
//! Delayed DES prevents premature switching to LES in boundary layers by
//! modifying the length scale computation to account for the wall distance.
//!
//! ## IDDES (Improved DDES)
//!
//! IDDES further improves the shielding function and adds a wall-modeled LES
//! capability for very fine grids near walls.
//!
//! ## References
//!
//! - Spalart, P. R., et al. (1997). Comments on the feasibility of LES for wings.
//! - Spalart, P. R., et al. (2006). A new version of detached-eddy simulation.
//! - Shur, M. L., et al. (2008). A hybrid RANS-LES approach with delayed-DES.
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

// Dependencies: RANS model integration, wall distance computation, length scale consistency
// Mathematical Foundation: Spalart et al. (1997) DES97, Shur et al. (2008) IDDES

mod length_scale;

use super::boundary_conditions::TurbulenceBoundaryCondition;
use super::spalart_allmaras::SpalartAllmaras;
use super::traits::LESTurbulenceModel;
use nalgebra::{DMatrix, Vector2};
use std::f64;

#[cfg(feature = "gpu")]
use cfd_core::compute::gpu::turbulence_compute::GpuTurbulenceCompute;

/// DES model variants
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DESVariant {
    /// Original DES97
    DES97,
    /// Delayed DES (DDES)
    DDES,
    /// Improved DDES (IDDES)
    IDDES,
}

/// DES configuration
#[derive(Debug, Clone)]
pub struct DESConfig {
    /// DES variant to use
    pub variant: DESVariant,
    /// DES constant (C_DES)
    pub des_constant: f64,
    /// Maximum SGS viscosity ratio
    pub max_sgs_ratio: f64,
    /// Molecular viscosity (constant)
    pub rans_viscosity: f64,
    /// Enable GPU acceleration
    pub use_gpu: bool,
}

impl Default for DESConfig {
    fn default() -> Self {
        Self {
            variant: DESVariant::DDES,
            des_constant: 0.65,
            max_sgs_ratio: 0.5,   // Prevent excessive SGS viscosity
            rans_viscosity: 1e-5, // Default molecular viscosity
            use_gpu: false,
        }
    }
}

/// Detached Eddy Simulation model
#[derive(Debug)]
pub struct DetachedEddySimulation {
    /// DES configuration
    pub(super) config: DESConfig,
    /// SGS viscosity field (LES mode)
    pub(super) sgs_viscosity: DMatrix<f64>,
    /// Modified turbulent viscosity field (RANS state variable)
    pub(super) nu_tilde: DMatrix<f64>,
    /// Underlying Spalart-Allmaras RANS model
    pub(super) spalart_allmaras: SpalartAllmaras<f64>,
    /// DES length scale field
    pub(super) des_length_scale: DMatrix<f64>,
    /// Wall distance field (for shielding functions)
    pub(super) wall_distance: DMatrix<f64>,
    /// Reusable buffer for velocity field adaptation (SoA to AoS)
    velocity_buffer: Vec<Vector2<f64>>,
    /// GPU compute manager (if GPU acceleration is enabled)
    #[cfg(feature = "gpu")]
    gpu_compute: Option<GpuTurbulenceCompute>,
    /// Grid spacing in x-direction
    pub(super) dx: f64,
    /// Grid spacing in y-direction
    pub(super) dy: f64,
}

impl DetachedEddySimulation {
    /// Create a new DES model
    pub fn new(
        nx: usize,
        ny: usize,
        dx: f64,
        dy: f64,
        config: DESConfig,
        boundaries: &[(&str, TurbulenceBoundaryCondition<f64>)],
    ) -> Self {
        let mut wall_distance = DMatrix::from_element(nx, ny, f64::MAX);

        let has_boundaries = !boundaries.is_empty();

        // If no boundaries provided, assume all walls (fallback)
        let process_west = !has_boundaries
            || boundaries.iter().any(|(n, bc)| {
                *n == "west" && matches!(bc, TurbulenceBoundaryCondition::Wall { .. })
            });
        let process_east = !has_boundaries
            || boundaries.iter().any(|(n, bc)| {
                *n == "east" && matches!(bc, TurbulenceBoundaryCondition::Wall { .. })
            });
        let process_south = !has_boundaries
            || boundaries.iter().any(|(n, bc)| {
                *n == "south" && matches!(bc, TurbulenceBoundaryCondition::Wall { .. })
            });
        let process_north = !has_boundaries
            || boundaries.iter().any(|(n, bc)| {
                *n == "north" && matches!(bc, TurbulenceBoundaryCondition::Wall { .. })
            });

        for i in 0..nx {
            for j in 0..ny {
                let mut d = f64::MAX;

                if process_west {
                    d = d.min((i as f64 + 0.5) * dx);
                }
                if process_east {
                    d = d.min((nx as f64 - 1.0 - i as f64 + 0.5) * dx);
                }
                if process_south {
                    d = d.min((j as f64 + 0.5) * dy);
                }
                if process_north {
                    d = d.min((ny as f64 - 1.0 - j as f64 + 0.5) * dy);
                }

                wall_distance[(i, j)] = d;
            }
        }

        #[cfg(feature = "gpu")]
        let gpu_compute = if config.use_gpu {
            match GpuTurbulenceCompute::new() {
                Ok(compute) => Some(compute),
                Err(e) => {
                    tracing::warn!("Failed to initialize GPU compute for DES: {}", e);
                    None
                }
            }
        } else {
            None
        };

        let spalart_allmaras = SpalartAllmaras::new(nx, ny);

        Self {
            config,
            sgs_viscosity: DMatrix::zeros(nx, ny),
            nu_tilde: DMatrix::zeros(nx, ny),
            spalart_allmaras,
            des_length_scale: DMatrix::zeros(nx, ny),
            wall_distance,
            velocity_buffer: vec![Vector2::zeros(); nx * ny],
            #[cfg(feature = "gpu")]
            gpu_compute,
            dx,
            dy,
        }
    }
}

impl LESTurbulenceModel for DetachedEddySimulation {
    fn update(
        &mut self,
        velocity_u: &DMatrix<f64>,
        velocity_v: &DMatrix<f64>,
        _pressure: &DMatrix<f64>,
        _density: f64,
        _viscosity: f64,
        dt: f64,
        dx: f64,
        dy: f64,
    ) -> cfd_core::error::Result<()> {
        self.dx = dx;
        self.dy = dy;

        // Populate velocity buffer (SoA -> AoS adaptation)
        let nx = velocity_u.nrows();
        let ny = velocity_v.ncols();
        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;
                self.velocity_buffer[idx] = Vector2::new(velocity_u[(i, j)], velocity_v[(i, j)]);
            }
        }

        // Compute DES length scale
        self.des_length_scale = self.compute_des_length_scale(velocity_u, velocity_v, dx, dy);

        // Update nu_tilde using Spalart-Allmaras solver with DES length scale override
        // This effectively implements the DES formulation where the destruction term
        // uses min(d, C_DES * Delta) instead of d.
        self.spalart_allmaras.update_with_distance(
            self.nu_tilde.as_mut_slice(),
            &self.velocity_buffer,
            self.config.rans_viscosity, // Use config value as molecular viscosity
            self.des_length_scale.as_slice(),
            dx,
            dy,
            dt,
        )?;

        // Update SGS/Turbulent viscosity field
        for j in 0..ny {
            for i in 0..nx {
                let nu_t = self
                    .spalart_allmaras
                    .eddy_viscosity(self.nu_tilde[(i, j)], self.config.rans_viscosity);

                // Limit SGS viscosity if needed (though SA usually stable)
                let max_visc = self.config.max_sgs_ratio * self.config.rans_viscosity;
                self.sgs_viscosity[(i, j)] = nu_t.min(max_visc);
            }
        }

        Ok(())
    }

    fn get_viscosity(&self, i: usize, j: usize) -> f64 {
        // Return total viscosity (molecular + SGS)
        self.config.rans_viscosity + self.sgs_viscosity[(i, j)]
    }

    fn get_turbulent_viscosity_field(&self) -> &DMatrix<f64> {
        // Return SGS viscosity field (turbulent viscosity only)
        // The total viscosity = molecular + turbulent is computed in get_viscosity()
        &self.sgs_viscosity
    }

    fn get_turbulent_kinetic_energy(&self, i: usize, j: usize) -> f64 {
        // For DES, TKE is typically obtained from the underlying RANS model
        // Since this implementation doesn't include RANS coupling, we estimate
        // TKE from the SGS viscosity using the relation: k ≈ (ν_sgs / C_k)^{2/3}
        // where C_k is a constant related to the SGS model

        let nu_sgs = self.sgs_viscosity[(i, j)];
        if nu_sgs > 0.0 {
            // Estimate k from SGS viscosity: k ≈ (ν_sgs / C_k)^{2/3}
            // Using C_k ≈ 0.1 (typical value for Smagorinsky-based models)
            let c_k = 0.1;
            (nu_sgs / c_k).powf(2.0 / 3.0)
        } else {
            0.0
        }
    }

    fn get_dissipation_rate(&self, i: usize, j: usize) -> f64 {
        // For DES, dissipation rate ε is related to TKE and length scale: ε = k^{3/2} / l
        // This provides a physically consistent estimate

        let k = self.get_turbulent_kinetic_energy(i, j);
        let l_des = self.des_length_scale[(i, j)];

        if k > 0.0 && l_des > 0.0 {
            k.powf(1.5) / l_des
        } else {
            0.0
        }
    }

    fn boundary_condition_update(
        &mut self,
        _boundary_manager: &super::boundary_conditions::TurbulenceBoundaryManager<f64>,
    ) -> cfd_core::error::Result<()> {
        // DES typically uses homogeneous boundary conditions
        // or recycling/rescaling methods for inflow
        Ok(())
    }

    fn get_model_name(&self) -> &str {
        match self.config.variant {
            DESVariant::DES97 => "DES97",
            DESVariant::DDES => "Delayed DES",
            DESVariant::IDDES => "Improved DDES",
        }
    }

    fn get_model_constants(&self) -> Vec<(&str, f64)> {
        let constants = vec![
            ("DES Constant", self.config.des_constant),
            ("Max SGS Ratio", self.config.max_sgs_ratio),
        ];

        // DES-specific constants only (RANS constants would come from separate model)

        constants
    }
}

// Additional methods for DES-specific functionality
impl DetachedEddySimulation {
    /// Get the DES length scale field
    pub fn get_des_length_scale_field(&self) -> &DMatrix<f64> {
        &self.des_length_scale
    }

    /// Get the SGS viscosity field
    pub fn get_sgs_viscosity_field(&self) -> &DMatrix<f64> {
        &self.sgs_viscosity
    }

    /// Check if a point is in LES mode (DES length scale active)
    pub fn is_les_mode(&self, i: usize, j: usize) -> bool {
        // LES mode when DES length scale is smaller than grid scale
        let des_length = self.des_length_scale[(i, j)];

        // Use local grid scale (dx, dy) or Δ definition consistent with DES formulation.
        // For DES97/DDES, Δ = max(dx, dy)
        let grid_scale = self.dx.max(self.dy);

        des_length < grid_scale
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_fields(nx: usize, ny: usize) -> (DMatrix<f64>, DMatrix<f64>, DMatrix<f64>) {
        let mut velocity_u = DMatrix::zeros(nx, ny);
        let velocity_v = DMatrix::zeros(nx, ny);
        let pressure = DMatrix::zeros(nx, ny);

        // Simple shear flow
        for i in 0..nx {
            for j in 0..ny {
                velocity_u[(i, j)] = (j as f64) * 0.1;
            }
        }

        (velocity_u, velocity_v, pressure)
    }

    #[test]
    fn test_des_creation() {
        let config = DESConfig::default();
        let des = DetachedEddySimulation::new(10, 10, 0.1, 0.1, config, &[]);

        assert_eq!(des.sgs_viscosity.nrows(), 10);
        assert_eq!(des.sgs_viscosity.ncols(), 10);
        assert_eq!(des.des_length_scale.nrows(), 10);
        assert_eq!(des.des_length_scale.ncols(), 10);
    }

    #[test]
    fn test_des_variants() {
        let variants = vec![DESVariant::DES97, DESVariant::DDES, DESVariant::IDDES];

        for variant in variants {
            let config = DESConfig {
                variant,
                ..Default::default()
            };
            let des = DetachedEddySimulation::new(10, 10, 0.1, 0.1, config, &[]);

            let name = des.get_model_name();
            assert!(name.contains("DES"));
        }
    }

    #[test]
    fn test_des_length_scale_computation() {
        let config = DESConfig {
            variant: DESVariant::DES97, // Use simpler DES97 for testing
            ..DESConfig::default()
        };
        let des = DetachedEddySimulation::new(10, 10, 0.1, 0.1, config, &[]);

        // Create dummy velocity fields
        let velocity_u = DMatrix::from_element(10, 10, 1.0);
        let velocity_v = DMatrix::from_element(10, 10, 0.5);

        let length_scale = des.compute_des_length_scale(&velocity_u, &velocity_v, 0.1, 0.1);

        // Check dimensions
        assert_eq!(length_scale.nrows(), 10);
        assert_eq!(length_scale.ncols(), 10);

        // Length scale should be positive and reasonable
        assert!(length_scale.iter().all(|&l| l > 0.0 && l.is_finite()));
    }

    #[test]
    fn test_des_model_update() {
        let config = DESConfig::default();
        let mut des = DetachedEddySimulation::new(10, 10, 0.1, 0.1, config, &[]);
        let (velocity_u, velocity_v, pressure) = create_test_fields(10, 10);

        let result = des.update(
            &velocity_u,
            &velocity_v,
            &pressure,
            1.0,
            0.01,
            0.001,
            0.1,
            0.1,
        );
        assert!(result.is_ok());
    }

    #[test]
    fn test_des_mode_detection() {
        let config = DESConfig::default();
        let des = DetachedEddySimulation::new(10, 10, 0.1, 0.1, config, &[]);

        // Test a point - should work without panicking
        let is_les = des.is_les_mode(5, 5);
        // Result depends on computed fields, just check it doesn't crash
        let _ = is_les;
    }

    #[test]
    fn test_des_model_constants() {
        let config = DESConfig {
            des_constant: 0.7,
            max_sgs_ratio: 0.6,
            ..Default::default()
        };
        let des = DetachedEddySimulation::new(10, 10, 0.1, 0.1, config, &[]);

        let constants = des.get_model_constants();
        assert!(!constants.is_empty());

        // Should include DES-specific constants
        assert!(constants.iter().any(|(name, _)| *name == "DES Constant"));
        assert!(constants.iter().any(|(name, _)| *name == "Max SGS Ratio"));
    }

    #[test]
    fn test_iddes_length_scale_computation() {
        let config = DESConfig {
            variant: DESVariant::IDDES,
            ..Default::default()
        };
        let des = DetachedEddySimulation::new(10, 10, 0.1, 0.1, config, &[]);

        let velocity_u = DMatrix::from_element(10, 10, 1.0);
        let velocity_v = DMatrix::from_element(10, 10, 0.5);

        let length_scale = des.compute_des_length_scale(&velocity_u, &velocity_v, 0.1, 0.1);

        assert_eq!(length_scale.nrows(), 10);
        assert_eq!(length_scale.ncols(), 10);
        assert!(length_scale.iter().all(|&l| l > 0.0 && l.is_finite()));
    }

    #[test]
    fn test_ddes_shielding_behavior() {
        // Setup DDES
        let config = DESConfig {
            variant: DESVariant::DDES,
            des_constant: 0.65,
            ..Default::default()
        };

        let mut des = DetachedEddySimulation::new(10, 10, 1.0, 1.0, config, &[]);

        // Manually overwrite wall distance for a specific point (5, 5)
        // Set d_w = 1.0.
        des.wall_distance[(5, 5)] = 1.0;

        // l_RANS = 1.0.
        // l_LES = 0.65 * 1.0 = 0.65.
        // l_RANS > l_LES, so shielding logic enters.

        // Case 1: High Eddy Viscosity (Shielding Active -> RANS mode)
        // nu_tilde needs to be high enough.
        des.nu_tilde[(5, 5)] = 10.0;

        // Create velocity field that gives S=1.0 at (5,5).
        // Simple shear: u = y. du/dy = 1.
        let mut velocity_u = DMatrix::zeros(10, 10);
        for j in 0..10 {
            for i in 0..10 {
                velocity_u[(i, j)] = j as f64;
            }
        }
        let velocity_v = DMatrix::zeros(10, 10);

        let length_scale = des.compute_des_length_scale(&velocity_u, &velocity_v, 1.0, 1.0);

        let l = length_scale[(5, 5)];
        // Expect RANS length scale (1.0) because of shielding
        assert!(
            (l - 1.0).abs() < 1e-3,
            "With high viscosity, DDES should shield and return RANS length (1.0), got {}",
            l
        );

        // Case 2: Zero Eddy Viscosity (Shielding Inactive -> LES mode)
        des.nu_tilde[(5, 5)] = 0.0;
        let length_scale_les = des.compute_des_length_scale(&velocity_u, &velocity_v, 1.0, 1.0);
        let l_les = length_scale_les[(5, 5)];

        // Expect LES length scale (0.65)
        assert!(
            (l_les - 0.65).abs() < 1e-3,
            "With zero viscosity, DDES should return LES length (0.65), got {}",
            l_les
        );
    }
}
