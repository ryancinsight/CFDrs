use super::config::DESConfig;
use super::super::boundary_conditions::TurbulenceBoundaryCondition;
use super::super::spalart_allmaras::SpalartAllmaras;
use cfd_core::error::Error;
use leto::geometry::Vector2;
use leto::Array2;
use std::f64;

/// Detached Eddy Simulation model
#[derive(Debug)]
pub struct DetachedEddySimulation {
    /// DES configuration
    pub(super) config: DESConfig,
    /// SGS viscosity field (LES mode)
    pub(super) sgs_viscosity: Array2<f64>,
    /// Modified turbulent viscosity field (RANS state variable)
    pub(super) nu_tilde: Array2<f64>,
    /// Underlying Spalart-Allmaras RANS model
    pub(super) spalart_allmaras: SpalartAllmaras<f64>,
    /// DES length scale field
    pub(super) des_length_scale: Array2<f64>,
    /// Wall distance field (for shielding functions)
    pub(super) wall_distance: Array2<f64>,
    /// Reusable buffer for velocity field adaptation (SoA to AoS)
    pub(super) velocity_buffer: Vec<Vector2<f64>>,
    /// Grid spacing in x-direction
    pub(super) dx: f64,
    /// Grid spacing in y-direction
    pub(super) dy: f64,
}

impl DetachedEddySimulation {
    /// Create a new DES model.
    ///
    /// # Panics
    /// Panics if any invariant is violated (see [`Self::try_new`]).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nx: usize,
        ny: usize,
        dx: f64,
        dy: f64,
        config: DESConfig,
        boundaries: &[(&str, TurbulenceBoundaryCondition<f64>)],
    ) -> Self {
        Self::try_new(nx, ny, dx, dy, config, boundaries).unwrap_or_else(|error| {
            panic!("DetachedEddySimulation::new called with invalid arguments: {error}");
        })
    }

    /// Create a new DES model with invariant validation.
    ///
    /// Validates:
    /// - `nx >= 1`, `ny >= 1`.
    /// - `dx > 0`, `dy > 0` (DES length scale is `Δ = max(dx, dy, dy,
    ///   d_w)`; zero or non-finite spacing collapses it).
    /// - `des_constant > 0` and finite (the RANS-LES shielding threshold
    ///   `r̃_d = C_DES · Δ_w / Δ_max` requires `C_DES > 0`).
    /// - `rans_viscosity > 0` and finite.
    ///
    /// # Errors
    /// Returns `Error::InvalidConfiguration` if any invariant is violated.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        nx: usize,
        ny: usize,
        dx: f64,
        dy: f64,
        config: DESConfig,
        boundaries: &[(&str, TurbulenceBoundaryCondition<f64>)],
    ) -> cfd_core::error::Result<Self> {
        if nx == 0 {
            return Err(Error::InvalidConfiguration(
                "DetachedEddySimulation::try_new: nx must be at least 1".to_string(),
            ));
        }
        if ny == 0 {
            return Err(Error::InvalidConfiguration(
                "DetachedEddySimulation::try_new: ny must be at least 1".to_string(),
            ));
        }
        if !dx.is_finite() || dx <= 0.0 {
            return Err(Error::InvalidConfiguration(format!(
                "DetachedEddySimulation::try_new: dx must be finite and positive, got {dx:?}"
            )));
        }
        if !dy.is_finite() || dy <= 0.0 {
            return Err(Error::InvalidConfiguration(format!(
                "DetachedEddySimulation::try_new: dy must be finite and positive, got {dy:?}"
            )));
        }
        if !config.des_constant.is_finite() || config.des_constant <= 0.0 {
            return Err(Error::InvalidConfiguration(format!(
                "DetachedEddySimulation::try_new: des_constant must be finite and positive, got {:?}",
                config.des_constant
            )));
        }
        if !config.rans_viscosity.is_finite() || config.rans_viscosity <= 0.0 {
            return Err(Error::InvalidConfiguration(format!(
                "DetachedEddySimulation::try_new: rans_viscosity must be finite and positive, got {:?}",
                config.rans_viscosity
            )));
        }
        let mut wall_distance = Array2::from_elem([nx, ny], f64::MAX);

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

                wall_distance[[i, j]] = d;
            }
        }

        let spalart_allmaras = SpalartAllmaras::new(nx, ny);

        Ok(Self {
            config,
            sgs_viscosity: Array2::zeros([nx, ny]),
            nu_tilde: Array2::zeros([nx, ny]),
            spalart_allmaras,
            des_length_scale: Array2::zeros([nx, ny]),
            wall_distance,
            velocity_buffer: vec![Vector2::zeros(); nx * ny],
            dx,
            dy,
        })
    }
}
// Additional methods for DES-specific functionality
impl DetachedEddySimulation {
    /// Get the DES length scale field
    pub fn get_des_length_scale_field(&self) -> &Array2<f64> {
        &self.des_length_scale
    }

    /// Get the SGS viscosity field
    pub fn get_sgs_viscosity_field(&self) -> &Array2<f64> {
        &self.sgs_viscosity
    }

    /// Check if a point is in LES mode (DES length scale active)
    pub fn is_les_mode(&self, i: usize, j: usize) -> bool {
        // LES mode when DES length scale is smaller than grid scale
        let des_length = self.des_length_scale[[i, j]];

        // Use local grid scale (dx, dy) or Δ definition consistent with DES formulation.
        // For DES97/DDES, Δ = max(dx, dy)
        let grid_scale = self.dx.max(self.dy);

        des_length < grid_scale
    }
}
