//! SIMPLEC and PIMPLE solver — struct definition and construction
//!
//! # Mathematical Foundations
//!
//! ## SIMPLEC Algorithm (Van Doormaal & Raithby, 1984)
//!
//! Addresses the pressure-velocity decoupling issue in incompressible Navier-Stokes:
//!
//! `∂u/∂t − ν∇²u + (u·∇)u = −∇p/ρ + f` (momentum)
//! `∇·u = 0` (continuity)
//!
//! ## PIMPLE Algorithm (Issa, 1986)
//!
//! Merges PISO outer correctors with SIMPLE inner corrections for transient flows.
//!
//! # Theorem (Rhie-Chow Consistency — Rhie & Chow 1983)
//!
//! On colocated grids with positive momentum diagonals, face velocities must be
//! interpolated consistently with the momentum equation stencil to suppress
//! checkerboard pressure modes.
//! See [`interpolation`] for the complete proof.
//!
//! # Theorem (SIMPLEC Convergence — Van Doormaal & Raithby 1984)
//!
//! Under diagonal dominance of the momentum matrix and a solved pressure
//! correction, SIMPLEC converges linearly at a rate controlled by `α_u` and
//! `α_p`. See [`algorithms`] for the full derivation.
//!
//! ## References
//! - Patankar & Spalding (1972). *Int. J. Heat Mass Transfer*, 15(10), 1787–1806.
//! - Van Doormaal & Raithby (1984). *Numerical Heat Transfer*, 7(2), 147–163.
//! - Issa (1986). *J. Comput. Phys.*, 62(1), 40–65.
//! - Rhie & Chow (1983). *AIAA Journal*, 21(11), 1525–1532.
//!

use super::config::{AlgorithmType, SimplecPimpleConfig};
use crate::fields::Field2D;
use crate::grid::StructuredGrid2D;
use crate::physics::MomentumSolver;
use crate::pressure_velocity::{PressureCorrectionSolver, RhieChowInterpolation};
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};

/// SIMPLEC/PIMPLE pressure-velocity coupling solver
///
/// This solver extends the basic SIMPLE algorithm with:
/// - SIMPLEC: Consistent pressure-velocity coupling using Rhie-Chow interpolation
/// - PIMPLE: Merged PISO-SIMPLE for better convergence in transient flows
pub struct SimplecPimpleSolver<T: RealField + Copy> {
    pub(super) config: SimplecPimpleConfig<T>,
    pub(super) grid: StructuredGrid2D<T>,
    pub(super) momentum_solver: MomentumSolver<T>,
    pub(super) pressure_solver: PressureCorrectionSolver<T>,
    pub(super) rhie_chow: Option<RhieChowInterpolation<T>>,
    pub(super) iterations: usize,
    pub(super) _u_face_cache: std::cell::RefCell<Option<crate::grid::array2d::Array2D<T>>>,
    pub(super) _v_face_cache: std::cell::RefCell<Option<crate::grid::array2d::Array2D<T>>>,
    pub(super) _d_x_cache: std::cell::RefCell<Option<crate::grid::array2d::Array2D<T>>>,
    pub(super) _d_y_cache: std::cell::RefCell<Option<crate::grid::array2d::Array2D<T>>>,
    pub(super) _vel_field_cache: std::cell::RefCell<Option<crate::fields::Field2D<nalgebra::Vector2<T>>>>,
    pub(super) _cons_vel_cache: std::cell::RefCell<Option<crate::grid::array2d::Array2D<nalgebra::Vector2<T>>>>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::LowerExp>
    SimplecPimpleSolver<T>
{
    /// Create new SIMPLEC/PIMPLE solver
    pub fn new(
        grid: StructuredGrid2D<T>,
        config: SimplecPimpleConfig<T>,
    ) -> cfd_core::error::Result<Self> {
        config.validate()?;

        // Convert SimplecPimple convection scheme to Momentum convection scheme
        let momentum_convection = match config.convection_scheme {
            crate::schemes::SpatialScheme::FirstOrderUpwind => {
                crate::physics::momentum::ConvectionScheme::Upwind
            }
            crate::schemes::SpatialScheme::QuadraticUpstreamInterpolation => {
                crate::physics::momentum::ConvectionScheme::DeferredCorrectionQuick {
                    relaxation_factor: 0.7,
                }
            }
            crate::schemes::SpatialScheme::Muscl => {
                crate::physics::momentum::ConvectionScheme::TvdVanLeer {
                    relaxation_factor: 0.7,
                }
            }
            _ => crate::physics::momentum::ConvectionScheme::Upwind, // Fallback for Central/SecondOrder/WENO
        };

        let mut momentum_solver = MomentumSolver::new(&grid);
        momentum_solver.set_convection_scheme(momentum_convection);
        momentum_solver.set_velocity_relaxation(config.alpha_u);

        // Create pressure solver using configuration
        let pressure_solver =
            PressureCorrectionSolver::new(grid.clone(), config.pressure_linear_solver)?;

        let rhie_chow = if config.use_rhie_chow {
            let mut rhie_chow = RhieChowInterpolation::new(&grid);
            // Initialize Rhie-Chow coefficients with reasonable defaults
            // This prevents issues in the first iteration
            let dx_f64 = grid.dx.to_f64().unwrap_or(1.0);
            let dy_f64 = grid.dy.to_f64().unwrap_or(1.0);
            let default_ap_f64 = 1.0 / (dx_f64 * dx_f64 + dy_f64 * dy_f64);
            let default_ap = T::from_f64(default_ap_f64).unwrap_or_else(|| T::one());
            let ap_u = Field2D::new(grid.nx, grid.ny, default_ap);
            let ap_v = Field2D::new(grid.nx, grid.ny, default_ap);
            rhie_chow.update_u_coefficients(&ap_u);
            rhie_chow.update_v_coefficients(&ap_v);
            Some(rhie_chow)
        } else {
            None
        };

        Ok(Self {
            config,
            grid,
            momentum_solver,
            pressure_solver,
            rhie_chow,
            iterations: 0,
            _u_face_cache: std::cell::RefCell::new(None),
            _v_face_cache: std::cell::RefCell::new(None),
            _d_x_cache: std::cell::RefCell::new(None),
            _d_y_cache: std::cell::RefCell::new(None),
            _vel_field_cache: std::cell::RefCell::new(None),
            _cons_vel_cache: std::cell::RefCell::new(None),
        })
    }

    /// Get algorithm type
    pub fn algorithm(&self) -> AlgorithmType {
        self.config.algorithm
    }

    /// Set boundary condition
    pub fn set_boundary(
        &mut self,
        name: String,
        bc: cfd_core::physics::boundary::BoundaryCondition<T>,
    ) {
        self.momentum_solver.set_boundary(name, bc);
    }

    /// Get current iteration count
    pub fn iterations(&self) -> usize {
        self.iterations
    }

    /// Reset iteration counter
    pub fn reset_iterations(&mut self) {
        self.iterations = 0;
    }
}
