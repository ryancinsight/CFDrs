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
use crate::grid::array2d::Array2D;
use crate::grid::StructuredGrid2D;
use crate::physics::MomentumSolver;
use crate::pressure_velocity::{PressureCorrectionSolver, RhieChowInterpolation};
use crate::scalar;
use cfd_core::error::Error;
use cfd_core::CfdScalar;
use eunomia::{FloatElement, NumericElement};
use leto::geometry::Vector2;

type PressureLayer = Box<[(usize, usize)]>;

struct SolidPressureLayers {
    mask: Box<[bool]>,
    layers: Box<[PressureLayer]>,
}

/// SIMPLEC/PIMPLE pressure-velocity coupling solver
///
/// This solver extends the basic SIMPLE algorithm with:
/// - SIMPLEC: Consistent pressure-velocity coupling using Rhie-Chow interpolation
/// - PIMPLE: Merged PISO-SIMPLE for better convergence in transient flows
pub struct SimplecPimpleSolver<T: CfdScalar + Copy> {
    pub(super) config: SimplecPimpleConfig<T>,
    pub(super) grid: StructuredGrid2D<T>,
    pub(super) momentum_solver: MomentumSolver<T>,
    pub(super) pressure_solver: PressureCorrectionSolver<T>,
    pub(super) rhie_chow: Option<RhieChowInterpolation<T>>,
    pub(super) iterations: usize,
    pub(super) u_star_workspace: Array2D<Vector2<T>>,
    pub(super) u_corrected_workspace: Array2D<Vector2<T>>,
    pub(super) p_workspace: Array2D<T>,
    pub(super) p_correction_workspace: Array2D<T>,
    pub(super) _u_face_cache: std::cell::RefCell<Option<crate::grid::array2d::Array2D<T>>>,
    pub(super) _v_face_cache: std::cell::RefCell<Option<crate::grid::array2d::Array2D<T>>>,
    pub(super) _d_x_cache: std::cell::RefCell<Option<crate::grid::array2d::Array2D<T>>>,
    pub(super) _d_y_cache: std::cell::RefCell<Option<crate::grid::array2d::Array2D<T>>>,
    pub(super) _vel_field_cache: std::cell::RefCell<Option<crate::fields::Field2D<Vector2<T>>>>,
    pub(super) _cons_vel_cache:
        std::cell::RefCell<Option<crate::grid::array2d::Array2D<Vector2<T>>>>,
    _solid_pressure_layers: std::cell::RefCell<Option<SolidPressureLayers>>,
}

impl<T: CfdScalar + Copy + std::fmt::LowerExp + FloatElement> SimplecPimpleSolver<T> {
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
            crate::schemes::SpatialScheme::SecondOrderUpwind => {
                crate::physics::momentum::ConvectionScheme::SecondOrderUpwind {
                    relaxation_factor: 0.7,
                }
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
            crate::schemes::SpatialScheme::WenoZ5 => {
                crate::physics::momentum::ConvectionScheme::WenoZ {
                    relaxation_factor: 0.7,
                }
            }
            crate::schemes::SpatialScheme::CentralDifference
            | crate::schemes::SpatialScheme::Weno5
            | crate::schemes::SpatialScheme::Weno9
            | crate::schemes::SpatialScheme::FourthOrderCentral => {
                return Err(Error::InvalidConfiguration(
                    "Selected convection scheme is not supported by SIMPLEC/PIMPLE".into(),
                ));
            }
        };

        let mut momentum_solver = MomentumSolver::new(&grid);
        momentum_solver.set_convection_scheme(momentum_convection);
        momentum_solver.set_velocity_relaxation(config.alpha_u);
        // Eisenstat-Walker-style inexact coupling: the linear momentum solve
        // is one order tighter than the SIMPLE nonlinear residual contract.
        // This prevents the default 1e-8 linear target from over-solving a
        // configuration whose outer tolerance is materially looser.
        momentum_solver
            .set_linear_solver_tolerance(config.tolerance * <T as FloatElement>::from_f64(0.1));

        // Create pressure solver using configuration
        let pressure_solver =
            PressureCorrectionSolver::new(grid.clone(), config.pressure_linear_solver)?;
        let grid_nx = grid.nx;
        let grid_ny = grid.ny;

        let rhie_chow = if config.use_rhie_chow {
            let mut rhie_chow = RhieChowInterpolation::new(&grid);
            // Initialize Rhie-Chow coefficients with reasonable defaults
            // This prevents issues in the first iteration
            let dx_f64 = NumericElement::to_f64(grid.dx);
            let dy_f64 = NumericElement::to_f64(grid.dy);
            let default_ap_f64 = 1.0 / (dx_f64 * dx_f64 + dy_f64 * dy_f64);
            let default_ap = <T as FloatElement>::from_f64(default_ap_f64);
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
            u_star_workspace: Array2D::new(grid_nx, grid_ny, Vector2::zeros()),
            u_corrected_workspace: Array2D::new(grid_nx, grid_ny, Vector2::zeros()),
            p_workspace: Array2D::new(grid_nx, grid_ny, scalar::zero()),
            p_correction_workspace: Array2D::new(grid_nx, grid_ny, scalar::zero()),
            _u_face_cache: std::cell::RefCell::new(None),
            _v_face_cache: std::cell::RefCell::new(None),
            _d_x_cache: std::cell::RefCell::new(None),
            _d_y_cache: std::cell::RefCell::new(None),
            _vel_field_cache: std::cell::RefCell::new(None),
            _cons_vel_cache: std::cell::RefCell::new(None),
            _solid_pressure_layers: std::cell::RefCell::new(None),
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

    /// Extrapolate pressure field to solid cells
    pub(super) fn extrapolate_pressure_to_solids(
        &self,
        fields: &mut crate::fields::SimulationFields<T>,
    ) {
        let nx = self.grid.nx;
        let ny = self.grid.ny;

        let mask = fields.mask.as_slice();
        let mut layers_cache = self._solid_pressure_layers.borrow_mut();
        if layers_cache
            .as_ref()
            .is_none_or(|cached| cached.mask.as_ref() != mask)
        {
            let mut valid = mask.to_vec();
            let mut queued = vec![false; nx * ny];
            let mut frontier = Vec::new();
            for j in 0..ny {
                for i in 0..nx {
                    if valid[i * ny + j] {
                        frontier.push((i, j));
                    }
                }
            }

            let mut layers = Vec::new();
            // Discover synchronous distance layers once per mask. A layer is
            // updated from the fluid cells and every earlier layer, matching
            // the former full-grid sweep exactly.
            while !frontier.is_empty() {
                let mut candidates = Vec::new();
                for &(i, j) in &frontier {
                    let neighbors = [
                        (i.checked_sub(1), Some(j)),
                        (i.checked_add(1).filter(|&next_i| next_i < nx), Some(j)),
                        (Some(i), j.checked_sub(1)),
                        (Some(i), j.checked_add(1).filter(|&next_j| next_j < ny)),
                    ];
                    for (maybe_i, maybe_j) in neighbors {
                        let (Some(candidate_i), Some(candidate_j)) = (maybe_i, maybe_j) else {
                            continue;
                        };
                        let candidate_idx = candidate_i * ny + candidate_j;
                        if !valid[candidate_idx] && !queued[candidate_idx] {
                            queued[candidate_idx] = true;
                            candidates.push((candidate_i, candidate_j));
                        }
                    }
                }

                if candidates.is_empty() {
                    break;
                }

                candidates.sort_unstable_by_key(|&(i, j)| j * nx + i);
                layers.push(candidates.clone().into_boxed_slice());
                for &(i, j) in &candidates {
                    valid[i * ny + j] = true;
                }
                frontier = candidates;
            }

            *layers_cache = Some(SolidPressureLayers {
                mask: mask.to_vec().into_boxed_slice(),
                layers: layers.into_boxed_slice(),
            });
        }

        let cached = layers_cache
            .as_ref()
            .expect("invariant: pressure layers are initialized above");
        let mut valid = mask.to_vec();
        for layer in &cached.layers {
            let mut updates = Vec::with_capacity(layer.len());
            for &(i, j) in layer {
                let mut sum: T = scalar::zero();
                let mut count = 0;

                if i > 0 && valid[(i - 1) * ny + j] {
                    sum += fields.p.at(i - 1, j);
                    count += 1;
                }
                if i < nx - 1 && valid[(i + 1) * ny + j] {
                    sum += fields.p.at(i + 1, j);
                    count += 1;
                }
                if j > 0 && valid[i * ny + j - 1] {
                    sum += fields.p.at(i, j - 1);
                    count += 1;
                }
                if j < ny - 1 && valid[i * ny + j + 1] {
                    sum += fields.p.at(i, j + 1);
                    count += 1;
                }

                if count > 0 {
                    updates.push(((i, j), sum / scalar::from_usize(count)));
                }
            }

            for ((i, j), value) in updates {
                fields.p.set(i, j, value);
                valid[i * ny + j] = true;
            }
        }
    }
}
