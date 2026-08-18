use crate::grid::array2d::Array2D;
use crate::scalar;
use crate::solvers::ns_fvm::boundary::BloodModel;
use crate::solvers::ns_fvm::config::SIMPLEConfig;
use crate::solvers::ns_fvm::field::FlowField2D;
use crate::solvers::ns_fvm::grid::StaggeredGrid2D;
use cfd_core::CfdScalar;
use eunomia::FloatElement;

#[derive(Clone, Copy)]
pub(super) enum InletProfile {
    Parabolic,
}

/// 2D Navier-Stokes FVM solver with SIMPLE pressure-velocity coupling.
///
/// Used as the numerical engine by geometry-specific pass-through solvers
/// (bifurcation, Venturi, serpentine, etc.).
pub struct NavierStokesSolver2D<T: CfdScalar + Copy + FloatElement> {
    /// Computational grid
    pub grid: StaggeredGrid2D<T>,
    /// Flow field (u, v, p, μ, γ̇)
    pub field: FlowField2D<T>,
    /// Blood rheology model
    pub blood: BloodModel<T>,
    /// Fluid density [kg/m³]
    pub density: T,
    /// SIMPLE configuration
    pub config: SIMPLEConfig<T>,
    /// Central coefficient storage for u-momentum
    pub(super) a_p_u: Array2D<T>,
    /// Central coefficient storage for v-momentum
    pub(super) a_p_v: Array2D<T>,
    /// Snapshot of the staggered u field reused across momentum sweeps.
    pub(super) u_old_workspace: Array2D<T>,
    /// Snapshot of the staggered v field reused across momentum sweeps.
    pub(super) v_old_workspace: Array2D<T>,
    /// Pressure Poisson coefficient workspace for east-west face diffusion.
    pub(super) pressure_poisson_d_u: Array2D<T>,
    /// Pressure Poisson coefficient workspace for north-south face diffusion.
    pub(super) pressure_poisson_d_v: Array2D<T>,
    /// Pressure correction workspace reused across SIMPLE iterations.
    pub(super) pressure_poisson_p_prime: Array2D<T>,
    /// Pressure correction right-hand-side workspace reused across SIMPLE iterations.
    pub(super) pressure_poisson_rhs: Array2D<T>,
    /// East pressure-correction coefficient workspace reused across SOR sweeps.
    pub(super) pressure_poisson_a_e: Array2D<T>,
    /// West pressure-correction coefficient workspace reused across SOR sweeps.
    pub(super) pressure_poisson_a_w: Array2D<T>,
    /// North pressure-correction coefficient workspace reused across SOR sweeps.
    pub(super) pressure_poisson_a_n: Array2D<T>,
    /// South pressure-correction coefficient workspace reused across SOR sweeps.
    pub(super) pressure_poisson_a_s: Array2D<T>,
    /// Diagonal pressure-correction coefficient workspace reused across SOR sweeps.
    pub(super) pressure_poisson_a_p: Array2D<T>,
    /// Inlet profile selected by the geometry-specific solve wrapper.
    pub(super) inlet_profile: InletProfile,
    /// Optional k-omega SST turbulence model.  When Some, the solver
    /// computes turbulent viscosity nu_t each iteration and adds it to
    /// the molecular viscosity in the momentum equation diffusion terms.
    pub(super) turbulence: Option<TurbulenceCoupling<T>>,
}

/// Turbulence model coupling state for the SIMPLE solver.
pub(super) struct TurbulenceCoupling<T: CfdScalar + Copy> {
    /// k-omega SST model.
    pub(super) model: crate::physics::turbulence::k_omega_sst::KOmegaSSTModel<T>,
    /// Turbulent kinetic energy at cell centers \[nx]\[ny].
    pub(super) k: Vec<T>,
    /// Specific dissipation rate at cell centers \[nx]\[ny].
    pub(super) omega: Vec<T>,
    /// Update interval (every N SIMPLE iterations).
    pub(super) update_interval: usize,
}

impl<T: CfdScalar + Copy + FloatElement> NavierStokesSolver2D<T> {
    /// Create a new solver.
    pub fn new(
        grid: StaggeredGrid2D<T>,
        blood: BloodModel<T>,
        density: T,
        config: SIMPLEConfig<T>,
    ) -> Self {
        let field = FlowField2D::<T>::new(grid.nx, grid.ny);
        let one: T = scalar::one();
        let zero: T = scalar::zero();
        let a_p_u = Array2D::new(grid.nx + 1, grid.ny, one);
        let a_p_v = Array2D::new(grid.nx, grid.ny + 1, one);
        let u_old_workspace = Array2D::new(grid.nx + 1, grid.ny, zero);
        let v_old_workspace = Array2D::new(grid.nx, grid.ny + 1, zero);
        let pressure_poisson_d_u = Array2D::new(grid.nx + 1, grid.ny, zero);
        let pressure_poisson_d_v = Array2D::new(grid.nx, grid.ny + 1, zero);
        let pressure_poisson_p_prime = Array2D::new(grid.nx, grid.ny, zero);
        let pressure_poisson_rhs = Array2D::new(grid.nx, grid.ny, zero);
        let pressure_poisson_a_e = Array2D::new(grid.nx, grid.ny, zero);
        let pressure_poisson_a_w = Array2D::new(grid.nx, grid.ny, zero);
        let pressure_poisson_a_n = Array2D::new(grid.nx, grid.ny, zero);
        let pressure_poisson_a_s = Array2D::new(grid.nx, grid.ny, zero);
        let pressure_poisson_a_p = Array2D::new(grid.nx, grid.ny, zero);
        Self {
            grid,
            field,
            blood,
            density,
            config,
            a_p_u,
            a_p_v,
            u_old_workspace,
            v_old_workspace,
            pressure_poisson_d_u,
            pressure_poisson_d_v,
            pressure_poisson_p_prime,
            pressure_poisson_rhs,
            pressure_poisson_a_e,
            pressure_poisson_a_w,
            pressure_poisson_a_n,
            pressure_poisson_a_s,
            pressure_poisson_a_p,
            inlet_profile: InletProfile::Uniform,
            turbulence: None,
        }
    }

    /// Enable k-omega SST turbulence modeling for high-Re flows.
    ///
    /// When enabled, the solver computes turbulent viscosity nu_t at
    /// each iteration and adds it to the molecular viscosity in the
    /// momentum equation.  This is needed for venturi throat flows
    /// at Re > 2000 where the laminar assumption breaks down.
    pub fn enable_turbulence(&mut self) {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let size = nx * ny;
        // Initial k and omega from free-stream turbulence intensity ~1%.
        let u_ref: T = <T as FloatElement>::from_f64(0.1);
        let ti: T = <T as FloatElement>::from_f64(0.01);
        let k_init = <T as FloatElement>::from_f64(1.5) * (u_ref * ti) * (u_ref * ti);
        let omega_init = k_init / <T as FloatElement>::from_f64(0.001);
        self.turbulence = Some(TurbulenceCoupling {
            model: crate::physics::turbulence::k_omega_sst::KOmegaSSTModel::new(nx, ny),
            k: vec![k_init; size],
            omega: vec![omega_init; size],
            update_interval: 5,
        });
    }

    /// Initialise viscosity field from the blood model's apparent viscosity at
    /// the reference shear rate.  Using μ_∞ (high-shear limit) here would
    /// under-estimate viscosity by ~30 %, making the initial velocity field too
    /// fast and slowing SIMPLE convergence for non-Newtonian blood.
    pub fn initialize_viscosity(&mut self) {
        let mu_init = match &self.blood {
            BloodModel::Casson(m) => m.apparent_viscosity(m.reference_shear_rate.into_base()),
            BloodModel::CarreauYasuda(m) => {
                m.apparent_viscosity(m.reference_shear_rate.into_base())
            }
            BloodModel::Newtonian(mu) => *mu,
        };
        for i in 0..self.grid.nx {
            for j in 0..self.grid.ny {
                self.field.mu[(i, j)] = mu_init;
            }
        }
    }
}
