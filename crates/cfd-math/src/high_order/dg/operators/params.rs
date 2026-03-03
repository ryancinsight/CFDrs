//! Parameter types, enums, and builder API for DG operators.

use super::super::{FluxParams, FluxType, LimiterParams, LimiterType};

/// Type of boundary condition
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundaryCondition {
    /// Dirichlet boundary condition (fixed value)
    Dirichlet(f64),
    /// Neumann boundary condition (fixed derivative)
    Neumann(f64),
    /// Periodic boundary condition
    Periodic,
    /// Outflow boundary condition
    Outflow,
    /// Reflective boundary condition
    Reflective,
}

/// Type of numerical flux for boundary conditions
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryFlux {
    /// Central flux
    Central,
    /// Upwind flux
    Upwind,
    /// Lax-Friedrichs flux
    LaxFriedrichs,
    /// HLL flux
    HLL,
}

/// Parameters for DG operators
#[derive(Debug, Clone)]
pub struct DGOperatorParams {
    /// Type of numerical flux for volume integrals
    pub volume_flux: FluxType,
    /// Type of numerical flux for surface integrals
    pub surface_flux: FluxType,
    /// Type of limiter
    pub limiter: LimiterType,
    /// Parameters for the limiter
    pub limiter_params: LimiterParams,
    /// Parameters for the flux
    pub flux_params: FluxParams,
    /// Whether to use strong form (true) or weak form (false)
    pub strong_form: bool,
    /// Whether to use the divergence form (true) or the non-conservative form (false)
    pub divergence_form: bool,
    /// Whether to use the local Lax-Friedrichs flux
    pub use_lax_friedrichs: bool,
    /// Lax-Friedrichs parameter
    pub alpha: f64,
    /// CFL number for time stepping
    pub cfl: f64,
    /// Tolerance for time stepping
    pub tolerance: f64,
    /// Maximum number of iterations
    pub max_iter: usize,
}

impl Default for DGOperatorParams {
    fn default() -> Self {
        Self {
            volume_flux: FluxType::Central,
            surface_flux: FluxType::LaxFriedrichs,
            limiter: LimiterType::Minmod,
            limiter_params: LimiterParams::default(),
            flux_params: FluxParams::default(),
            strong_form: true,
            divergence_form: true,
            use_lax_friedrichs: true,
            alpha: 1.0,
            cfl: 0.1,
            tolerance: 1e-10,
            max_iter: 1000,
        }
    }
}

impl DGOperatorParams {
    /// Create a new set of DG operator parameters
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the volume flux type
    pub fn with_volume_flux(mut self, flux_type: FluxType) -> Self {
        self.volume_flux = flux_type;
        self
    }

    /// Set the surface flux type
    pub fn with_surface_flux(mut self, flux_type: FluxType) -> Self {
        self.surface_flux = flux_type;
        self
    }

    /// Set the limiter type
    pub fn with_limiter(mut self, limiter: LimiterType) -> Self {
        self.limiter = limiter;
        self
    }

    /// Set the limiter parameters
    pub fn with_limiter_params(mut self, params: LimiterParams) -> Self {
        self.limiter_params = params;
        self
    }

    /// Set the flux parameters
    pub fn with_flux_params(mut self, params: FluxParams) -> Self {
        self.flux_params = params;
        self
    }

    /// Set whether to use the strong form
    pub fn with_strong_form(mut self, strong_form: bool) -> Self {
        self.strong_form = strong_form;
        self
    }

    /// Set whether to use the divergence form
    pub fn with_divergence_form(mut self, divergence_form: bool) -> Self {
        self.divergence_form = divergence_form;
        self
    }

    /// Set whether to use the Lax-Friedrichs flux
    pub fn with_lax_friedrichs(mut self, use_lax_friedrichs: bool) -> Self {
        self.use_lax_friedrichs = use_lax_friedrichs;
        self
    }

    /// Set the Lax-Friedrichs parameter
    pub fn with_alpha(mut self, alpha: f64) -> Self {
        self.alpha = alpha;
        self
    }

    /// Set the CFL number
    pub fn with_cfl(mut self, cfl: f64) -> Self {
        self.cfl = cfl;
        self
    }

    /// Set the tolerance
    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance;
        self
    }

    /// Set the maximum number of iterations
    pub fn with_max_iter(mut self, max_iter: usize) -> Self {
        self.max_iter = max_iter;
        self
    }
}
