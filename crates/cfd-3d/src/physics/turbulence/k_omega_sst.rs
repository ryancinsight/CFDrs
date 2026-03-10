//! k-omega SST (Shear Stress Transport) turbulence model (Menter 1994).
//!
//! # Theorem -- k-omega SST Blending (Menter 1994)
//!
//! Blends k-omega (wall layer) with k-epsilon (free stream) via F1 = tanh(arg1^4).
//! Eddy viscosity limiter: nu_t = a1*k / max(a1*omega, S*F2).
//!
//! ## References
//!
//! - Menter, F.R. (1994). Two-equation eddy-viscosity turbulence models.
//!   AIAA J. 32(8):1598-1605.
//! - Menter, Kuntz & Langtry (2003). Ten years of industrial experience.
//!   Turbulence, Heat and Mass Transfer 4, Begell House.

use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;
use nalgebra::RealField;
use num_traits::FromPrimitive;

use super::constants::{SST_A1, SST_BETA_STAR};

/// State fields for the k-omega SST model.
#[derive(Debug, Clone)]
pub struct KOmegaSSTState<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Turbulent kinetic energy k [m^2/s^2].
    pub k: Vec<T>,
    /// Specific dissipation rate omega [1/s].
    pub omega: Vec<T>,
    /// Per-point wall distance d [m].
    pub wall_distance: Vec<T>,
}

/// k-omega SST turbulence model (Menter 1994).
///
/// Blends k-omega near walls and k-epsilon in free-stream.  The eddy viscosity
/// limiter nu_t = a1*k/max(a1*omega, S*F2) suppresses overprediction in APG.
#[derive(Debug, Clone)]
pub struct KOmegaSSTModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// a1 = 0.31 (eddy viscosity limiter constant).
    pub a1: T,
    /// beta_star = 0.09 (k-equation destruction constant).
    pub beta_star: T,
    /// Kinematic viscosity nu [m^2/s] (used for F1, F2 blending).
    pub nu: T,
    /// Current model state (k, omega, wall distance).
    pub state: Option<KOmegaSSTState<T>>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> KOmegaSSTModel<T> {
    /// Create a k-omega SST model with given kinematic viscosity.
    pub fn new(nu: T) -> Self {
        Self {
            a1: <T as FromPrimitive>::from_f64(SST_A1)
                .expect("SST_A1 is an IEEE 754 representable f64 constant"),
            beta_star: <T as FromPrimitive>::from_f64(SST_BETA_STAR)
                .expect("SST_BETA_STAR is an IEEE 754 representable f64 constant"),
            nu,
            state: None,
        }
    }

    /// Initialise the SST state with turbulence intensity and length scale.
    pub fn initialize_state(
        &mut self,
        flow_field: &FlowField<T>,
        turbulence_intensity: T,
        length_scale: T,
        wall_distances: Vec<T>,
    ) {
        let n = flow_field.velocity.components.len();
        let three_half = <T as FromPrimitive>::from_f64(1.5)
            .expect("1.5 is an IEEE 754 representable f64 constant");
        let c_mu_quarter = <T as FromPrimitive>::from_f64(0.09_f64.sqrt())
            .expect("sqrt(0.09) is an IEEE 754 representable f64 constant");

        let k_field: Vec<T> = flow_field
            .velocity
            .components
            .iter()
            .map(|v| {
                let u_mag = v.norm();
                three_half * num_traits::Float::powi(u_mag * turbulence_intensity, 2)
            })
            .collect();

        let omega_field: Vec<T> = k_field
            .iter()
            .map(|&k| num_traits::Float::sqrt(k) / (c_mu_quarter * length_scale))
            .collect();

        let wall_dist = if wall_distances.len() == n {
            wall_distances
        } else {
            vec![T::one(); n]
        };

        self.state = Some(KOmegaSSTState {
            k: k_field,
            omega: omega_field,
            wall_distance: wall_dist,
        });
    }

    /// Initialise with exact k and omega fields.
    pub fn initialize_exact(&mut self, k: Vec<T>, omega: Vec<T>, wall_distances: Vec<T>) {
        self.state = Some(KOmegaSSTState {
            k,
            omega,
            wall_distance: wall_distances,
        });
    }

    /// Compute strain rate magnitude |S| at grid point (i, j, k).
    fn strain_rate(&self, flow: &FlowField<T>, i: usize, j: usize, k: usize) -> T {
        let (nx, ny, nz) = flow.velocity.dimensions;
        let two = <T as FromPrimitive>::from_f64(2.0)
            .expect("2.0 is representable in all IEEE 754 types");
        let h = T::one();
        let mut s_sq = T::zero();
        if i > 0 && i < nx - 1 {
            if let (Some(vp), Some(vm)) = (
                flow.velocity.get(i + 1, j, k),
                flow.velocity.get(i - 1, j, k),
            ) {
                let s = (vp.x - vm.x) / (two * h);
                s_sq += s * s;
            }
        }
        if j > 0 && j < ny - 1 {
            if let (Some(vp), Some(vm)) = (
                flow.velocity.get(i, j + 1, k),
                flow.velocity.get(i, j - 1, k),
            ) {
                let s = (vp.y - vm.y) / (two * h);
                s_sq += s * s;
            }
        }
        if k > 0 && k < nz - 1 {
            if let (Some(vp), Some(vm)) = (
                flow.velocity.get(i, j, k + 1),
                flow.velocity.get(i, j, k - 1),
            ) {
                let s = (vp.z - vm.z) / (two * h);
                s_sq += s * s;
            }
        }
        num_traits::Float::sqrt(two * s_sq)
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> TurbulenceModel<T>
    for KOmegaSSTModel<T>
{
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let n = nx * ny * nz;
        let eps = <T as FromPrimitive>::from_f64(1e-15)
            .expect("1e-15 is an IEEE 754 representable f64 constant");
        match &self.state {
            None => vec![T::zero(); n],
            Some(state) => {
                let mut viscosity = Vec::with_capacity(n);
                for idx in 0..n {
                    let k = if idx < state.k.len() {
                        state.k[idx]
                    } else {
                        T::zero()
                    };
                    let om = if idx < state.omega.len() {
                        state.omega[idx]
                    } else {
                        T::one()
                    };
                    let d = if idx < state.wall_distance.len() {
                        state.wall_distance[idx]
                    } else {
                        T::one()
                    };
                    let arg2 = num_traits::Float::max(
                        <T as FromPrimitive>::from_f64(2.0)
                            .expect("2.0 is representable in all IEEE 754 types")
                            * num_traits::Float::sqrt(k)
                            / (self.beta_star * om * d + eps),
                        <T as FromPrimitive>::from_f64(500.0)
                            .expect("500.0 is representable in all IEEE 754 types")
                            * self.nu
                            / (d * d * om + eps),
                    );
                    let f2 = num_traits::Float::tanh(arg2 * arg2);
                    let ii = idx % nx;
                    let jj = (idx / nx) % ny;
                    let kk = idx / (nx * ny);
                    let s_mag = self.strain_rate(flow_field, ii, jj, kk);
                    let denom = num_traits::Float::max(self.a1 * om, s_mag * f2 + eps);
                    viscosity.push(self.a1 * k / denom);
                }
                viscosity
            }
        }
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let n = flow_field.velocity.components.len();
        match &self.state {
            None => vec![T::zero(); n],
            Some(state) => state.k.clone(),
        }
    }

    fn name(&self) -> &'static str {
        "k-omega-SST"
    }
}
