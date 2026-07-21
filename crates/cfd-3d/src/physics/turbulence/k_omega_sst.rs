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
use eunomia::{FloatElement, NumericElement};

use super::constants::{SST_A1, SST_BETA_STAR};
use super::field_ops::{strain_magnitude, velocity_gradient_tensor};

/// State fields for the k-omega SST model.
#[derive(Debug, Clone)]
pub struct KOmegaSSTState<T: cfd_mesh::domain::core::Scalar + FloatElement> {
    /// Turbulent kinetic energy k [m^2/s^2].
    pub k: Vec<T>,
    /// Specific dissipation rate omega [1/s].
    pub omega: Vec<T>,
    /// Per-point wall distance d \[m].
    pub wall_distance: Vec<T>,
}

/// k-omega SST turbulence model (Menter 1994).
///
/// Blends k-omega near walls and k-epsilon in free-stream.  The eddy viscosity
/// limiter nu_t = a1*k/max(a1*omega, S*F2) suppresses overprediction in APG.
#[derive(Debug, Clone)]
pub struct KOmegaSSTModel<T: cfd_mesh::domain::core::Scalar + FloatElement> {
    /// a1 = 0.31 (eddy viscosity limiter constant).
    pub a1: T,
    /// beta_star = 0.09 (k-equation destruction constant).
    pub beta_star: T,
    /// Kinematic viscosity nu [m^2/s] (used for F1, F2 blending).
    pub nu: T,
    /// Physical grid spacing in the x direction \[m].
    pub dx: T,
    /// Physical grid spacing in the y direction \[m].
    pub dy: T,
    /// Physical grid spacing in the z direction \[m].
    pub dz: T,
    /// Current model state (k, omega, wall distance).
    pub state: Option<KOmegaSSTState<T>>,
}

impl<T: cfd_mesh::domain::core::Scalar + FloatElement> KOmegaSSTModel<T> {
    /// Create a k-omega SST model with given kinematic viscosity.
    pub fn new(nu: T) -> Self {
        Self {
            a1: <T as FloatElement>::from_f64(SST_A1),
            beta_star: <T as FloatElement>::from_f64(SST_BETA_STAR),
            nu,
            dx: T::ONE,
            dy: T::ONE,
            dz: T::ONE,
            state: None,
        }
    }

    /// Create a k-omega SST model with the physical grid spacing.
    pub fn with_grid_spacing(nu: T, dx: T, dy: T, dz: T) -> Self {
        Self {
            a1: <T as FloatElement>::from_f64(SST_A1),
            beta_star: <T as FloatElement>::from_f64(SST_BETA_STAR),
            nu,
            dx,
            dy,
            dz,
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
        assert_eq!(
            wall_distances.len(),
            n,
            "SST wall_distances must match the flow-field size"
        );
        let three_half = <T as FloatElement>::from_f64(1.5);
        let four = T::ONE + T::ONE + T::ONE + T::ONE;
        let c_mu = <T as FloatElement>::from_f64(
            cfd_core::physics::constants::physics::turbulence::K_EPSILON_C_MU,
        );
        let c_mu_fourth_root = <T as FloatElement>::powf(c_mu, T::ONE / four);

        let k_field: Vec<T> = flow_field
            .velocity
            .components
            .iter()
            .map(|v| {
                let u_mag = <T as NumericElement>::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
                three_half * <T as FloatElement>::powi(u_mag * turbulence_intensity, 2)
            })
            .collect();

        let omega_field: Vec<T> = k_field
            .iter()
            .map(|&k| <T as NumericElement>::sqrt(k) / (c_mu_fourth_root * length_scale))
            .collect();

        self.state = Some(KOmegaSSTState {
            k: k_field,
            omega: omega_field,
            wall_distance: wall_distances,
        });
    }

    /// Initialise with exact k and omega fields.
    pub fn initialize_exact(&mut self, k: Vec<T>, omega: Vec<T>, wall_distances: Vec<T>) {
        assert_eq!(k.len(), omega.len(), "SST k and omega lengths must match");
        assert_eq!(
            k.len(),
            wall_distances.len(),
            "SST wall_distance length must match k and omega"
        );
        self.state = Some(KOmegaSSTState {
            k,
            omega,
            wall_distance: wall_distances,
        });
    }

    /// Compute strain rate magnitude |S| at grid point (i, j, k).
    fn strain_rate(&self, flow: &FlowField<T>, i: usize, j: usize, k: usize) -> T {
        let (nx, ny, nz) = flow.velocity.dimensions;
        let gradient = velocity_gradient_tensor(
            &flow.velocity.components,
            nx,
            ny,
            nz,
            i,
            j,
            k,
            self.dx,
            self.dy,
            self.dz,
        );
        strain_magnitude(&gradient)
    }
}

impl<T: cfd_mesh::domain::core::Scalar + FloatElement> TurbulenceModel<T> for KOmegaSSTModel<T> {
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let n = nx * ny * nz;
        let eps = <T as FloatElement>::from_f64(1e-15);
        match &self.state {
            None => vec![T::ZERO; n],
            Some(state) => {
                assert_eq!(state.k.len(), n, "SST k must match the flow-field size");
                assert_eq!(
                    state.omega.len(),
                    n,
                    "SST omega must match the flow-field size"
                );
                assert_eq!(
                    state.wall_distance.len(),
                    n,
                    "SST wall_distance must match the flow-field size"
                );
                let mut viscosity = Vec::with_capacity(n);
                for idx in 0..n {
                    let k = state.k[idx];
                    let om = state.omega[idx];
                    let d = state.wall_distance[idx];
                    let arg2 = <T as NumericElement>::max_scalar(
                        <T as FloatElement>::from_f64(2.0) * <T as NumericElement>::sqrt(k)
                            / (self.beta_star * om * d + eps),
                        <T as FloatElement>::from_f64(500.0) * self.nu / (d * d * om + eps),
                    );
                    let f2 = <T as FloatElement>::tanh(arg2 * arg2);
                    let ii = idx % nx;
                    let jj = (idx / nx) % ny;
                    let kk = idx / (nx * ny);
                    let s_mag = self.strain_rate(flow_field, ii, jj, kk);
                    let denom = <T as NumericElement>::max_scalar(self.a1 * om, s_mag * f2 + eps);
                    viscosity.push(self.a1 * k / denom);
                }
                viscosity
            }
        }
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let n = flow_field.velocity.components.len();
        match &self.state {
            None => vec![T::ZERO; n],
            Some(state) => {
                assert_eq!(state.k.len(), n, "SST k must match the flow-field size");
                state.k.clone()
            }
        }
    }

    fn name(&self) -> &'static str {
        "k-omega-SST"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;
    use cfd_core::physics::fluid_dynamics::TurbulenceModel;
    use leto::geometry::Vector3;

    fn fill_velocity_field<F>(flow: &mut FlowField<f64>, mut generator: F)
    where
        F: FnMut(f64, f64, f64) -> Vector3<f64>,
    {
        let (nx, ny, nz) = flow.velocity.dimensions;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = k * nx * ny + j * nx + i;
                    flow.velocity.components[idx] = generator(i as f64, j as f64, k as f64);
                }
            }
        }
    }

    #[test]
    fn simple_shear_uses_full_strain_tensor() {
        let dx = 1.0;
        let dy = 2.0;
        let dz = 3.0;
        let nu = 1.0e-6;
        let mut flow = FlowField::<f64>::new(3, 3, 3);
        fill_velocity_field(&mut flow, |_, j, _| Vector3::new(j * dy, 0.0, 0.0));

        let n = flow.velocity.components.len();
        let mut model = KOmegaSSTModel::<f64>::with_grid_spacing(nu, dx, dy, dz);
        model.initialize_exact(vec![1.0; n], vec![1.0; n], vec![1.0; n]);

        let viscosity = model.turbulent_viscosity(&flow);
        let eps = 1.0e-15;
        let arg2 = (2.0 * 1.0_f64.sqrt() / (0.09 * 1.0 * 1.0 + eps)).max(500.0 * nu / (1.0 + eps));
        let f2 = (arg2 * arg2).tanh();
        let expected = 0.31_f64 * 1.0_f64 / (0.31_f64 * 1.0_f64).max(f2 + eps);

        assert_relative_eq!(viscosity[13], expected, epsilon = 1e-12);
        assert!(viscosity[13] < 1.0);
    }

    #[test]
    fn initialize_state_uses_cmu_fourth_root_in_omega() {
        let mut flow = FlowField::<f64>::new(1, 1, 1);
        fill_velocity_field(&mut flow, |_, _, _| Vector3::new(4.0, 0.0, 0.0));

        let nu = 1.0e-6;
        let turbulence_intensity = 0.1;
        let length_scale = 0.2;
        let mut model = KOmegaSSTModel::<f64>::with_grid_spacing(nu, 1.0, 1.0, 1.0);
        model.initialize_state(&flow, turbulence_intensity, length_scale, vec![0.5]);

        let state = model
            .state
            .as_ref()
            .expect("initialize_state must populate SST state");
        let expected_k = 1.5 * (4.0 * turbulence_intensity).powi(2);
        let c_mu = cfd_core::physics::constants::physics::turbulence::K_EPSILON_C_MU;
        let expected_omega = expected_k.sqrt() / (c_mu.powf(0.25) * length_scale);

        assert_relative_eq!(state.k[0], expected_k, epsilon = 1e-12);
        assert_relative_eq!(state.omega[0], expected_omega, epsilon = 1e-12);
    }
}
