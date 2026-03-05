//! Prandtl mixing length turbulence model (Prandtl 1925).
//!
//! # Theorem — Prandtl Mixing Length Hypothesis (Prandtl 1925)
//!
//! The turbulent shear stress in wall-bounded flow is modelled as:
//!
//! ```text
//! τ = ρ lₘ² |∂u/∂y|²
//! νₜ = lₘ² |∂u/∂y|
//! ```
//!
//! where $l_m$ is the mixing length. Near a wall, $l_m = \kappa y$
//! ($\kappa = 0.41$, von Kármán constant), transitioning to a constant value
//! in the outer region.
//!
//! **Proof sketch.** Prandtl's analogy with molecular transport postulates
//! turbulent eddies traverse a "mixing length" $l_m$ before fully exchanging
//! momentum. The turbulent velocity fluctuation scale is then $u' \sim l_m
//! |\partial \bar{u}/\partial y|$, giving $\overline{u'v'} \sim l_m^2
//! |\partial \bar{u}/\partial y|^2$, from which $\nu_t = l_m^2 |\partial
//! \bar{u}/\partial y|$ follows by dimensional consistency with the Reynolds
//! stress.
//!
//! **Reference:** Prandtl, L. (1925). "Über die ausgebildete Turbulenz."
//! *Z. Angew. Math. Mech.* 5:136–139.

use cfd_core::physics::fluid_dynamics::fields::{FlowField, VelocityField};
use cfd_core::physics::fluid_dynamics::TurbulenceModel;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Prandtl mixing length turbulence model.
///
/// See [module-level documentation](self) for the theorem and proof sketch.
#[derive(Debug, Clone)]
pub struct MixingLengthModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Mixing length scale lₘ [m]
    pub length_scale: T,
    /// von Kármán constant κ = 0.41 (dimensionless)
    pub kappa: T,
    /// Physical finite-difference step for velocity gradient computation [m].
    ///
    /// For grid-resolved computations, set to the physical cell spacing
    /// Δ = (dx·dy·dz)^(1/3) via [`MixingLengthModel::with_filter_width`].
    /// The default (`new()`) sets this equal to `length_scale`.
    pub filter_width: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> MixingLengthModel<T> {
    /// Create a new mixing length model.
    ///
    /// Sets `filter_width = length_scale` (appropriate when grid spacing and
    /// mixing length are comparable).  For physically correct gradient
    /// computation use [`MixingLengthModel::with_filter_width`].
    pub fn new(length_scale: T) -> Self {
        let kappa = <T as FromPrimitive>::from_f64(
            cfd_core::physics::constants::physics::fluid::VON_KARMAN,
        )
        .expect("VON_KARMAN is an IEEE 754 representable f64 constant");
        Self {
            length_scale,
            kappa,
            filter_width: length_scale,
        }
    }

    /// Create a mixing length model with the physically correct FD step.
    ///
    /// Computes Δ = (dx·dy·dz)^(1/3) per Deardorff (1970) as the finite-
    /// difference spacing used to evaluate velocity gradients.
    ///
    /// # Arguments
    /// * `length_scale` — Prandtl mixing length lₘ [m]
    /// * `dx`, `dy`, `dz` — physical cell dimensions [m]
    pub fn with_filter_width(length_scale: T, dx: T, dy: T, dz: T) -> Self {
        let kappa = <T as FromPrimitive>::from_f64(
            cfd_core::physics::constants::physics::fluid::VON_KARMAN,
        )
        .expect("VON_KARMAN is an IEEE 754 representable f64 constant");
        let one_third = T::one() / (T::one() + T::one() + T::one());
        let filter_width = num_traits::Float::powf(dx * dy * dz, one_third);
        Self {
            length_scale,
            kappa,
            filter_width,
        }
    }

    /// Calculate velocity gradient magnitude using the stored physical filter width.
    ///
    /// Central differences with spacing `self.filter_width` (= physical cell size Δ).
    #[allow(clippy::similar_names)] // CFD derivatives use standard notation
    fn calculate_velocity_gradient(&self, velocity: &VelocityField<T>, idx: usize) -> T {
        let (nx, ny, nz) = velocity.dimensions;
        let i = idx % nx;
        let j = (idx / nx) % ny;
        let k = idx / (nx * ny);

        // Physically correct FD spacing — see filter_width field documentation.
        let delta = self.filter_width;
        let two = T::one() + T::one();

        let mut grad_u_sq = T::zero();

        // Calculate ∂u/∂y for wall-bounded flows (primary gradient)
        if j > 0 && j < ny - 1 {
            if let (Some(u_jp), Some(u_jm)) = (velocity.get(i, j + 1, k), velocity.get(i, j - 1, k))
            {
                let dudy = (u_jp.x - u_jm.x) / (two * delta);
                grad_u_sq = dudy * dudy;
            }
        }

        // Add other gradient components for 3D flows
        if i > 0 && i < nx - 1 {
            if let (Some(u_ip), Some(u_im)) = (velocity.get(i + 1, j, k), velocity.get(i - 1, j, k))
            {
                #[allow(clippy::similar_names)] // CFD derivatives use standard notation
                {
                    let dudx = (u_ip.x - u_im.x) / (two * delta);
                    let dvdx = (u_ip.y - u_im.y) / (two * delta);
                    let dwdx = (u_ip.z - u_im.z) / (two * delta);
                    grad_u_sq = grad_u_sq + dudx * dudx + dvdx * dvdx + dwdx * dwdx;
                }
            }
        }

        if k > 0 && k < nz - 1 {
            if let (Some(u_kp), Some(u_km)) = (velocity.get(i, j, k + 1), velocity.get(i, j, k - 1))
            {
                #[allow(clippy::similar_names)] // CFD derivatives use standard notation
                {
                    let dudz = (u_kp.x - u_km.x) / (two * delta);
                    let dvdz = (u_kp.y - u_km.y) / (two * delta);
                    let dwdz = (u_kp.z - u_km.z) / (two * delta);
                    grad_u_sq = grad_u_sq + dudz * dudz + dvdz * dvdz + dwdz * dwdz;
                }
            }
        }

        num_traits::Float::sqrt(grad_u_sq)
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> TurbulenceModel<T>
    for MixingLengthModel<T>
{
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        // νₜ = l² * |∂u/∂y| (Prandtl's mixing length hypothesis)
        flow_field
            .velocity
            .components
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let grad_u = self.calculate_velocity_gradient(&flow_field.velocity, idx);
                self.length_scale * self.length_scale * grad_u
            })
            .collect()
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        // Estimate TKE from mixing length and velocity gradient
        // k ≈ (l * |∂u/∂y|)²
        flow_field
            .velocity
            .components
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let grad_u = self.calculate_velocity_gradient(&flow_field.velocity, idx);
                let l_grad_u = self.length_scale * grad_u;
                l_grad_u * l_grad_u
            })
            .collect()
    }

    fn name(&self) -> &'static str {
        "MixingLength"
    }
}
