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

use super::field_ops::derivative_y;

/// Prandtl mixing length turbulence model.
///
/// See [module-level documentation](self) for the theorem and proof sketch.
#[derive(Debug, Clone)]
pub struct MixingLengthModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Mixing length scale lₘ [m]
    pub length_scale: T,
    /// von Kármán constant κ = 0.41 (dimensionless)
    pub kappa: T,
    /// Wall-normal finite-difference step for velocity gradient computation [m].
    ///
    /// For grid-resolved computations, set this to the physical wall-normal
    /// cell spacing `dy`. The default (`new()`) sets this equal to
    /// `length_scale` so pre-existing callers retain the same dimensional
    /// scale when no grid information is available.
    pub wall_normal_spacing: T,
    /// Physical LES filter width Δ = (dx·dy·dz)^(1/3) [m].
    pub filter_width: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> MixingLengthModel<T> {
    /// Create a new mixing length model.
    ///
    /// Sets `wall_normal_spacing = length_scale` and `filter_width =
    /// length_scale` (appropriate when grid spacing and mixing length are
    /// comparable). For physically correct gradient computation use
    /// [`MixingLengthModel::with_filter_width`].
    pub fn new(length_scale: T) -> Self {
        let kappa = <T as FromPrimitive>::from_f64(
            cfd_core::physics::constants::physics::fluid::VON_KARMAN,
        )
        .expect("VON_KARMAN is an IEEE 754 representable f64 constant");
        Self {
            length_scale,
            kappa,
            wall_normal_spacing: length_scale,
            filter_width: length_scale,
        }
    }

    /// Create a mixing length model with the physically correct grid spacing.
    ///
    /// Computes Δ = (dx·dy·dz)^(1/3) per Deardorff (1970) as the resolved
    /// LES filter width, while the wall-normal derivative uses the physical
    /// `dy` spacing directly.
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
            wall_normal_spacing: dy,
            filter_width,
        }
    }

    /// Calculate the wall-normal shear rate using the stored physical spacing.
    ///
    /// The Prandtl mixing-length closure is wall-bounded, so the primary
    /// gradient is the streamwise velocity variation in the wall-normal
    /// direction.  Boundary values are evaluated with one-sided differences by
    /// the shared field helper.
    #[allow(clippy::similar_names)] // CFD derivatives use standard notation
    fn wall_normal_shear_rate_at(&self, velocity: &VelocityField<T>, idx: usize) -> T {
        let (nx, ny, _nz) = velocity.dimensions;
        let i = idx % nx;
        let j = (idx / nx) % ny;
        let k = idx / (nx * ny);

        let gradient = derivative_y(
            &velocity.components,
            nx,
            ny,
            i,
            j,
            k,
            self.wall_normal_spacing,
        );
        num_traits::Float::abs(gradient.x)
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
                let grad_u = self.wall_normal_shear_rate_at(&flow_field.velocity, idx);
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
                let grad_u = self.wall_normal_shear_rate_at(&flow_field.velocity, idx);
                let l_grad_u = self.length_scale * grad_u;
                l_grad_u * l_grad_u
            })
            .collect()
    }

    fn name(&self) -> &'static str {
        "MixingLength"
    }
}
