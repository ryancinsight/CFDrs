//! Turbulence modeling abstractions and implementations
//!
//! Provides trait-based turbulence model implementations including
//! LES (Smagorinsky, dynamic Smagorinsky, Sigma, Vreman, AMD) and
//! RANS (k-epsilon, Mixing Length) models.
//!
//! Moved from cfd-core to cfd-3d as these models are specific to 3D flow fields.
//!
//! # Theorem — Smagorinsky SGS Model (Smagorinsky 1963)
//!
//! The subgrid-scale eddy viscosity in the Smagorinsky LES model is
//!
//! ```text
//! ν_t = (C_s Δ)² |S|
//! ```
//!
//! where $C_s \in [0.1, 0.2]$ is the Smagorinsky constant, $\Delta$ is the
//! filter width (grid spacing), and $|S| = \sqrt{2 S_{ij} S_{ij}}$ is the
//! resolved strain-rate magnitude.
//!
//! **Proof sketch.** Dimensional analysis of the SGS stress tensor $\tau_{ij}
//! = -2 \nu_t \bar{S}_{ij}$ combined with Kolmogorov’s equilibrium hypothesis
//! (dissipation rate $\epsilon \sim \nu_t |S|^2 / \Delta^2$) gives $\nu_t
//! \propto \Delta^2 |S|$.
//!
//! # Theorem — Kolmogorov –5/3 Spectrum (Kolmogorov 1941)
//!
//! In the inertial subrange of homogeneous isotropic turbulence, the
//! energy spectrum follows
//!
//! ```text
//! E(k) = C_K ε^{2/3} k^{-5/3}
//! ```
//!
//! where $C_K \approx 1.5$ is the Kolmogorov constant. The Smagorinsky model
//! effectively assumes this spectral form to derive $C_s$ from $C_K$.
//!
//! **Reference:** Smagorinsky, J., "General circulation experiments with the
//! primitive equations", Mon. Wea. Rev. 91, 1963, pp. 99–164.

pub mod amd;
pub mod constants;
pub mod des;
pub mod dynamic_smagorinsky;
pub(crate) mod field_ops;
pub mod k_epsilon;
pub mod k_omega_sst;
pub mod mixing_length;
pub mod sigma;
pub mod spalart_allmaras;
pub mod vreman;
pub mod wall_functions;

pub use amd::AnisotropicMinimumDissipationModel;
pub use dynamic_smagorinsky::DynamicSmagorinskyModel;
pub use k_epsilon::{KEpsilonConstants, KEpsilonModel, KEpsilonState};
pub use mixing_length::MixingLengthModel;

use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::TurbulenceModel;
use nalgebra::RealField;
use num_traits::FromPrimitive;

use self::field_ops::{strain_magnitude, velocity_gradient_tensor};

/// Smagorinsky subgrid-scale model for LES
#[derive(Debug, Clone)]
pub struct SmagorinskyModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Smagorinsky constant (typically 0.1-0.2)
    pub cs: T,
    /// Base Smagorinsky constant for dynamic model
    pub cs_base: T,
    /// Physical grid spacing in the x direction [m].
    pub dx: T,
    /// Physical grid spacing in the y direction [m].
    pub dy: T,
    /// Physical grid spacing in the z direction [m].
    pub dz: T,
    /// Physical LES filter width Δ = (dx·dy·dz)^(1/3) [m].
    ///
    /// # Theorem — Geometric Mean Filter Width (Deardorff 1970)
    ///
    /// For anisotropic Cartesian cells, Deardorff (1970) showed that the
    /// geometric mean Δ = (δx·δy·δz)^(1/3) minimises aliasing of the
    /// resolved-to-subgrid energy transfer and correctly represents the
    /// local smallest resolved scale regardless of cell aspect ratio.
    ///
    /// **Proof sketch.** The subgrid-scale stress τ_ij = −2 νₜ S̄ᵢⱼ requires
    /// the LES filter scale to represent the physical grid cut-off.  For a
    /// unit cube (Δ = 1), the filter couples the SGS model to cell count
    /// rather than physical size, making νₜ wrong by orders of magnitude when
    /// actual grid spacing differs from 1.
    ///
    /// **Reference**: Deardorff, J.W. (1970). "A numerical study of three-
    /// dimensional turbulent channel flow at large Reynolds numbers."
    /// *J. Fluid Mech.* 41:453–480.
    pub filter_width: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + num_traits::Float>
    SmagorinskyModel<T>
{
    /// Create a new Smagorinsky model with standard constant.
    ///
    /// Uses `filter_width = T::one()` (unit-cube default).  For physically
    /// correct LES viscosity, use [`SmagorinskyModel::with_filter_width`]
    /// and pass the actual cell dimensions.
    pub fn new(cs: T) -> Self {
        Self {
            cs,
            cs_base: cs,
            dx: T::one(),
            dy: T::one(),
            dz: T::one(),
            filter_width: T::one(),
        }
    }

    /// Create a Smagorinsky model with the physically correct filter width.
    ///
    /// Computes Δ = (dx·dy·dz)^(1/3) per Deardorff (1970).
    ///
    /// # Arguments
    /// * `cs` — Smagorinsky constant (typically 0.10–0.17 for channel flow)
    /// * `dx`, `dy`, `dz` — physical cell dimensions [m]
    ///
    /// # Example
    /// ```rust,ignore
    /// // 50 µm × 50 µm × 100 µm cell: Δ ≈ 62.9 µm
    /// let smag = SmagorinskyModel::with_filter_width(0.1, 50e-6, 50e-6, 100e-6);
    /// ```
    pub fn with_filter_width(cs: T, dx: T, dy: T, dz: T) -> Self {
        let one_third = T::one() / (T::one() + T::one() + T::one());
        let filter_width = num_traits::Float::powf(dx * dy * dz, one_third);
        Self {
            cs,
            cs_base: cs,
            dx,
            dy,
            dz,
            filter_width,
        }
    }

    /// Calculate strain rate magnitude at a grid point
    fn calculate_strain_rate_at_point(
        &self,
        flow_field: &FlowField<T>,
        i: usize,
        j: usize,
        k: usize,
    ) -> T {
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let gradient = velocity_gradient_tensor(
            &flow_field.velocity.components,
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

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + num_traits::Float>
    TurbulenceModel<T> for SmagorinskyModel<T>
{
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        // Physically correct filter width Δ = (dx·dy·dz)^(1/3); stored at construction time.
        // Bug fix: the previous `Δ = 1/nx` used grid count instead of physical cell size,
        // making νₜ wrong by orders of magnitude for any domain that isn't a unit cube
        // subdivided into exactly nx uniform cells (Deardorff 1970, Pope 2000 §13.2).
        let delta = self.filter_width;
        let prefactor = self.cs * self.cs * delta * delta;

        let mut viscosity = Vec::with_capacity(nx * ny * nz);

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let strain_rate = self.calculate_strain_rate_at_point(flow_field, i, j, k);
                    viscosity.push(prefactor * strain_rate);
                }
            }
        }
        viscosity
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        // Estimate SGS TKE from strain rate: k ≈ (Cs · Δ · |S|)²
        // Uses the physically correct filter width (see turbulent_viscosity bug note).
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let delta = self.filter_width;
        let cs_delta = self.cs * delta;

        let mut tke = Vec::with_capacity(nx * ny * nz);

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let strain_rate = self.calculate_strain_rate_at_point(flow_field, i, j, k);
                    let cs_delta_s = cs_delta * strain_rate;
                    tke.push(cs_delta_s * cs_delta_s);
                }
            }
        }
        tke
    }

    fn name(&self) -> &'static str {
        "Smagorinsky"
    }
}
