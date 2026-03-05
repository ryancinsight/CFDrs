//! Turbulence modeling abstractions and implementations
//!
//! Provides trait-based turbulence model implementations including
//! LES (Smagorinsky) and RANS (k-epsilon, Mixing Length) models.
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

pub mod constants;
pub mod des;
pub mod k_epsilon;
pub mod k_omega_sst;
pub mod mixing_length;
pub mod sigma;
pub mod spalart_allmaras;
pub mod vreman;
pub mod wall_functions;

pub use k_epsilon::{KEpsilonConstants, KEpsilonModel, KEpsilonState};
pub use mixing_length::MixingLengthModel;

use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::TurbulenceModel;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Smagorinsky subgrid-scale model for LES
#[derive(Debug, Clone)]
pub struct SmagorinskyModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Smagorinsky constant (typically 0.1-0.2)
    pub cs: T,
    /// Base Smagorinsky constant for dynamic model
    pub cs_base: T,
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

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> SmagorinskyModel<T> {
    /// Create a new Smagorinsky model with standard constant.
    ///
    /// Uses `filter_width = T::one()` (unit-cube default).  For physically
    /// correct LES viscosity, use [`SmagorinskyModel::with_filter_width`]
    /// and pass the actual cell dimensions.
    pub fn new(cs: T) -> Self {
        Self {
            cs,
            cs_base: cs,
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
        delta: T,
    ) -> T {
        let (nx, ny, nz) = flow_field.velocity.dimensions;

        // Calculate strain rate tensor components
        let mut s11 = T::zero();
        let mut s22 = T::zero();
        let mut s33 = T::zero();

        // Central differences for velocity gradients
        if i > 0 && i < nx - 1 {
            if let (Some(u_plus), Some(u_minus)) = (
                flow_field.velocity.get(i + 1, j, k),
                flow_field.velocity.get(i - 1, j, k),
            ) {
                let two = T::one() + T::one();
                s11 = (u_plus.x - u_minus.x) / (two * delta);
            }
        }

        if j > 0 && j < ny - 1 {
            if let (Some(v_plus), Some(v_minus)) = (
                flow_field.velocity.get(i, j + 1, k),
                flow_field.velocity.get(i, j - 1, k),
            ) {
                let two = T::one() + T::one();
                s22 = (v_plus.y - v_minus.y) / (two * delta);
            }
        }

        if k > 0 && k < nz - 1 {
            if let (Some(w_plus), Some(w_minus)) = (
                flow_field.velocity.get(i, j, k + 1),
                flow_field.velocity.get(i, j, k - 1),
            ) {
                let two = T::one() + T::one();
                s33 = (w_plus.z - w_minus.z) / (two * delta);
            }
        }

        // Off-diagonal components: Sij = 0.5 * (∂ui/∂xj + ∂uj/∂xi)
        let mut s12 = T::zero();
        let mut s13 = T::zero();
        let mut s23 = T::zero();

        // S12 = 0.5 * (∂u/∂y + ∂v/∂x)
        if i > 0 && i < nx - 1 && j > 0 && j < ny - 1 {
            if let (Some(u_jp), Some(u_jm), Some(v_ip), Some(v_im)) = (
                flow_field.velocity.get(i, j + 1, k),
                flow_field.velocity.get(i, j - 1, k),
                flow_field.velocity.get(i + 1, j, k),
                flow_field.velocity.get(i - 1, j, k),
            ) {
                let two = T::one() + T::one();
                let du_dy = (u_jp.x - u_jm.x) / (two * delta);
                let dv_dx = (v_ip.y - v_im.y) / (two * delta);
                s12 = (du_dy + dv_dx) / two;
            }
        }

        // S13 = 0.5 * (∂u/∂z + ∂w/∂x)
        if i > 0 && i < nx - 1 && k > 0 && k < nz - 1 {
            if let (Some(u_kp), Some(u_km), Some(w_ip), Some(w_im)) = (
                flow_field.velocity.get(i, j, k + 1),
                flow_field.velocity.get(i, j, k - 1),
                flow_field.velocity.get(i + 1, j, k),
                flow_field.velocity.get(i - 1, j, k),
            ) {
                let two = T::one() + T::one();
                let du_dz = (u_kp.x - u_km.x) / (two * delta);
                let dw_dx = (w_ip.z - w_im.z) / (two * delta);
                s13 = (du_dz + dw_dx) / two;
            }
        }

        // S23 = 0.5 * (∂v/∂z + ∂w/∂y)
        if j > 0 && j < ny - 1 && k > 0 && k < nz - 1 {
            if let (Some(v_kp), Some(v_km), Some(w_jp), Some(w_jm)) = (
                flow_field.velocity.get(i, j, k + 1),
                flow_field.velocity.get(i, j, k - 1),
                flow_field.velocity.get(i, j + 1, k),
                flow_field.velocity.get(i, j - 1, k),
            ) {
                let two = T::one() + T::one();
                let dv_dz = (v_kp.y - v_km.y) / (two * delta);
                let dw_dy = (w_jp.z - w_jm.z) / (two * delta);
                s23 = (dv_dz + dw_dy) / two;
            }
        }

        // Strain rate magnitude: |S| = sqrt(2 * Sij * Sij)
        let two = T::one() + T::one();
        let s_mag_sq =
            two * (s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23));
        num_traits::Float::sqrt(s_mag_sq)
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> TurbulenceModel<T>
    for SmagorinskyModel<T>
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
                    let strain_rate =
                        self.calculate_strain_rate_at_point(flow_field, i, j, k, delta);
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
        let cs_delta = self.cs * self.filter_width;

        let mut tke = Vec::with_capacity(nx * ny * nz);

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let strain_rate =
                        self.calculate_strain_rate_at_point(flow_field, i, j, k, delta);
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
