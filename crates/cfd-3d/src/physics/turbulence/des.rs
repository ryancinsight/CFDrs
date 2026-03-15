//! Detached-Eddy Simulation (DES) hybrid RANS-LES model.
//!
//! # Theorem — DES Length Scale (Spalart et al. 1997)
//!
//! DES modifies the RANS destruction length scale by replacing the wall
//! distance d with:
//!
//! ```text
//! l_DES = min(d, C_DES · Δ_max)
//! ```
//!
//! where Δ_max = max(dx, dy, dz) is the largest cell dimension and
//! C_DES = 0.65 (calibrated from homogeneous turbulence).
//!
//! **Physical interpretation.** Near walls (d < C_DES · Δ), l_DES = d and the
//! model reduces to RANS (Spalart-Allmaras).  Away from walls
//! (d > C_DES · Δ), l_DES = C_DES · Δ and the model operates as LES.
//! The transition is governed by Δ_max / d, providing automatic mode switching.
//!
//! **Proof sketch.** For d ≫ C_DES Δ, the SA production–destruction balance
//! reduces to ε ~ ν̃ · |S| / (C_DES Δ)², giving νₜ ~ (C_DES Δ)² |S|,
//! identical to the Smagorinsky SGS model.  (Shur et al. 1999, §2.)
//!
//! ## Algorithm
//!
//! ```text
//! For each grid point:
//!   1. Evaluate wall distance d_i (provided at construction via wall_distances).
//!   2. Compute Δ_max_i = max(dx, dy, dz) (provided at construction).
//!   3. l_DES_i = min(d_i, C_DES · Δ_max_i)
//!   4. Compute RANS νₜ using SA model with modified destruction length l_DES.
//!      (Approximated here by Smagorinsky with Δ = l_DES for tractability.)
//! ```
//!
//! ## References
//!
//! - Spalart, P.R., Jou, W.H., Strelets, M. & Allmaras, S.R. (1997).
//!   "Comments on the feasibility of LES for wings, and on a hybrid RANS/LES
//!   approach." *Advances in DNS/LES* 1:4–8.
//! - Shur, M., Spalart, P.R., Strelets, M. & Travin, A. (1999). "Detached-
//!   eddy simulation of an airfoil at high angle of attack." *IUTAM Symp.*

use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;
use nalgebra::RealField;
use num_traits::FromPrimitive;

use super::constants::{DES_C_DES, SMAGORINSKY_CS_DEFAULT};

/// DES hybrid RANS-LES model (Spalart et al. 1997).
///
/// Blends RANS behaviour near walls with LES behaviour away from walls
/// by selecting the minimum of the RANS wall-distance and LES filter scale.
#[derive(Debug, Clone)]
pub struct DESModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> {
    /// DES constant C_DES = 0.65.
    pub c_des: T,
    /// Wall distance for each grid point (order: x-innermost, then y, then z).
    ///
    /// Must have length == nx·ny·nz.  Set via [`DESModel::with_wall_distances`].
    /// A uniform large value (e.g. 1.0 m) produces pure LES everywhere.
    pub wall_distances: Vec<T>,
    /// Maximum cell dimension Δ_max = max(dx, dy, dz) for each grid point.
    ///
    /// For uniform Cartesian grids, this is a constant.  For non-uniform grids
    /// supply the per-cell value.
    pub delta_max: Vec<T>,
    /// Smagorinsky constant C_s used for the LES sub-model away from walls.
    pub cs: T,
    /// Whether to use Delayed DES (DDES) shielding (Spalart et al. 2006).
    ///
    /// When `true`, the length scale is computed using the DDES shielding
    /// function f_d, which prevents premature LES activation inside attached
    /// boundary layers (modeled stress depletion).
    ///
    /// When `false` (default), the standard DES length scale
    /// `l = min(d, C_DES * delta)` is used.
    ///
    /// **Reference**: Spalart et al. (2006), *Theor. Comput. Fluid Dyn.* 20:181-195.
    pub use_ddes: bool,
    /// Molecular kinematic viscosity [m²/s] for DDES shielding computation.
    /// Only used when `use_ddes` is true. Default: 1.5e-5 (air at 20°C).
    pub molecular_viscosity: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> DESModel<T> {
    /// Create a DES model with uniform wall distance and filter width.
    ///
    /// # Arguments
    /// * `n_points` — number of grid points (nx·ny·nz)
    /// * `uniform_wall_distance` — wall distance [m] applied uniformly (use a
    ///   large value for pure LES, 0 for pure RANS)
    /// * `dx`, `dy`, `dz` — uniform cell dimensions [m]
    pub fn new(n_points: usize, uniform_wall_distance: T, dx: T, dy: T, dz: T) -> Self {
        let c_des = <T as FromPrimitive>::from_f64(DES_C_DES)
            .expect("DES_C_DES is an IEEE 754 representable f64 constant");
        let delta_max = num_traits::Float::max(dx, num_traits::Float::max(dy, dz));
        Self {
            c_des,
            wall_distances: vec![uniform_wall_distance; n_points],
            delta_max: vec![delta_max; n_points],
            cs: <T as FromPrimitive>::from_f64(SMAGORINSKY_CS_DEFAULT)
                .expect("SMAGORINSKY_CS_DEFAULT is an IEEE 754 representable f64 constant"),
            use_ddes: false,
            molecular_viscosity: <T as FromPrimitive>::from_f64(1.5e-5)
                .expect("1.5e-5 is an IEEE 754 representable f64 constant"),
        }
    }

    /// Create a DES model with per-point wall distances and delta_max values.
    pub fn with_wall_distances(wall_distances: Vec<T>, delta_max: Vec<T>, cs: T) -> Self {
        let c_des = <T as FromPrimitive>::from_f64(DES_C_DES)
            .expect("DES_C_DES is an IEEE 754 representable f64 constant");
        Self {
            c_des,
            wall_distances,
            delta_max,
            cs,
            use_ddes: false,
            molecular_viscosity: <T as FromPrimitive>::from_f64(1.5e-5)
                .expect("1.5e-5 is an IEEE 754 representable f64 constant"),
        }
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> TurbulenceModel<T>
    for DESModel<T>
{
    /// Compute DES eddy viscosity.
    ///
    /// At each grid point, uses l_DES = min(d, C_DES·Δ_max) as the effective
    /// length scale for the Smagorinsky LES sub-model.
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let n = nx * ny * nz;
        let mut viscosity = Vec::with_capacity(n);
        let two = <T as FromPrimitive>::from_f64(2.0)
            .expect("2.0 is representable in all IEEE 754 types");

        for idx in 0..n {
            let d = if idx < self.wall_distances.len() {
                self.wall_distances[idx]
            } else {
                T::one()
            };
            let delta = if idx < self.delta_max.len() {
                self.delta_max[idx]
            } else {
                T::one()
            };

            let i = idx % nx;
            let j = (idx / nx) % ny;
            let k = idx / (nx * ny);

            // Compute strain rate and vorticity magnitudes for both
            // the eddy viscosity formula and DDES shielding (if enabled).
            // We use the grid spacing (delta) for finite differences here,
            // rather than l_des, to get physical velocity gradients.
            let mut s_sq = T::zero();
            let mut w_sq = T::zero();
            if i > 0 && i < nx - 1 {
                if let (Some(vp), Some(vm)) = (
                    flow_field.velocity.get(i + 1, j, k),
                    flow_field.velocity.get(i - 1, j, k),
                ) {
                    let s11 = (vp.x - vm.x) / (two * delta);
                    let s12 = (vp.y - vm.y) / (two * delta);
                    let s13 = (vp.z - vm.z) / (two * delta);
                    s_sq += s11 * s11;
                    // Cross terms for vorticity
                    w_sq += s12 * s12 + s13 * s13;
                }
            }
            if j > 0 && j < ny - 1 {
                if let (Some(vp), Some(vm)) = (
                    flow_field.velocity.get(i, j + 1, k),
                    flow_field.velocity.get(i, j - 1, k),
                ) {
                    let s22 = (vp.y - vm.y) / (two * delta);
                    let s21 = (vp.x - vm.x) / (two * delta);
                    let s23 = (vp.z - vm.z) / (two * delta);
                    s_sq += s22 * s22;
                    w_sq += s21 * s21 + s23 * s23;
                }
            }
            if k > 0 && k < nz - 1 {
                if let (Some(vp), Some(vm)) = (
                    flow_field.velocity.get(i, j, k + 1),
                    flow_field.velocity.get(i, j, k - 1),
                ) {
                    let s33 = (vp.z - vm.z) / (two * delta);
                    let s31 = (vp.x - vm.x) / (two * delta);
                    let s32 = (vp.y - vm.y) / (two * delta);
                    s_sq += s33 * s33;
                    w_sq += s31 * s31 + s32 * s32;
                }
            }

            // Compute the effective length scale
            let l_des = if self.use_ddes {
                // DDES shielding: use f_d to delay LES activation in boundary layers.
                // Previous iteration's nu_t is approximated using standard DES l for this point.
                let l_standard = num_traits::Float::min(d, self.c_des * delta);
                let s_mag_approx = num_traits::Float::sqrt(two * s_sq);
                let nu_t_approx = self.cs * self.cs * l_standard * l_standard * s_mag_approx;

                // Convert to f64 for the shielding function
                let nu_t_f64 = nu_t_approx.to_f64().unwrap_or(0.0);
                let nu_f64 = self.molecular_viscosity.to_f64().unwrap_or(1.5e-5);
                let d_f64 = d.to_f64().unwrap_or(1.0);
                let strain_f64 = num_traits::Float::sqrt(s_sq).to_f64().unwrap_or(0.0);
                let vort_f64 = num_traits::Float::sqrt(w_sq).to_f64().unwrap_or(0.0);

                let fd = ddes_shielding(nu_t_f64, nu_f64, d_f64, strain_f64, vort_f64);
                let c_des_f64 = self.c_des.to_f64().unwrap_or(DES_C_DES);
                let delta_f64 = delta.to_f64().unwrap_or(1.0);
                let l_ddes_f64 = ddes_length_scale(d_f64, c_des_f64, delta_f64, fd);
                <T as FromPrimitive>::from_f64(l_ddes_f64)
                    .expect("DDES length scale is a finite f64 value")
            } else {
                // Standard DES: l = min(d, C_DES · Δ_max)
                num_traits::Float::min(d, self.c_des * delta)
            };

            // Compute strain rate magnitude using l_des as spacing (original formulation)
            let mut s_sq_local = T::zero();
            if i > 0 && i < nx - 1 {
                if let (Some(vp), Some(vm)) = (
                    flow_field.velocity.get(i + 1, j, k),
                    flow_field.velocity.get(i - 1, j, k),
                ) {
                    let s11 = (vp.x - vm.x) / (two * l_des);
                    s_sq_local += s11 * s11;
                }
            }
            if j > 0 && j < ny - 1 {
                if let (Some(vp), Some(vm)) = (
                    flow_field.velocity.get(i, j + 1, k),
                    flow_field.velocity.get(i, j - 1, k),
                ) {
                    let s22 = (vp.y - vm.y) / (two * l_des);
                    s_sq_local += s22 * s22;
                }
            }
            if k > 0 && k < nz - 1 {
                if let (Some(vp), Some(vm)) = (
                    flow_field.velocity.get(i, j, k + 1),
                    flow_field.velocity.get(i, j, k - 1),
                ) {
                    let s33 = (vp.z - vm.z) / (two * l_des);
                    s_sq_local += s33 * s33;
                }
            }
            let s_mag = num_traits::Float::sqrt(two * s_sq_local);
            viscosity.push(self.cs * self.cs * l_des * l_des * s_mag);
        }
        viscosity
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        self.turbulent_viscosity(flow_field)
    }

    fn name(&self) -> &'static str {
        "DES"
    }
}

// ── DDES (Delayed Detached-Eddy Simulation) ─────────────────────────────────

/// Spalart et al. (2006) Delayed Detached-Eddy Simulation (DDES) shielding function.
///
/// ## Theorem — DDES Shielding (Spalart et al. 2006)
///
/// The standard DES model switches from RANS to LES based on a simple
/// length-scale comparison: l = min(d, C_DES·Δ). This can cause premature
/// LES activation inside attached boundary layers ("modeled stress depletion").
///
/// The DDES modification introduces a shielding function f_d that delays
/// the LES activation:
///
/// ```text
/// l_DDES = d − f_d · max(0, d − C_DES · Δ)
/// ```
///
/// where the shielding function is:
///
/// ```text
/// f_d = 1 − tanh((C_d1 · r_d)^(C_d2))
/// ```
///
/// and the RANS indicator r_d is:
///
/// ```text
/// r_d = (ν_t + ν) / (κ² · d² · √(0.5·(S² + Ω²)))
/// ```
///
/// where κ = 0.41, d is wall distance, S is strain rate, Ω is vorticity.
///
/// **Physical basis**: In the boundary layer, r_d ≈ 1 (RANS regime) →
/// f_d ≈ 0 → l_DDES ≈ d (pure RANS). Away from walls, r_d → 0 →
/// f_d → 1 → l_DDES ≈ C_DES·Δ (pure LES). The smooth transition
/// prevents grid-induced separation.
///
/// Constants: C_d1 = 8.0, C_d2 = 3.0 (Spalart et al. 2006).
///
/// **Reference**: Spalart, P.R., Deck, S., Shur, M.L., Squires, K.D.,
/// Strelets, M.Kh. & Travin, A. (2006). "A New Version of Detached-Eddy
/// Simulation, Resistant to Ambiguous Grid Densities",
/// *Theor. Comput. Fluid Dyn.* 20:181-195.

/// DDES constant C_d1 = 8.0 (Spalart et al. 2006).
pub const CD1: f64 = 8.0;
/// DDES constant C_d2 = 3.0 (Spalart et al. 2006).
pub const CD2: f64 = 3.0;

/// Compute the DDES shielding function f_d.
///
/// Returns f_d ∈ [0, 1]:
/// - f_d ≈ 0 in RANS regions (boundary layer) → preserves RANS
/// - f_d ≈ 1 in LES regions (separated flow) → activates LES
pub fn ddes_shielding(
    nu_t: f64,          // turbulent kinematic viscosity [m²/s]
    nu: f64,            // molecular kinematic viscosity [m²/s]
    wall_distance: f64, // distance to nearest wall [m]
    strain_rate: f64,   // strain rate magnitude |S| [s⁻¹]
    vorticity: f64,     // vorticity magnitude |Ω| [s⁻¹]
) -> f64 {
    let kappa = 0.41;
    let denominator = kappa * kappa * wall_distance * wall_distance
        * (0.5 * (strain_rate * strain_rate + vorticity * vorticity))
            .sqrt()
            .max(1e-30);
    let r_d = (nu_t + nu) / denominator;
    1.0 - (CD1 * r_d).powf(CD2).tanh()
}

/// Compute DDES length scale.
///
/// Returns l_DDES = d − f_d · max(0, d − C_DES · Δ).
/// When f_d = 0 (RANS region), l = d.
/// When f_d = 1 and d > C_DES·Δ, l = C_DES·Δ (LES region).
pub fn ddes_length_scale(wall_distance: f64, c_des: f64, delta: f64, f_d: f64) -> f64 {
    wall_distance - f_d * (wall_distance - c_des * delta).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ddes_shielding_rans_region() {
        // Near wall with moderate turbulence: high r_d → f_d ≈ 0
        let nu_t = 1e-5; // small eddy viscosity near wall
        let nu = 1.5e-5; // molecular viscosity (air)
        let d = 0.001; // 1 mm from wall
        let strain = 50.0;
        let vorticity = 50.0;
        let fd = ddes_shielding(nu_t, nu, d, strain, vorticity);
        assert!(
            fd < 0.1,
            "f_d should be near 0 in RANS region (near wall), got {fd}"
        );
    }

    #[test]
    fn test_ddes_shielding_les_region() {
        // Far from wall with strong shear: low r_d → f_d ≈ 1
        let nu_t = 1e-4;
        let nu = 1.5e-5;
        let d = 1.0; // 1 m from wall
        let strain = 1000.0;
        let vorticity = 1000.0;
        let fd = ddes_shielding(nu_t, nu, d, strain, vorticity);
        assert!(
            fd > 0.9,
            "f_d should be near 1 in LES region (far from wall), got {fd}"
        );
    }

    #[test]
    fn test_ddes_shielding_bounded() {
        // f_d must always be in [0, 1] across a wide parameter range
        let params = [
            (1e-6, 1e-6, 0.0001, 1.0, 1.0),
            (1e-3, 1.5e-5, 0.01, 100.0, 100.0),
            (1e-1, 1.5e-5, 1.0, 1000.0, 1000.0),
            (1e-4, 1e-6, 0.1, 10.0, 10.0),
            (1e-2, 1e-5, 0.5, 500.0, 200.0),
        ];
        for (nu_t, nu, d, s, w) in params {
            let fd = ddes_shielding(nu_t, nu, d, s, w);
            assert!(
                (0.0..=1.0).contains(&fd),
                "f_d must be in [0,1], got {fd} for nu_t={nu_t}, nu={nu}, d={d}, S={s}, W={w}"
            );
        }
    }

    #[test]
    fn test_ddes_length_scale_rans_limit() {
        // f_d = 0 → l_DDES = d (pure RANS)
        let d = 0.005;
        let c_des = 0.65;
        let delta = 0.01;
        let l = ddes_length_scale(d, c_des, delta, 0.0);
        assert!(
            (l - d).abs() < 1e-15,
            "With f_d=0, l_DDES should equal d={d}, got {l}"
        );
    }

    #[test]
    fn test_ddes_length_scale_les_limit() {
        // f_d = 1 and d > C_DES·Δ → l_DDES = C_DES·Δ (pure LES)
        let d = 1.0;
        let c_des = 0.65;
        let delta = 0.01;
        let expected = c_des * delta;
        let l = ddes_length_scale(d, c_des, delta, 1.0);
        assert!(
            (l - expected).abs() < 1e-15,
            "With f_d=1 and d > C_DES*Δ, l_DDES should equal C_DES*Δ={expected}, got {l}"
        );
    }

    #[test]
    fn test_des_default_not_ddes() {
        let model = DESModel::<f64>::new(8, 0.5, 0.01, 0.01, 0.01);
        assert!(
            !model.use_ddes,
            "DES model must default to standard DES (use_ddes = false)"
        );
    }

    #[test]
    fn test_des_with_wall_distances_default_not_ddes() {
        let model = DESModel::<f64>::with_wall_distances(
            vec![0.5; 8],
            vec![0.01; 8],
            0.1,
        );
        assert!(
            !model.use_ddes,
            "with_wall_distances must default to standard DES (use_ddes = false)"
        );
    }

    #[test]
    fn test_des_ddes_flag_toggleable() {
        use cfd_core::physics::fluid_dynamics::fields::FlowField;

        let n = 27; // 3x3x3
        let mut model = DESModel::<f64>::new(n, 0.5, 0.1, 0.1, 0.1);

        // Create a flow field with some velocity gradients
        let mut flow = FlowField::<f64>::new(3, 3, 3);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    if let Some(v) = flow.velocity.get_mut(i, j, k) {
                        *v = nalgebra::Vector3::new(
                            (i as f64) * 0.5,
                            (j as f64) * 0.3,
                            (k as f64) * 0.1,
                        );
                    }
                }
            }
        }

        // Compute with standard DES
        model.use_ddes = false;
        let visc_des: Vec<f64> = TurbulenceModel::turbulent_viscosity(&model, &flow);

        // Compute with DDES
        model.use_ddes = true;
        let visc_ddes: Vec<f64> = TurbulenceModel::turbulent_viscosity(&model, &flow);

        // The results should differ because DDES uses a different length scale
        // (At least at some interior points where gradients are non-zero)
        let any_differ = visc_des
            .iter()
            .zip(visc_ddes.iter())
            .any(|(a, b)| (a - b).abs() > 1e-20);

        // For uniform wall distance of 0.5 and delta of 0.1, C_DES*delta = 0.065.
        // Standard DES: l = min(0.5, 0.065) = 0.065
        // DDES with f_d that depends on local flow: l_DDES varies.
        // They should produce different results at interior points.
        assert!(
            any_differ,
            "DDES should produce different eddy viscosity than standard DES for non-trivial flow"
        );
    }
}
