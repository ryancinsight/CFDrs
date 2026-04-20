//! Interface reconstruction methods for VOF
//!
//! # Theorem — PLIC Volume Conservation (Youngs 1982, Scardovelli & Zaleski 2000)
//!
//! Given a volume fraction field $\alpha(x)$ on a structured grid, the
//! Piecewise Linear Interface Construction (PLIC) method reconstructs a planar
//! interface in each cell $\Omega_i$ as a hyperplane
//!
//! ```text
//! n · x = C
//! ```
//!
//! where $\mathbf{n}$ is the outward interface normal estimated from a
//! dominant-axis height function when the interface is graph-like, or from
//! `-∇α / |∇α|` as a Youngs fallback, and $C$ is chosen so that the clipped
//! volume matches $\alpha_i |\Omega_i|$ exactly (to machine precision). This
//! ensures global volume conservation.
//!
//! **Proof sketch.** For a convex cell with normal $\mathbf{n}$ fixed, the
//! clipped volume $V(C)$ is a monotonically increasing, piecewise-polynomial
//! function of $C$. Inverting $V(C) = \alpha \cdot |\Omega|$ via the closed-form
//! Pilliod & Puckett (2004) formula gives the unique $C$.
//!
//! ## Mathematical Foundation
//! ```math
//! Γᵢ = { x ∈ Ω_i : n̂ᵢ · x = Cᵢ }
//! ```
//!
//! where `n̂ᵢ` is the unit outward normal and `Cᵢ` is chosen so that
//!
//! ```math
//! ∫_{Ω_i ∩ {n̂·x ≤ C}} dV = αᵢ |Ω_i|
//! ```
//!
//! **Normal Estimation** (Youngs 1984; Cummins et al. 2005): The normal is
//! estimated from a dominant-axis height function when a local graph exists,
//! otherwise from the discrete gradient of α:
//!
//! ```math
//! n̂ᵢ ≈ -∇α / |∇α|  (mixed finite-difference gradient)
//! ```
//!
//! **Plane Constant** `Cᵢ`: Obtained by binary bisection of the monotone function
//! `V(C) = ∫_{n·x≤C} dV`, using the analytical formula of Scardovelli & Zaleski
//! (2000) (Eqs. 2.34–2.38).
//!
//! **Accuracy**:
//! - Interface position: O(h²) for smooth interfaces
//! - Normal direction: O(h) for Youngs' gradient stencil
//!
//! ### Curvature Computation Theorem (continuum surface force)
//!
//! **Statement**: The interface curvature κ = −∇·n̂ computed from the
//! interface normal field **n̂** satisfies:
//!
//! ```math
//! κ = −(∂n̂_x/∂x + ∂n̂_y/∂y + ∂n̂_z/∂z)
//! ```
//!
//! **Reference**: Brackbill, J.U., Kothe, D.B. & Zemach, C. (1992).
//! "A continuum method for modeling surface tension".
//! J. Comput. Phys. 100:335–354.

use super::config::{constants, VofConfig, VOF_EPSILON, VOF_INTERFACE_LOWER, VOF_INTERFACE_UPPER};
use super::plic_geometry::volume_under_plane_3d;
use super::solver::VofSolver;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;
use std::cmp::Ordering;

// ── Height-Function Normal Estimation ────────────────────────────────────────

#[derive(Debug)]
struct DirectionalHeightCache<T> {
    x: Vec<T>,
    y: Vec<T>,
    z: Vec<T>,
}

#[derive(Clone, Copy, Debug)]
struct HybridSelection<T> {
    reference: Vector3<T>,
    axis: Option<usize>,
}

/// Dominant-axis threshold for the hybrid VOF interface reconstruction.
///
/// The value is kept above 1/sqrt(3) so a genuinely non-graph-like 3D cell
/// can still fall back to the Youngs gradient path.
const HYBRID_AXIS_DOMINANCE: f64 = 0.70;

/// Height-function normal estimation for PLIC interface reconstruction
/// (Cummins, Francois & Kothe 2005).
///
/// ## Theorem — Height-Function Normal Estimation
///
/// The height function H_i in column i is defined as the sum of volume
/// fractions in the column:
///
/// ```text
/// H_i = Σ_j α_{i,j} · Δy_j
/// ```
///
/// The interface normal components are then:
///
/// ```text
/// n_x = −(H_{i+1} − H_{i−1}) / (2Δx)
/// n_y = 1
/// ```
///
/// normalized: n = (n_x, n_y) / |n|
///
/// This achieves O(Δx²) accuracy in the normal direction, compared to
/// O(Δx) for the Youngs (1982) gradient method. The improvement comes
/// from summing volume fractions over an entire column (integrating out
/// one spatial direction), which filters out cell-level noise.
///
/// For 3D, the height function extends to a 3×3 column stencil with
/// derivatives in both x and z directions:
///
/// ```text
/// n_x = −∂H/∂x,  n_z = −∂H/∂z,  n_y = 1
/// ```
///
/// **Reference**: Cummins, S.J., Francois, M.M. & Kothe, D.B. (2005).
/// "Estimating curvature from volume fractions",
/// *Computers & Structures* 83(6-7):425-434.
#[must_use]
pub fn height_function_normal_2d(
    alpha: &[Vec<f64>], // volume fraction field [nx][ny]
    i: usize,
    _j: usize,
    dx: f64,
    dy: f64,
) -> [f64; 2] {
    let nx = alpha.len();
    let ny = if nx > 0 { alpha[0].len() } else { 0 };

    // Height function: sum volume fractions in column i
    let column_height = |col: usize| -> f64 { (0..ny).map(|j| alpha[col][j] * dy).sum() };

    // Central difference for normal using clamped boundary indices
    let i_left = if i > 0 { i - 1 } else { 0 };
    let i_right = if i + 1 < nx { i + 1 } else { nx - 1 };

    let h_left = column_height(i_left);
    let h_right = column_height(i_right);

    let nx_comp = -(h_right - h_left) / (2.0 * dx);
    let ny_comp = 1.0;
    let mag = (nx_comp * nx_comp + ny_comp * ny_comp).sqrt();

    [nx_comp / mag, ny_comp / mag]
}

/// Height-function normal estimation for 3D PLIC interface reconstruction.
///
/// Extends the 2D height-function method to 3D by computing column sums
/// along the y-axis and taking central differences in both x and z
/// directions. The normal is:
///
/// ```text
/// n_x = −∂H/∂x,  n_z = −∂H/∂z,  n_y = 1
/// ```
///
/// normalized to unit length. This achieves O(Δx²) accuracy in 3D.
#[must_use]
pub fn height_function_normal_3d(
    alpha: &[f64], // volume fraction field (flat array, indexed [i*ny*nz + j*nz + k])
    nx: usize,
    ny: usize,
    nz: usize,
    i: usize,
    _j: usize,
    k: usize,
    dx: f64,
    dy: f64,
    dz: f64,
) -> [f64; 3] {
    // Height function: sum volume fractions along y-axis for column (ci, ck)
    let column_height = |ci: usize, ck: usize| -> f64 {
        (0..ny)
            .map(|cj| alpha[ci * ny * nz + cj * nz + ck] * dy)
            .sum()
    };

    let i_left = if i > 0 { i - 1 } else { 0 };
    let i_right = if i + 1 < nx { i + 1 } else { nx - 1 };
    let k_left = if k > 0 { k - 1 } else { 0 };
    let k_right = if k + 1 < nz { k + 1 } else { nz - 1 };

    let dh_dx = -(column_height(i_right, k) - column_height(i_left, k)) / (2.0 * dx);
    let dh_dz = -(column_height(i, k_right) - column_height(i, k_left)) / (2.0 * dz);
    let ny_comp = 1.0;

    let mag = (dh_dx * dh_dx + ny_comp * ny_comp + dh_dz * dh_dz).sqrt();

    [dh_dx / mag, ny_comp / mag, dh_dz / mag]
}

/// Compute the Youngs (1982) gradient-based normal for a 2D volume fraction
/// field, for comparison with the height-function method.
///
/// Uses central differences: n = −∇α / |∇α|.
#[must_use]
pub fn youngs_normal_2d(alpha: &[Vec<f64>], i: usize, j: usize, dx: f64, dy: f64) -> [f64; 2] {
    let nx = alpha.len();
    let ny = if nx > 0 { alpha[0].len() } else { 0 };

    let i_left = if i > 0 { i - 1 } else { 0 };
    let i_right = if i + 1 < nx { i + 1 } else { nx - 1 };
    let j_left = if j > 0 { j - 1 } else { 0 };
    let j_right = if j + 1 < ny { j + 1 } else { ny - 1 };

    let da_dx = (alpha[i_right][j] - alpha[i_left][j]) / (2.0 * dx);
    let da_dy = (alpha[i][j_right] - alpha[i][j_left]) / (2.0 * dy);

    let mag = (da_dx * da_dx + da_dy * da_dy).sqrt();
    if mag < 1e-15 {
        return [0.0, 1.0];
    }

    // Normal points from fluid 1 to fluid 0: n = -∇α / |∇α|
    [-da_dx / mag, -da_dy / mag]
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>
    DirectionalHeightCache<T>
{
    fn build(solver: &VofSolver<T>) -> Self {
        let mut x_heights = vec![T::zero(); solver.ny * solver.nz];
        let mut y_heights = vec![T::zero(); solver.nx * solver.nz];
        let mut z_heights = vec![T::zero(); solver.nx * solver.ny];

        for k in 0..solver.nz {
            for j in 0..solver.ny {
                let mut height = T::zero();
                for i in 0..solver.nx {
                    height += solver.alpha[solver.index(i, j, k)] * solver.dx;
                }
                x_heights[k * solver.ny + j] = height;
            }
        }

        for k in 0..solver.nz {
            for i in 0..solver.nx {
                let mut height = T::zero();
                for j in 0..solver.ny {
                    height += solver.alpha[solver.index(i, j, k)] * solver.dy;
                }
                y_heights[k * solver.nx + i] = height;
            }
        }

        for j in 0..solver.ny {
            for i in 0..solver.nx {
                let mut height = T::zero();
                for k in 0..solver.nz {
                    height += solver.alpha[solver.index(i, j, k)] * solver.dz;
                }
                z_heights[j * solver.nx + i] = height;
            }
        }

        Self {
            x: x_heights,
            y: y_heights,
            z: z_heights,
        }
    }
}

/// Cache blocking parameters for optimized 3D traversal
const CACHE_BLOCK_SIZE_I: usize = 8;
const CACHE_BLOCK_SIZE_J: usize = 8;
const CACHE_BLOCK_SIZE_K: usize = 8;

use serde::{Deserialize, Serialize};

/// Interface reconstruction strategy
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum InterfaceReconstruction {
    /// Piecewise Linear Interface Calculation (PLIC)
    ///
    /// Uses Youngs' gradient-based normal estimation and the Scardovelli-Zaleski
    /// analytical volume formula for the plane constant.
    PLIC,
    /// Simple gradient-based reconstruction (fastest, least accurate)
    Gradient,
}

impl InterfaceReconstruction {
    /// Create reconstruction method based on configuration
    #[must_use]
    pub fn create(config: &VofConfig) -> Self {
        config.reconstruction_method
    }

    /// Reconstruct interface normals and curvature
    pub fn reconstruct<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        self,
        solver: &mut VofSolver<T>,
    ) {
        let height_cache = if matches!(self, Self::PLIC) {
            Some(DirectionalHeightCache::build(solver))
        } else {
            None
        };

        self.calculate_normals(solver, height_cache.as_ref());
        self.calculate_curvature(solver, height_cache.as_ref());
    }

    /// Calculate interface normal vectors using gradient of volume fraction.
    ///
    /// Uses cache blocking for improved memory access patterns on 3D grids.
    fn calculate_normals<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        self,
        solver: &mut VofSolver<T>,
        height_cache: Option<&DirectionalHeightCache<T>>,
    ) {
        let interface_lower = <T as FromPrimitive>::from_f64(VOF_INTERFACE_LOWER)
            .expect("VOF_INTERFACE_LOWER is an IEEE 754 representable f64 constant");
        let interface_upper = <T as FromPrimitive>::from_f64(VOF_INTERFACE_UPPER)
            .expect("VOF_INTERFACE_UPPER is an IEEE 754 representable f64 constant");

        for k_block in (1..solver.nz - 1).step_by(CACHE_BLOCK_SIZE_K) {
            for j_block in (1..solver.ny - 1).step_by(CACHE_BLOCK_SIZE_J) {
                for i_block in (1..solver.nx - 1).step_by(CACHE_BLOCK_SIZE_I) {
                    let k_end = (k_block + CACHE_BLOCK_SIZE_K).min(solver.nz - 1);
                    let j_end = (j_block + CACHE_BLOCK_SIZE_J).min(solver.ny - 1);
                    let i_end = (i_block + CACHE_BLOCK_SIZE_I).min(solver.nx - 1);

                    for k in k_block..k_end {
                        for j in j_block..j_end {
                            for i in i_block..i_end {
                                let idx = solver.index(i, j, k);
                                let alpha = solver.alpha[idx];

                                if alpha > interface_lower && alpha < interface_upper {
                                    match self {
                                        Self::PLIC => {
                                            let cache = height_cache.expect(
                                                "directional height cache must exist for PLIC",
                                            );
                                            let (normal, _) =
                                                self.plic_reconstruction(solver, cache, i, j, k);
                                            solver.normals[idx] = normal;
                                        }
                                        Self::Gradient => {
                                            let normal = self.calculate_gradient(solver, i, j, k);
                                            let epsilon = <T as FromPrimitive>::from_f64(
                                                VOF_EPSILON,
                                            )
                                            .expect("VOF_EPSILON is an IEEE 754 representable f64 constant");
                                            if normal.norm() > epsilon {
                                                solver.normals[idx] = normal.normalize();
                                            } else {
                                                solver.normals[idx] =
                                                    Vector3::new(T::zero(), T::zero(), T::one());
                                            }
                                        }
                                    }
                                } else {
                                    solver.normals[idx] = Vector3::zeros();
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// Determine whether the local interface is graph-like enough to justify
    /// a directional height function.
    fn hybrid_selection<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        self,
        solver: &VofSolver<T>,
        i: usize,
        j: usize,
        k: usize,
    ) -> HybridSelection<T> {
        let reference = self.calculate_gradient(solver, i, j, k);
        let epsilon = <T as FromPrimitive>::from_f64(VOF_EPSILON)
            .expect("VOF_EPSILON is an IEEE 754 representable f64 constant");
        let reference_norm = reference.norm();

        if reference_norm <= epsilon {
            return HybridSelection {
                reference,
                axis: None,
            };
        }

        let abs_components = [
            num_traits::Float::abs(reference.x),
            num_traits::Float::abs(reference.y),
            num_traits::Float::abs(reference.z),
        ];
        let mut axes = [1usize, 0, 2];
        axes.sort_by(|lhs, rhs| {
            abs_components[*rhs]
                .partial_cmp(&abs_components[*lhs])
                .unwrap_or(Ordering::Equal)
        });

        let axis_dominance = <T as FromPrimitive>::from_f64(HYBRID_AXIS_DOMINANCE)
            .expect("hybrid axis dominance threshold is representable in IEEE 754");
        let axis = if abs_components[axes[0]] / reference_norm >= axis_dominance {
            Some(axes[0])
        } else {
            None
        };

        HybridSelection { reference, axis }
    }

    /// Calculate interface normal from the volume-fraction gradient using
    /// Youngs' mixed finite-difference stencil.
    fn calculate_gradient<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        self,
        solver: &VofSolver<T>,
        i: usize,
        j: usize,
        k: usize,
    ) -> Vector3<T> {
        let two = <T as FromPrimitive>::from_f64(2.0)
            .expect("2.0 is representable in all IEEE 754 types");

        let dx = (solver.alpha[solver.index(i + 1, j, k)]
            - solver.alpha[solver.index(i - 1, j, k)])
            / (two * solver.dx);
        let dy = (solver.alpha[solver.index(i, j + 1, k)]
            - solver.alpha[solver.index(i, j - 1, k)])
            / (two * solver.dy);
        let dz = (solver.alpha[solver.index(i, j, k + 1)]
            - solver.alpha[solver.index(i, j, k - 1)])
            / (two * solver.dz);

        // The interface normal points from fluid 1 to fluid 0, which is the
        // negative volume-fraction gradient for the convention used here.
        Vector3::new(-dx, -dy, -dz)
    }

    /// Hybrid normal reconstruction: prefer a directional height function when
    /// the local gradient is graph-like, otherwise fall back to Youngs.
    fn hybrid_normal<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        self,
        solver: &VofSolver<T>,
        cache: &DirectionalHeightCache<T>,
        i: usize,
        j: usize,
        k: usize,
    ) -> Vector3<T> {
        let selection = self.hybrid_selection(solver, i, j, k);
        let Some(axis) = selection.axis else {
            return selection.reference;
        };

        let epsilon = <T as FromPrimitive>::from_f64(VOF_EPSILON)
            .expect("VOF_EPSILON is an IEEE 754 representable f64 constant");
        let candidate = Self::directional_height_normal(solver, cache, axis, i, j, k);
        if candidate.norm() <= epsilon {
            return selection.reference;
        }

        if candidate.dot(&selection.reference) < T::zero() {
            -candidate
        } else {
            candidate
        }
    }

    /// Directional height-function normal using a cached column integral.
    fn directional_height_normal<
        T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
    >(
        solver: &VofSolver<T>,
        cache: &DirectionalHeightCache<T>,
        axis: usize,
        i: usize,
        j: usize,
        k: usize,
    ) -> Vector3<T> {
        let two = <T as FromPrimitive>::from_f64(2.0)
            .expect("2.0 is representable in all IEEE 754 types");

        match axis {
            0 => {
                let jm = j - 1;
                let jp = j + 1;
                let km = k - 1;
                let kp = k + 1;
                let h_jm = cache.x[k * solver.ny + jm];
                let h_jp = cache.x[k * solver.ny + jp];
                let h_km = cache.x[km * solver.ny + j];
                let h_kp = cache.x[kp * solver.ny + j];
                let dh_dy = -(h_jp - h_jm) / (two * solver.dy);
                let dh_dz = -(h_kp - h_km) / (two * solver.dz);
                Vector3::new(T::one(), dh_dy, dh_dz)
            }
            1 => {
                let im = i - 1;
                let ip = i + 1;
                let km = k - 1;
                let kp = k + 1;
                let h_im = cache.y[k * solver.nx + im];
                let h_ip = cache.y[k * solver.nx + ip];
                let h_km = cache.y[km * solver.nx + i];
                let h_kp = cache.y[kp * solver.nx + i];
                let dh_dx = -(h_ip - h_im) / (two * solver.dx);
                let dh_dz = -(h_kp - h_km) / (two * solver.dz);
                Vector3::new(dh_dx, T::one(), dh_dz)
            }
            _ => {
                let im = i - 1;
                let ip = i + 1;
                let jm = j - 1;
                let jp = j + 1;
                let h_im = cache.z[j * solver.nx + im];
                let h_ip = cache.z[j * solver.nx + ip];
                let h_jm = cache.z[jm * solver.nx + i];
                let h_jp = cache.z[jp * solver.nx + i];
                let dh_dx = -(h_ip - h_im) / (two * solver.dx);
                let dh_dy = -(h_jp - h_jm) / (two * solver.dy);
                Vector3::new(dh_dx, dh_dy, T::one())
            }
        }
    }

    /// PLIC reconstruction using Youngs' gradient normal and the Scardovelli-Zaleski
    /// analytical volume formula.
    ///
    /// ## Algorithm
    ///
    /// 1. Estimate interface normal from the mixed-difference gradient of α.
    /// 2. Bisect the monotone function `V(C)` to find `C` such that
    ///    `V(n̂, C, dx, dy, dz) = α_i dx dy dz`.
    ///
    /// **Volume formula**: The exact formula of Scardovelli & Zaleski (2000)
    /// Eqs. 2.34–2.38 is used, implemented in `volume_under_plane_3d`.
    fn plic_reconstruction<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        self,
        solver: &VofSolver<T>,
        cache: &DirectionalHeightCache<T>,
        i: usize,
        j: usize,
        k: usize,
    ) -> (Vector3<T>, T) {
        // 1. Interface normal from a dominant-axis height function when it is
        // graph-like, otherwise use the Youngs gradient fallback.
        let mut normal = self.hybrid_normal(solver, cache, i, j, k);
        let epsilon = <T as FromPrimitive>::from_f64(VOF_EPSILON)
            .expect("VOF_EPSILON is an IEEE 754 representable f64 constant");

        if normal.norm() > epsilon {
            normal = normal.normalize();
        } else {
            normal = Vector3::new(T::zero(), T::zero(), T::one());
        }

        // 2. Find plane constant that conserves volume
        let target = solver.alpha[solver.index(i, j, k)];
        let plane_constant =
            self.find_plane_constant(normal, target, solver.dx, solver.dy, solver.dz);

        (normal, plane_constant)
    }

    /// Directional height-function curvature reconstruction.
    ///
    /// Returns the curvature together with the unnormalised directional
    /// height normal so the caller can preserve the outward orientation.
    fn directional_height_curvature<
        T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
    >(
        solver: &VofSolver<T>,
        cache: &DirectionalHeightCache<T>,
        axis: usize,
        i: usize,
        j: usize,
        k: usize,
    ) -> (Vector3<T>, T) {
        let two = <T as FromPrimitive>::from_f64(2.0)
            .expect("2.0 is representable in all IEEE 754 types");
        let four = <T as FromPrimitive>::from_f64(4.0)
            .expect("4.0 is representable in all IEEE 754 types");

        match axis {
            0 => {
                let h = cache.x[k * solver.ny + j];
                let h_jm = cache.x[k * solver.ny + (j - 1)];
                let h_jp = cache.x[k * solver.ny + (j + 1)];
                let h_km = cache.x[(k - 1) * solver.ny + j];
                let h_kp = cache.x[(k + 1) * solver.ny + j];
                let h_jm_km = cache.x[(k - 1) * solver.ny + (j - 1)];
                let h_jm_kp = cache.x[(k + 1) * solver.ny + (j - 1)];
                let h_jp_km = cache.x[(k - 1) * solver.ny + (j + 1)];
                let h_jp_kp = cache.x[(k + 1) * solver.ny + (j + 1)];

                let p = (h_jp - h_jm) / (two * solver.dy);
                let q = (h_kp - h_km) / (two * solver.dz);
                let p_yy = (h_jp - two * h + h_jm) / (solver.dy * solver.dy);
                let q_zz = (h_kp - two * h + h_km) / (solver.dz * solver.dz);
                let p_yz = (h_jp_kp - h_jp_km - h_jm_kp + h_jm_km) / (four * solver.dy * solver.dz);
                let base = T::one() + p * p + q * q;
                let denom = base * num_traits::Float::sqrt(base);
                let curvature = ((T::one() + q * q) * p_yy - two * p * q * p_yz
                    + (T::one() + p * p) * q_zz)
                    / denom;

                (Vector3::new(T::one(), -p, -q), curvature)
            }
            1 => {
                let h = cache.y[k * solver.nx + i];
                let h_im = cache.y[k * solver.nx + (i - 1)];
                let h_ip = cache.y[k * solver.nx + (i + 1)];
                let h_km = cache.y[(k - 1) * solver.nx + i];
                let h_kp = cache.y[(k + 1) * solver.nx + i];
                let h_im_km = cache.y[(k - 1) * solver.nx + (i - 1)];
                let h_im_kp = cache.y[(k + 1) * solver.nx + (i - 1)];
                let h_ip_km = cache.y[(k - 1) * solver.nx + (i + 1)];
                let h_ip_kp = cache.y[(k + 1) * solver.nx + (i + 1)];

                let p = (h_ip - h_im) / (two * solver.dx);
                let q = (h_kp - h_km) / (two * solver.dz);
                let p_xx = (h_ip - two * h + h_im) / (solver.dx * solver.dx);
                let q_zz = (h_kp - two * h + h_km) / (solver.dz * solver.dz);
                let p_xz = (h_ip_kp - h_ip_km - h_im_kp + h_im_km) / (four * solver.dx * solver.dz);
                let base = T::one() + p * p + q * q;
                let denom = base * num_traits::Float::sqrt(base);
                let curvature = ((T::one() + q * q) * p_xx - two * p * q * p_xz
                    + (T::one() + p * p) * q_zz)
                    / denom;

                (Vector3::new(-p, T::one(), -q), curvature)
            }
            _ => {
                let h = cache.z[j * solver.nx + i];
                let h_im = cache.z[j * solver.nx + (i - 1)];
                let h_ip = cache.z[j * solver.nx + (i + 1)];
                let h_jm = cache.z[(j - 1) * solver.nx + i];
                let h_jp = cache.z[(j + 1) * solver.nx + i];
                let h_im_jm = cache.z[(j - 1) * solver.nx + (i - 1)];
                let h_im_jp = cache.z[(j + 1) * solver.nx + (i - 1)];
                let h_ip_jm = cache.z[(j - 1) * solver.nx + (i + 1)];
                let h_ip_jp = cache.z[(j + 1) * solver.nx + (i + 1)];

                let p = (h_ip - h_im) / (two * solver.dx);
                let q = (h_jp - h_jm) / (two * solver.dy);
                let p_xx = (h_ip - two * h + h_im) / (solver.dx * solver.dx);
                let q_yy = (h_jp - two * h + h_jm) / (solver.dy * solver.dy);
                let p_xy = (h_ip_jp - h_ip_jm - h_im_jp + h_im_jm) / (four * solver.dx * solver.dy);
                let base = T::one() + p * p + q * q;
                let denom = base * num_traits::Float::sqrt(base);
                let curvature = ((T::one() + q * q) * p_xx - two * p * q * p_xy
                    + (T::one() + p * p) * q_yy)
                    / denom;

                (Vector3::new(-p, -q, T::one()), curvature)
            }
        }
    }

    /// Find the PLIC plane constant C by bisection such that
    /// `V_fluid(n, C, dx, dy, dz) = alpha * dx * dy * dz`.
    ///
    /// The bisection terminates when the volume residual is smaller than
    /// `PLIC_TOLERANCE * dx * dy * dz`, with a finite-precision guard that
    /// returns the last midpoint if the interval stops shrinking.
    fn find_plane_constant<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        self,
        normal: Vector3<T>,
        target_volume: T,
        dx: T,
        dy: T,
        dz: T,
    ) -> T {
        use num_traits::Float;

        let cell_volume = dx * dy * dz;
        let mut c_min = T::zero();
        let mut c_max =
            Float::abs(normal.x) * dx + Float::abs(normal.y) * dy + Float::abs(normal.z) * dz;

        let tolerance = <T as FromPrimitive>::from_f64(constants::PLIC_TOLERANCE)
            .expect("PLIC_TOLERANCE is an IEEE 754 representable f64 constant");
        let volume_tolerance = tolerance * cell_volume;
        let half =
            <T as FromPrimitive>::from_f64(0.5).expect("0.5 is exactly representable in IEEE 754");

        loop {
            let c_mid = c_min + (c_max - c_min) * half;
            let volume = volume_under_plane_3d(normal, c_mid, dx, dy, dz);

            if Float::abs(volume - target_volume * cell_volume) < volume_tolerance {
                return c_mid;
            }

            if c_mid == c_min || c_mid == c_max {
                return c_mid;
            }

            if volume < target_volume * cell_volume {
                c_min = c_mid;
            } else {
                c_max = c_mid;
            }
        }
    }

    /// Calculate interface curvature from the hybrid normal field.
    ///
    /// Graph-like interface cells reuse the precomputed outward normal field
    /// to select the dominant height-function axis; steep or ambiguous cells
    /// fall back to the divergence of the pre-computed normal field:
    ///
    /// `κ = −∇·n̂ = −(∂n̂_x/∂x + ∂n̂_y/∂y + ∂n̂_z/∂z)`
    fn calculate_curvature<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        self,
        solver: &mut VofSolver<T>,
        height_cache: Option<&DirectionalHeightCache<T>>,
    ) {
        let two = <T as FromPrimitive>::from_f64(2.0)
            .expect("2.0 is representable in all IEEE 754 types");
        let interface_lower = <T as FromPrimitive>::from_f64(VOF_INTERFACE_LOWER)
            .expect("VOF_INTERFACE_LOWER is an IEEE 754 representable f64 constant");
        let interface_upper = <T as FromPrimitive>::from_f64(VOF_INTERFACE_UPPER)
            .expect("VOF_INTERFACE_UPPER is an IEEE 754 representable f64 constant");

        for k_block in (1..solver.nz - 1).step_by(CACHE_BLOCK_SIZE_K) {
            for j_block in (1..solver.ny - 1).step_by(CACHE_BLOCK_SIZE_J) {
                for i_block in (1..solver.nx - 1).step_by(CACHE_BLOCK_SIZE_I) {
                    let k_end = (k_block + CACHE_BLOCK_SIZE_K).min(solver.nz - 1);
                    let j_end = (j_block + CACHE_BLOCK_SIZE_J).min(solver.ny - 1);
                    let i_end = (i_block + CACHE_BLOCK_SIZE_I).min(solver.nx - 1);

                    for k in k_block..k_end {
                        for j in j_block..j_end {
                            for i in i_block..i_end {
                                let idx = solver.index(i, j, k);
                                let alpha = solver.alpha[idx];

                                if alpha > interface_lower && alpha < interface_upper {
                                    if let Some(cache) = height_cache {
                                        let normal = solver.normals[idx];
                                        if let Some(axis) = Self::dominant_axis_from_normal(&normal)
                                        {
                                            let (candidate, curvature) =
                                                Self::directional_height_curvature(
                                                    solver, cache, axis, i, j, k,
                                                );
                                            solver.curvature[idx] =
                                                if candidate.dot(&normal) < T::zero() {
                                                    -curvature
                                                } else {
                                                    curvature
                                                };
                                            continue;
                                        }
                                    }

                                    let idx_xm = solver.index(i - 1, j, k);
                                    let idx_xp = solver.index(i + 1, j, k);
                                    let idx_ym = solver.index(i, j - 1, k);
                                    let idx_yp = solver.index(i, j + 1, k);
                                    let idx_zm = solver.index(i, j, k - 1);
                                    let idx_zp = solver.index(i, j, k + 1);

                                    let dnx_dx = (solver.normals[idx_xp].x
                                        - solver.normals[idx_xm].x)
                                        / (two * solver.dx);
                                    let dny_dy = (solver.normals[idx_yp].y
                                        - solver.normals[idx_ym].y)
                                        / (two * solver.dy);
                                    let dnz_dz = (solver.normals[idx_zp].z
                                        - solver.normals[idx_zm].z)
                                        / (two * solver.dz);

                                    solver.curvature[idx] = -(dnx_dx + dny_dy + dnz_dz);
                                } else {
                                    solver.curvature[idx] = T::zero();
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// Determine whether a precomputed normal is graph-like enough to justify
    /// a directional height-function curvature estimate.
    fn dominant_axis_from_normal<
        T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
    >(
        normal: &Vector3<T>,
    ) -> Option<usize> {
        let abs_components = [
            num_traits::Float::abs(normal.x),
            num_traits::Float::abs(normal.y),
            num_traits::Float::abs(normal.z),
        ];
        let mut axes = [1usize, 0, 2];
        axes.sort_by(|lhs, rhs| {
            abs_components[*rhs]
                .partial_cmp(&abs_components[*lhs])
                .unwrap_or(Ordering::Equal)
        });

        let axis_dominance = <T as FromPrimitive>::from_f64(HYBRID_AXIS_DOMINANCE)
            .expect("hybrid axis dominance threshold is representable in IEEE 754");
        if abs_components[axes[0]] >= axis_dominance {
            Some(axes[0])
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Horizontal interface (α=1 below midpoint, α=0 above) should yield
    /// a vertical normal [0, 1].
    #[test]
    fn test_height_function_flat_interface() {
        let nx = 10;
        let ny = 10;
        let dx = 1.0;
        let dy = 1.0;

        // Fill: bottom half = 1.0, top half = 0.0
        let alpha: Vec<Vec<f64>> = (0..nx)
            .map(|_| {
                (0..ny)
                    .map(|j| if j < ny / 2 { 1.0 } else { 0.0 })
                    .collect()
            })
            .collect();

        // At the interface row, the normal should point upward: [0, 1]
        let n = height_function_normal_2d(&alpha, 5, ny / 2, dx, dy);
        assert!(
            n[0].abs() < 1e-10,
            "n_x should be ~0 for flat interface, got {}",
            n[0]
        );
        assert!(
            n[1] > 0.99,
            "n_y should be ~1 for flat interface, got {}",
            n[1]
        );
    }

    /// Youngs-gradient reconstruction should use the same outward normal
    /// convention as the height-function path.
    #[test]
    fn test_youngs_gradient_flat_interface_points_upward() {
        let nx = 10;
        let ny = 10;
        let dx = 1.0;
        let dy = 1.0;

        let alpha: Vec<Vec<f64>> = (0..nx)
            .map(|_| {
                (0..ny)
                    .map(|j| if j < ny / 2 { 1.0 } else { 0.0 })
                    .collect()
            })
            .collect();

        let n = youngs_normal_2d(&alpha, 5, ny / 2, dx, dy);
        assert!(n[0].abs() < 1e-10, "n_x should be ~0 for flat interface");
        assert!(
            n[1] > 0.99,
            "n_y should be ~1 for flat interface, got {}",
            n[1]
        );
    }

    /// The solver-facing gradient reconstruction should match the helper
    /// orientation and point from fluid 1 toward fluid 0.
    #[test]
    fn test_gradient_reconstruction_flat_interface_points_upward() {
        let nx = 10;
        let ny = 10;
        let nz = 3;
        let dx = 1.0;
        let dy = 1.0;
        let dz = 1.0;

        let config = VofConfig {
            reconstruction_method: InterfaceReconstruction::Gradient,
            ..VofConfig::default()
        };
        let mut solver = VofSolver::create(config, nx, ny, nz, dx, dy, dz);

        let mut alpha = vec![0.0; nx * ny * nz];
        for k in 0..nz {
            for j in 0..ny {
                let value = if j < ny / 2 {
                    1.0
                } else if j == ny / 2 {
                    0.5
                } else {
                    0.0
                };
                for i in 0..nx {
                    alpha[k * ny * nx + j * nx + i] = value;
                }
            }
        }
        solver
            .set_volume_fraction(alpha)
            .expect("volume fraction field should match solver dimensions");
        solver.reconstruct_interface();

        let idx = solver.linear_index(5, ny / 2, 1);
        let normal: nalgebra::Vector3<f64> = solver.normals()[idx];
        assert!(
            normal.x.abs() < 1e-10_f64,
            "n_x should be ~0 for flat interface"
        );
        assert!(
            normal.y > 0.99,
            "n_y should be ~1 for flat interface, got {}",
            normal.y
        );
    }

    /// PLIC reconstruction should select the dominant-axis height function
    /// for a vertical interface and keep the outward normal orientation.
    #[test]
    fn test_plic_reconstruction_vertical_interface_points_right() {
        let nx = 10;
        let ny = 8;
        let nz = 6;
        let dx = 1.0;
        let dy = 1.0;
        let dz = 1.0;

        let config = VofConfig {
            reconstruction_method: InterfaceReconstruction::PLIC,
            ..VofConfig::default()
        };
        let mut solver = VofSolver::create(config, nx, ny, nz, dx, dy, dz);

        let mut alpha = vec![0.0; nx * ny * nz];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let value = if i < nx / 2 {
                        1.0
                    } else if i == nx / 2 {
                        0.5
                    } else {
                        0.0
                    };
                    alpha[k * ny * nx + j * nx + i] = value;
                }
            }
        }

        solver
            .set_volume_fraction(alpha)
            .expect("volume fraction field should match solver dimensions");
        solver.reconstruct_interface();

        let idx = solver.linear_index(nx / 2, ny / 2, nz / 2);
        let normal: nalgebra::Vector3<f64> = solver.normals()[idx];
        assert!(
            normal.x > 0.99,
            "n_x should be ~1 for the vertical interface, got {}",
            normal.x
        );
        assert!(
            normal.y.abs() < 0.05,
            "n_y should be ~0 for the vertical interface, got {}",
            normal.y
        );
        assert!(
            normal.z.abs() < 0.05,
            "n_z should be ~0 for the vertical interface, got {}",
            normal.z
        );
        assert!(
            solver.curvature()[idx].abs() < 1e-10_f64,
            "curvature should be ~0 for a planar vertical interface, got {}",
            solver.curvature()[idx]
        );
    }

    /// A diagonal 3D interface should fail the dominant-axis test and fall
    /// back to the Youngs gradient path.
    #[test]
    fn test_plic_reconstruction_diagonal_interface_uses_gradient_fallback() {
        let nx = 6;
        let ny = 6;
        let nz = 6;
        let dx = 1.0;
        let dy = 1.0;
        let dz = 1.0;

        let config = VofConfig {
            reconstruction_method: InterfaceReconstruction::PLIC,
            ..VofConfig::default()
        };
        let mut solver = VofSolver::create(config, nx, ny, nz, dx, dy, dz);

        let mut alpha = vec![0.0; nx * ny * nz];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let sum = i + j + k;
                    let value = if sum < 9 {
                        1.0
                    } else if sum == 9 {
                        0.5
                    } else {
                        0.0
                    };
                    alpha[k * ny * nx + j * nx + i] = value;
                }
            }
        }

        solver
            .set_volume_fraction(alpha)
            .expect("volume fraction field should match solver dimensions");
        solver.reconstruct_interface();

        let idx = solver.linear_index(3, 3, 3);
        let selection = InterfaceReconstruction::PLIC.hybrid_selection(&solver, 3, 3, 3);
        assert!(
            selection.axis.is_none(),
            "diagonal interface should not be treated as graph-like"
        );

        let normal: nalgebra::Vector3<f64> = solver.normals()[idx];
        let expected = nalgebra::Vector3::new(1.0, 1.0, 1.0).normalize();
        assert!(
            (normal - expected).norm() < 0.05,
            "fallback normal should track the Youngs gradient, got {} vs {}",
            normal,
            expected
        );
        assert!(
            solver.curvature()[idx].abs() < 1e-10_f64,
            "planar diagonal interface should have ~zero curvature, got {}",
            solver.curvature()[idx]
        );
    }

    /// 45-degree tilted interface should give a normal close to
    /// [-0.707, 0.707] (pointing from the filled region toward the empty region).
    #[test]
    fn test_height_function_tilted_interface() {
        let nx = 20;
        let ny = 20;
        let dx = 1.0;
        let dy = 1.0;

        // Fill: α = 1 where j < i (below a 45° line from bottom-left)
        // This creates a tilted interface where column height H_i = i * dy
        let alpha: Vec<Vec<f64>> = (0..nx)
            .map(|i| (0..ny).map(|j| if j < i { 1.0 } else { 0.0 }).collect())
            .collect();

        // At a central cell the height-function should detect the tilt.
        // H_i = i * dy, so dH/dx = dy/dx = 1.0
        // n_x = -1.0, n_y = 1.0, normalised: [-0.707, 0.707]
        let n = height_function_normal_2d(&alpha, 10, 10, dx, dy);
        let expected_nx = -1.0_f64 / 2.0_f64.sqrt();
        let expected_ny = 1.0_f64 / 2.0_f64.sqrt();

        assert!(
            (n[0] - expected_nx).abs() < 0.05,
            "n_x should be ~{:.4}, got {:.4}",
            expected_nx,
            n[0]
        );
        assert!(
            (n[1] - expected_ny).abs() < 0.05,
            "n_y should be ~{:.4}, got {:.4}",
            expected_ny,
            n[1]
        );
    }

    /// For a smooth graph-like interface, the height-function normal should
    /// have lower angular error than the Youngs gradient normal.
    #[test]
    fn test_height_function_vs_youngs_accuracy() {
        let nx = 48;
        let ny = 40;
        let dx = 1.0;
        let dy = 1.0;

        // Smooth single-valued interface: y = f(x)
        let interface = |x: f64| 19.0 + 3.5 * (x / 7.5).sin();
        let interface_derivative = |x: f64| (3.5 / 7.5) * (x / 7.5).cos();

        let alpha: Vec<Vec<f64>> = (0..nx)
            .map(|i| {
                let x = i as f64 + 0.5;
                let y0 = interface(x);
                (0..ny)
                    .map(|j| {
                        let y = j as f64 + 0.5;
                        if y < y0 - 0.5 {
                            1.0
                        } else if y > y0 + 0.5 {
                            0.0
                        } else {
                            // Linear sub-cell interpolation for interface cells.
                            (y0 + 0.5 - y).clamp(0.0, 1.0)
                        }
                    })
                    .collect()
            })
            .collect();

        let mut hf_total_error = 0.0;
        let mut youngs_total_error = 0.0;
        let mut count = 0u32;

        // Sample interface cells (those near the graph interface).
        for i in 2..nx - 2 {
            for j in 2..ny - 2 {
                let a = alpha[i][j];
                if a > 0.05 && a < 0.95 {
                    let x = i as f64 + 0.5;
                    let slope = interface_derivative(x);
                    let norm = (1.0 + slope * slope).sqrt();
                    let exact = [-slope / norm, 1.0 / norm];

                    let hf_n = height_function_normal_2d(&alpha, i, j, dx, dy);
                    let y_n = youngs_normal_2d(&alpha, i, j, dx, dy);

                    // Angular error (dot product → angle)
                    let hf_dot = (hf_n[0] * exact[0] + hf_n[1] * exact[1]).clamp(-1.0, 1.0);
                    let y_dot = (y_n[0] * exact[0] + y_n[1] * exact[1]).clamp(-1.0, 1.0);

                    hf_total_error += hf_dot.acos();
                    youngs_total_error += y_dot.acos();
                    count += 1;
                }
            }
        }

        assert!(
            count > 10,
            "Should have enough interface cells, got {count}"
        );

        let hf_mean_error = hf_total_error / count as f64;
        let youngs_mean_error = youngs_total_error / count as f64;

        assert!(
            hf_mean_error < youngs_mean_error,
            "Height-function mean angular error ({:.4} rad) should be less than \
             Youngs gradient error ({:.4} rad)",
            hf_mean_error,
            youngs_mean_error
        );
    }

    /// PLIC plane bisection must scale its stopping tolerance with the cell
    /// size used in the geometric inversion.
    #[test]
    fn test_find_plane_constant_scales_with_cell_size() {
        let normal = Vector3::new(0.83_f64, 0.31, 0.46).normalize();
        let dx = 1.0e-4_f64;
        let dy = 2.0e-4_f64;
        let dz = 4.0e-4_f64;
        let target_volume_fraction = 0.25_f64;

        let plane_constant = InterfaceReconstruction::PLIC.find_plane_constant(
            normal,
            target_volume_fraction,
            dx,
            dy,
            dz,
        );
        let volume = volume_under_plane_3d(normal, plane_constant, dx, dy, dz);
        let target_volume = target_volume_fraction * dx * dy * dz;
        let tolerance = constants::PLIC_TOLERANCE * dx * dy * dz;

        assert!(
            (volume - target_volume).abs() <= tolerance,
            "bisection volume error must scale with the cell size"
        );
        assert!(plane_constant >= 0.0);
        assert!(plane_constant <= normal.x.abs() * dx + normal.y.abs() * dy + normal.z.abs() * dz);
    }
}
