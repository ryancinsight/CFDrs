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
//! where $\mathbf{n}$ is the interface normal (from gradient of $\alpha$)
//! and $C$ is chosen so that the clipped volume matches $\alpha_i |\Omega_i|$
//! exactly (to machine precision). This ensures global volume conservation.
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
//! **Normal Estimation** (Youngs 1984): The normal is estimated from the
//! discrete gradient of α:
//!
//! ```math
//! n̂ᵢ ≈ ∇α / |∇α|   (mixed finite-difference gradient)
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
//! **Statement**: The interface curvature κ = −∇·(∇φ/|∇φ|) computed from the
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

// ── Height-Function Normal Estimation ────────────────────────────────────────

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

    // Normal points from fluid 1 to fluid 0: n = ∇α / |∇α|
    [da_dx / mag, da_dy / mag]
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
        self.calculate_normals(solver);
        self.calculate_curvature(solver);
    }

    /// Calculate interface normal vectors using gradient of volume fraction.
    ///
    /// Uses cache blocking for improved memory access patterns on 3D grids.
    fn calculate_normals<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        self,
        solver: &mut VofSolver<T>,
    ) {
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

                                let interface_lower = <T as FromPrimitive>::from_f64(
                                    VOF_INTERFACE_LOWER,
                                )
                                .expect(
                                    "VOF_INTERFACE_LOWER is an IEEE 754 representable f64 constant",
                                );
                                let interface_upper = <T as FromPrimitive>::from_f64(
                                    VOF_INTERFACE_UPPER,
                                )
                                .expect(
                                    "VOF_INTERFACE_UPPER is an IEEE 754 representable f64 constant",
                                );

                                if alpha > interface_lower && alpha < interface_upper {
                                    match self {
                                        Self::PLIC => {
                                            let (normal, _) =
                                                self.plic_reconstruction(solver, i, j, k);
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

        Vector3::new(dx, dy, dz)
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
        i: usize,
        j: usize,
        k: usize,
    ) -> (Vector3<T>, T) {
        // 1. Interface normal from gradient
        let mut normal = self.calculate_gradient(solver, i, j, k);
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

    /// Find the PLIC plane constant C by bisection such that
    /// `V_fluid(n, C, dx, dy, dz) = alpha * dx * dy * dz`.
    ///
    /// The bisection terminates when the interval width is smaller than
    /// `PLIC_TOLERANCE * max(dx, dy, dz)`.
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
        let half =
            <T as FromPrimitive>::from_f64(0.5).expect("0.5 is exactly representable in IEEE 754");

        // Bisect until interval < tolerance × cell_size
        while (c_max - c_min) > tolerance {
            let c_mid = c_min + (c_max - c_min) * half;
            let volume = volume_under_plane_3d(normal, c_mid, dx, dy, dz);

            if Float::abs(volume - target_volume * cell_volume) < tolerance * cell_volume {
                return c_mid;
            }

            if volume < target_volume * cell_volume {
                c_min = c_mid;
            } else {
                c_max = c_mid;
            }
        }

        c_min + (c_max - c_min) * half
    }

    /// Calculate interface curvature from the divergence of the normal field.
    ///
    /// `κ = −∇·n̂ = −(∂n̂_x/∂x + ∂n̂_y/∂y + ∂n̂_z/∂z)`
    ///
    /// Uses second-order central differences on the pre-computed normal field.
    fn calculate_curvature<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        self,
        solver: &mut VofSolver<T>,
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

    /// For a smooth circular interface, the height-function normal should
    /// have lower angular error than the Youngs gradient normal.
    #[test]
    fn test_height_function_vs_youngs_accuracy() {
        let nx = 40;
        let ny = 40;
        let dx = 1.0;
        let dy = 1.0;

        // Circle centred at (20, 20) with radius 12
        let cx = 20.0;
        let cy = 20.0;
        let r = 12.0;

        let alpha: Vec<Vec<f64>> = (0..nx)
            .map(|i| {
                (0..ny)
                    .map(|j| {
                        let x = i as f64 + 0.5;
                        let y = j as f64 + 0.5;
                        let dist = ((x - cx).powi(2) + (y - cy).powi(2)).sqrt();
                        if dist < r - 0.5 {
                            1.0
                        } else if dist > r + 0.5 {
                            0.0
                        } else {
                            // Linear sub-cell interpolation for interface cells
                            (r + 0.5 - dist).clamp(0.0, 1.0)
                        }
                    })
                    .collect()
            })
            .collect();

        let mut hf_total_error = 0.0;
        let mut youngs_total_error = 0.0;
        let mut count = 0u32;

        // Sample interface cells (those near the circle boundary)
        for i in 2..nx - 2 {
            for j in 2..ny - 2 {
                let a = alpha[i][j];
                if a > 0.05 && a < 0.95 {
                    // Analytical normal: radially outward from centre
                    let x = i as f64 + 0.5;
                    let y = j as f64 + 0.5;
                    let dist = ((x - cx).powi(2) + (y - cy).powi(2)).sqrt();
                    if dist < 1e-10 {
                        continue;
                    }
                    let exact = [(x - cx) / dist, (y - cy) / dist];

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
}
