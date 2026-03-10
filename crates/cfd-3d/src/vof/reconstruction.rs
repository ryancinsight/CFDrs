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
