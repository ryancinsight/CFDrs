//! Interface reconstruction methods for VOF

use super::config::{constants, VofConfig, VOF_EPSILON, VOF_INTERFACE_LOWER, VOF_INTERFACE_UPPER};
use super::solver::VofSolver;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

/// Cache blocking parameters for optimized 3D traversal
const CACHE_BLOCK_SIZE_I: usize = 8;
const CACHE_BLOCK_SIZE_J: usize = 8;
const CACHE_BLOCK_SIZE_K: usize = 8;

/// Interface reconstruction strategy
pub struct InterfaceReconstruction {
    use_plic: bool,
}

impl InterfaceReconstruction {
    /// Create reconstruction method based on configuration
    #[must_use] pub fn create(config: &VofConfig) -> Self {
        Self {
            use_plic: config.use_plic,
        }
    }

    /// Reconstruct interface normals and curvature
    pub fn reconstruct<T: RealField + FromPrimitive + Copy>(&self, solver: &mut VofSolver<T>) {
        // Calculate interface normals
        self.calculate_normals(solver);

        // Calculate interface curvature
        self.calculate_curvature(solver);
    }

    /// Calculate interface normal vectors using gradient of volume fraction
    /// Uses cache blocking for improved memory access patterns
    fn calculate_normals<T: RealField + FromPrimitive + Copy>(&self, solver: &mut VofSolver<T>) {
        // Process domain in cache-friendly blocks
        for k_block in (1..solver.nz - 1).step_by(CACHE_BLOCK_SIZE_K) {
            for j_block in (1..solver.ny - 1).step_by(CACHE_BLOCK_SIZE_J) {
                for i_block in (1..solver.nx - 1).step_by(CACHE_BLOCK_SIZE_I) {
                    // Process all cells within the block
                    let k_end = (k_block + CACHE_BLOCK_SIZE_K).min(solver.nz - 1);
                    let j_end = (j_block + CACHE_BLOCK_SIZE_J).min(solver.ny - 1);
                    let i_end = (i_block + CACHE_BLOCK_SIZE_I).min(solver.nx - 1);

                    for k in k_block..k_end {
                        for j in j_block..j_end {
                            for i in i_block..i_end {
                                let idx = solver.index(i, j, k);

                                // Only calculate for interface cells
                                let alpha = solver.alpha[idx];
                                let interface_lower = T::from_f64(VOF_INTERFACE_LOWER)
                                    .expect("Failed to represent VOF_INTERFACE_LOWER constant");
                                let interface_upper = T::from_f64(VOF_INTERFACE_UPPER)
                                    .expect("Failed to represent VOF_INTERFACE_UPPER constant");

                                if alpha > interface_lower && alpha < interface_upper {
                                    if self.use_plic {
                                        let (normal, _) = self.plic_reconstruction(solver, i, j, k);
                                        solver.normals[idx] = normal;
                                    } else {
                                        // Simple gradient-based normal
                                        let normal = self.calculate_gradient(solver, i, j, k);
                                        let epsilon = T::from_f64(VOF_EPSILON)
                                            .expect("Failed to represent VOF_EPSILON constant");
                                        if normal.norm() > epsilon {
                                            solver.normals[idx] = normal.normalize();
                                        } else {
                                            solver.normals[idx] =
                                                Vector3::new(T::zero(), T::zero(), T::one());
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

    /// Calculate gradient of volume fraction using central differences
    fn calculate_gradient<T: RealField + FromPrimitive + Copy>(
        &self,
        solver: &VofSolver<T>,
        i: usize,
        j: usize,
        k: usize,
    ) -> Vector3<T> {
        let two = T::from_f64(2.0).expect("Failed to represent constant 2.0");

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

    /// PLIC reconstruction of interface
    ///
    /// This is a single-step PLIC reconstruction using Youngs' method for the normal.
    /// The normal is calculated from the gradient of the volume fraction field.
    /// Note: This is NOT an iterative Newton-Raphson solver - it's a direct calculation.
    fn plic_reconstruction<T: RealField + FromPrimitive + Copy>(
        &self,
        solver: &VofSolver<T>,
        i: usize,
        j: usize,
        k: usize,
    ) -> (Vector3<T>, T) {
        // 1. Calculate the interface normal using Youngs' gradient method
        let mut normal = self.calculate_gradient(solver, i, j, k);
        let epsilon = T::from_f64(VOF_EPSILON).expect("Failed to represent VOF_EPSILON constant");

        if normal.norm() > epsilon {
            normal = normal.normalize();
        } else {
            // Default to vertical interface if gradient is zero
            normal = Vector3::new(T::zero(), T::zero(), T::one());
        }

        // 2. Find the plane constant that conserves the volume fraction
        let target_volume = solver.alpha[solver.index(i, j, k)];
        let plane_constant =
            self.find_plane_constant(normal, target_volume, solver.dx, solver.dy, solver.dz);

        (normal, plane_constant)
    }

    /// Find plane constant for given normal and target volume using binary search
    ///
    /// The search terminates when the interval width is smaller than the tolerance,
    /// not after a fixed number of iterations.
    fn find_plane_constant<T: RealField + FromPrimitive + Copy>(
        &self,
        normal: Vector3<T>,
        target_volume: T,
        dx: T,
        dy: T,
        dz: T,
    ) -> T {
        // Binary search for plane constant
        let cell_volume = dx * dy * dz;
        let mut c_min = T::zero();
        let mut c_max = normal.x.abs() * dx + normal.y.abs() * dy + normal.z.abs() * dz;

        let tolerance = T::from_f64(constants::PLIC_TOLERANCE)
            .expect("Failed to represent PLIC_TOLERANCE constant");
        let half = T::from_f64(0.5).expect("Failed to represent constant 0.5");

        // Loop until the search interval is smaller than the tolerance
        while (c_max - c_min) > tolerance {
            let c_mid = c_min + (c_max - c_min) * half;
            let volume = self.calculate_volume_under_plane_3d(normal, c_mid, dx, dy, dz);

            if (volume - target_volume * cell_volume).abs() < tolerance * cell_volume {
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

    /// Calculate volume of fluid under a plane in a 3D cell
    ///
    /// This implements the full 3D analytical formula from Scardovelli & Zaleski (2000).
    /// The formula requires sorting the normal components and using different expressions
    /// based on their relative magnitudes.
    fn calculate_volume_under_plane_3d<T: RealField + FromPrimitive + Copy>(
        &self,
        normal: Vector3<T>,
        plane_constant: T,
        dx: T,
        dy: T,
        dz: T,
    ) -> T {
        #[allow(clippy::no_effect_underscore_binding)] // Context variables for Scardovelli & Zaleski formula
        {
        let cell_volume = dx * dy * dz;

        // Normalize the normal vector components by cell dimensions
        let mut n = [
            normal.x.abs() * dx,
            normal.y.abs() * dy,
            normal.z.abs() * dz,
        ];

        // Sort components in ascending order
        n.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let n1 = n[0];
        let n2 = n[1];
        let n3 = n[2];

        // Check for degenerate case
        let epsilon = T::from_f64(1e-10).expect("Failed to represent epsilon");
        if n3 < epsilon {
            let half = T::from_f64(0.5).expect("Failed to represent 0.5");
            return half * cell_volume;
        }

        // Normalized plane constant
        let alpha = plane_constant / (n1 + n2 + n3);

        if alpha <= T::zero() {
            return T::zero();
        } else if alpha >= T::one() {
            return cell_volume;
        }

        // Constants for the analytical formula
        let six = T::from_f64(6.0).expect("Failed to represent 6.0");
        let two = T::from_f64(2.0).expect("Failed to represent 2.0");
        let three = T::from_f64(3.0).expect("Failed to represent 3.0");

        // Compute volume fraction based on the region
        // These are the analytical formulas from Scardovelli & Zaleski (2000)
        let m1 = n1 / n3;
        let m2 = n2 / n3;
        let m12 = m1 + m2;

        let volume_fraction = if alpha <= m1 {
            // Region 1: Small alpha, pyramid shape
            alpha * alpha * alpha / (six * m1 * m2)
        } else if alpha <= m2 {
            // Region 2: Intermediate, pentahedron
            let _alpha_m1 = alpha - m1;
            (alpha * alpha * (three * m2 - alpha) + m1 * m1 * (alpha - three * m2))
                / (six * m1 * m2)
        } else if alpha <= T::one() - m12 {
            // Region 3: Central region, hexahedron
            (alpha * alpha * (three - two * alpha)
                + m1 * m1 * (three * alpha - T::one())
                + m2 * m2 * (three * alpha - T::one()))
                / (six * m1 * m2)
        } else if alpha <= T::one() - m1 {
            // Region 4: Mirror of region 2
            let one_minus_alpha = T::one() - alpha;
            T::one()
                - (one_minus_alpha * one_minus_alpha * one_minus_alpha) / (six * m1 * m2)
                - (m1 - one_minus_alpha) * (m1 - one_minus_alpha) * (m1 - one_minus_alpha)
                    / (six * m1 * m2)
        } else {
            // Region 5: Large alpha, inverted pyramid
            let one_minus_alpha = T::one() - alpha;
            T::one() - (one_minus_alpha * one_minus_alpha * one_minus_alpha) / (six * m1 * m2)
        };

        volume_fraction * cell_volume
        }
    }

    /// Calculate interface curvature from normals using cache blocking
    fn calculate_curvature<T: RealField + FromPrimitive + Copy>(&self, solver: &mut VofSolver<T>) {
        let two = T::from_f64(2.0).expect("Failed to represent constant 2.0");
        let interface_lower = T::from_f64(VOF_INTERFACE_LOWER)
            .expect("Failed to represent VOF_INTERFACE_LOWER constant");
        let interface_upper = T::from_f64(VOF_INTERFACE_UPPER)
            .expect("Failed to represent VOF_INTERFACE_UPPER constant");

        // Process domain in cache-friendly blocks
        for k_block in (1..solver.nz - 1).step_by(CACHE_BLOCK_SIZE_K) {
            for j_block in (1..solver.ny - 1).step_by(CACHE_BLOCK_SIZE_J) {
                for i_block in (1..solver.nx - 1).step_by(CACHE_BLOCK_SIZE_I) {
                    // Process all cells within the block
                    let k_end = (k_block + CACHE_BLOCK_SIZE_K).min(solver.nz - 1);
                    let j_end = (j_block + CACHE_BLOCK_SIZE_J).min(solver.ny - 1);
                    let i_end = (i_block + CACHE_BLOCK_SIZE_I).min(solver.nx - 1);

                    for k in k_block..k_end {
                        for j in j_block..j_end {
                            for i in i_block..i_end {
                                let idx = solver.index(i, j, k);

                                // Only calculate for interface cells
                                let alpha = solver.alpha[idx];
                                if alpha > interface_lower && alpha < interface_upper {
                                    // Calculate divergence of normal (curvature = -∇·n)
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
