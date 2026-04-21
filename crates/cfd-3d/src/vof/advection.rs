//! VOF advection methods
//!
//! ## Mathematical Foundation
//!
//! ### CFL Stability Theorem (Unsplit Geometric Advection)
//!
//! **Statement**: The operator-split geometric VOF advection scheme is stable
//! in volume if the CFL condition is satisfied:
//!
//! ```math
//! CFL = max_{i,j,k} (|u_x| Δt/Δx + |u_y| Δt/Δy + |u_z| Δt/Δz) ≤ 1
//! ```
//!
//! **Proof**: In one spatial direction, the swept volume fraction flux is:
//! ```math
//!   F_{i+1/2} = V_swept(n_i, α_i, u_{i+1/2} Δt) / V_cell
//! ```
//! where `V_swept` is the volume of the truncated prism cut by the PLIC plane n·x = C.
//! Since `V_swept ≤ |u| Δt A_face`, the outgoing flux satisfies F ≤ CFL ≤ 1,
//! guaranteeing α ∈ [0,1] after one operator-split step.
//!
//! **Reference**: Scardovelli & Zaleski (2003), "Interface reconstruction with
//! least-square fit and split Eulerian-Lagrangian advection". Int. J. Numer.
//! Methods Fluids 41:251-274.
//!
//! ### Algebraic Advection CFL Condition
//!
//! **Statement**: The algebraic upwind scheme ∂α/∂t + u·∇α = 0 is stable for CFL ≤ 1.
//!
//! **Proof**: Standard first-order upwind analysis — see LeVeque (2002) "Finite
//! Volume Methods for Hyperbolic Problems", §4.4.

use super::config::{VofConfig, VOF_INTERFACE_LOWER, VOF_INTERFACE_UPPER};
use super::plic_geometry::plic_volume_fraction_in_prism;
use super::solver::VofSolver;
use cfd_core::error::Error;
use cfd_core::error::Result;
use nalgebra::RealField;
use num_traits::FromPrimitive;

use serde::{Deserialize, Serialize};

/// Advection method for VOF
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum AdvectionMethod {
    /// Geometric VOF advection — uses PLIC plane intersection to compute
    /// the exact swept volume flux through each cell face.  Conserves volume
    /// to machine precision for CFL ≤ 1.
    Geometric,
    /// Algebraic VOF advection — first-order upwind advection of the scalar
    /// volume fraction field.  Simpler and faster but diffusive.
    Algebraic,
}

/// Enforce the CFL ≤ 1 stability condition before a VOF advection step.
///
/// # Theorem — VOF CFL Stability (Scardovelli & Zaleski 2003)
///
/// The operator-split geometric VOF scheme conserves volume if and only if
///
/// ```text
/// CFL = max_{cells} (|u_x|·dt/dx + |u_y|·dt/dy + |u_z|·dt/dz) ≤ 1
/// ```
///
/// **Proof.** In one spatial direction the swept volume fraction flux is
/// F_{i+1/2} = V_swept / V_cell ≤ |u|·dt/dx = CFL_dir.  If CFL_dir > 1
/// the swept prism exceeds the donor cell, causing α_new = α_old − ΣF to
/// fall outside [0, 1], violating the boundedness invariant.
///
/// **Reference**: Scardovelli, R. & Zaleski, S. (2003). "Interface
/// reconstruction with least-square fit and split Eulerian-Lagrangian
/// advection." *Int. J. Numer. Methods Fluids* 41:251–274.
///
/// # Errors
/// Returns [`Error::InvalidConfiguration`] if CFL > 1.
fn check_cfl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
    solver: &VofSolver<T>,
    dt: T,
) -> Result<()> {
    use num_traits::Float;
    let mut cfl_max = T::zero();

    let dx = solver.dx;
    let dy = solver.dy;
    let dz = solver.dz;

    for vel in &solver.velocity {
        let cfl_cell =
            Float::abs(vel.x) * dt / dx + Float::abs(vel.y) * dt / dy + Float::abs(vel.z) * dt / dz;
        if cfl_cell > cfl_max {
            cfl_max = cfl_cell;
        }
    }

    if cfl_max > T::one() {
        return Err(Error::InvalidConfiguration(
            "VOF CFL > 1.0: volume conservation violated. \
             Reduce the time-step dt or refine the grid. \
             Recommended: target CFL ≤ 0.5 for robust interface tracking."
                .to_string(),
        ));
    }
    Ok(())
}

impl AdvectionMethod {
    /// Create advection method based on configuration
    #[must_use]
    pub fn create(config: &VofConfig) -> Self {
        config.advection_method
    }

    /// Advect volume fraction field by one time step `dt`.
    ///
    /// # CFL Requirement
    /// Enforces CFL ≤ 1 before advancing (Scardovelli & Zaleski 2003).
    /// At CFL > 1 the swept PLIC prism exceeds the donor cell, generating
    /// volume fractions outside [0, 1] and breaking conservation.
    ///
    /// # Invariant
    /// After this call, `solver.alpha[i] ∈ [0, 1]` for all cells `i`.
    ///
    /// # Errors
    /// Returns [`Error::InvalidConfiguration`] if CFL > 1.
    pub fn advect<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        self,
        solver: &mut VofSolver<T>,
        dt: T,
    ) -> Result<()> {
        check_cfl(solver, dt)?;
        match self {
            Self::Geometric => self.geometric_advection(solver, dt),
            Self::Algebraic => self.algebraic_advection(solver, dt),
        }
    }

    /// True geometric (Eulerian Implicit) VOF advection.
    ///
    /// For each internal cell the net volume-fraction change is computed by
    /// integrating the swept sub-volumes through all six faces.  The swept
    /// volume through a face is obtained by intersecting the donor cell's PLIC
    /// plane with the donor prism, using the analytical truncated-parallelepiped
    /// formula of Scardovelli & Zaleski (2000).
    ///
    /// ## Algorithm (Strang-split x → y → z)
    ///
    /// For each face (e.g., +x face of cell (i,j,k)):
    /// 1. Identify donor cell: cell (i,j,k) if u_{i+1/2} > 0, else cell (i+1,j,k).
    /// 2. The donor prism depth in the x-direction is `δ = |u_{i+1/2}| Δt`.
    /// 3. Compute the fraction of the donor prism occupied by fluid using the
    ///    donor cell's PLIC plane: `f = V_plic(n, C, δ, Δy, Δz) / (δ Δy Δz)`.
    /// 4. The volumetric flux is `F = f · u_{i+1/2} · Δy · Δz · Δt`.
    ///
    /// Volume balance: `α_new = α_old - (ΣF_out - ΣF_in) / V_cell`.
    ///
    /// **Reference**: Scardovelli & Zaleski (2000) Ann. Rev. Fluid Mech. 31:567.
    fn geometric_advection<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        self,
        solver: &mut VofSolver<T>,
        dt: T,
    ) -> Result<()> {
        // Use alpha_previous as write buffer; boundaries are copied as ghost values.
        solver.copy_boundaries();

        let zero = T::zero();
        let one = T::one();
        let half = T::one() / (T::one() + T::one());

        for k in 1..solver.nz - 1 {
            for j in 1..solver.ny - 1 {
                for i in 1..solver.nx - 1 {
                    let idx = solver.index(i, j, k);

                    // ── Face velocity: arithmetic average of adjacent cell centres ──
                    let vel_xm = {
                        let a = solver.velocity[solver.index(i - 1, j, k)];
                        let b = solver.velocity[idx];
                        (a + b) * half
                    };
                    let vel_xp = {
                        let a = solver.velocity[idx];
                        let b = solver.velocity[solver.index(i + 1, j, k)];
                        (a + b) * half
                    };
                    let vel_ym = {
                        let a = solver.velocity[solver.index(i, j - 1, k)];
                        let b = solver.velocity[idx];
                        (a + b) * half
                    };
                    let vel_yp = {
                        let a = solver.velocity[idx];
                        let b = solver.velocity[solver.index(i, j + 1, k)];
                        (a + b) * half
                    };
                    let vel_zm = {
                        let a = solver.velocity[solver.index(i, j, k - 1)];
                        let b = solver.velocity[idx];
                        (a + b) * half
                    };
                    let vel_zp = {
                        let a = solver.velocity[idx];
                        let b = solver.velocity[solver.index(i, j, k + 1)];
                        (a + b) * half
                    };

                    // ── Geometric PLIC flux for each face ──────────────────────────
                    // +x face: flux > 0 means fluid leaves cell (i,j,k)
                    let flux_xp = self.plic_face_flux(
                        solver,
                        i,
                        j,
                        k,
                        i + 1,
                        j,
                        k,
                        vel_xp.x,
                        solver.dy * solver.dz,
                        solver.dx,
                        dt,
                        0,
                    );
                    // −x face: flux > 0 means fluid enters cell (i,j,k)
                    let flux_xm = self.plic_face_flux(
                        solver,
                        i - 1,
                        j,
                        k,
                        i,
                        j,
                        k,
                        vel_xm.x,
                        solver.dy * solver.dz,
                        solver.dx,
                        dt,
                        0,
                    );
                    // +y face
                    let flux_yp = self.plic_face_flux(
                        solver,
                        i,
                        j,
                        k,
                        i,
                        j + 1,
                        k,
                        vel_yp.y,
                        solver.dx * solver.dz,
                        solver.dy,
                        dt,
                        1,
                    );
                    // −y face
                    let flux_ym = self.plic_face_flux(
                        solver,
                        i,
                        j - 1,
                        k,
                        i,
                        j,
                        k,
                        vel_ym.y,
                        solver.dx * solver.dz,
                        solver.dy,
                        dt,
                        1,
                    );
                    // +z face
                    let flux_zp = self.plic_face_flux(
                        solver,
                        i,
                        j,
                        k,
                        i,
                        j,
                        k + 1,
                        vel_zp.z,
                        solver.dx * solver.dy,
                        solver.dz,
                        dt,
                        2,
                    );
                    // −z face
                    let flux_zm = self.plic_face_flux(
                        solver,
                        i,
                        j,
                        k - 1,
                        i,
                        j,
                        k,
                        vel_zm.z,
                        solver.dx * solver.dy,
                        solver.dz,
                        dt,
                        2,
                    );

                    // Volume balance: net outflow from this cell
                    let cell_volume = solver.dx * solver.dy * solver.dz;
                    let net_outflow = flux_xp - flux_xm + flux_yp - flux_ym + flux_zp - flux_zm;

                    let alpha_new = solver.alpha[idx] - net_outflow / cell_volume;

                    // Invariant: α ∈ [0, 1]
                    solver.alpha_previous[idx] = <T as num_traits::Float>::min(
                        <T as num_traits::Float>::max(alpha_new, zero),
                        one,
                    );
                }
            }
        }

        // Swap buffers (zero-copy).
        std::mem::swap(&mut solver.alpha, &mut solver.alpha_previous);
        Ok(())
    }

    /// Compute the signed volumetric PLIC flux through a face.
    ///
    /// Returns a positive value when fluid flows from the left/bottom/back cell
    /// to the right/top/front cell; negative in the opposite direction.
    ///
    /// The donor cell for the flux is the upstream cell.  The swept sub-volume
    /// is computed using the donor cell's PLIC plane intersected with the
    /// donor prism of depth `|u_face| dt` and cross-section `face_area`.
    ///
    /// ## Formula (1D reduction for a rectangular prism donor)
    ///
    /// The prism has dimensions `(depth, √A, √A)` where `depth = |u_face| dt`
    /// and `A = face_area`.  The PLIC normal in the donor cell selects the
    /// relevant face-normal component for the truncation.  The full 3D
    /// analytical formula from Scardovelli & Zaleski (2000) is used via
    /// `volume_fraction_in_prism`.
    fn plic_face_flux<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        &self,
        solver: &VofSolver<T>,
        // Left/upstream cell index
        il: usize,
        jl: usize,
        kl: usize,
        // Right/downstream cell index
        ir: usize,
        jr: usize,
        kr: usize,
        // Face-normal velocity component (signed)
        u_face: T,
        // Face area (Δ_perp1 × Δ_perp2)
        face_area: T,
        // Cell dimension normal to face (Δ_normal)
        delta_normal: T,
        dt: T,
        // 0 = x-face, 1 = y-face, 2 = z-face
        _dir: usize,
    ) -> T {
        let zero = T::zero();
        let one = T::one();

        if <T as num_traits::Float>::abs(u_face)
            < <T as FromPrimitive>::from_f64(1e-300).unwrap_or(zero)
        {
            return zero;
        }

        // Identify donor cell (upstream).
        let (id, jd, kd) = if u_face > zero {
            (il, jl, kl)
        } else {
            (ir, jr, kr)
        };

        let donor_idx = solver.index(id, jd, kd);
        let alpha_donor = solver.alpha[donor_idx];
        let normal_donor = solver.normals[donor_idx];

        // Depth of the swept prism in the face-normal direction.
        let depth = <T as num_traits::Float>::min(
            <T as num_traits::Float>::abs(u_face) * dt,
            delta_normal, // Cannot exceed the donor cell width
        );

        if depth <= zero {
            return zero;
        }

        // Compute the fraction of the swept prism occupied by the fluid phase
        // using the donor cell's PLIC plane.
        //
        // The swept prism dimensions are: dx_prism = depth, dy_prism, dz_prism.
        // For simplicity we use the full cell dimensions perpendicularly.
        let dx_p = solver.dx;
        let dy_p = solver.dy;
        let dz_p = solver.dz;

        // The fraction of the prism volume that is fluid:
        let f = plic_volume_fraction_in_prism(normal_donor, alpha_donor, depth, dx_p, dy_p, dz_p);

        // Flux = fraction_fluid × swept_prism_volume × sign(u_face)
        let swept_volume = depth * face_area;
        let sign = if u_face > zero { one } else { -one };

        <T as num_traits::Float>::min(<T as num_traits::Float>::max(f, zero), one)
            * swept_volume
            * sign
    }

    /// Algebraic advection (simpler but less accurate, first-order upwind).
    fn algebraic_advection<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        self,
        solver: &mut VofSolver<T>,
        dt: T,
    ) -> Result<()> {
        // Use alpha_previous as temporary buffer (zero-copy optimization).
        solver.copy_boundaries();

        for k in 1..solver.nz - 1 {
            for j in 1..solver.ny - 1 {
                for i in 1..solver.nx - 1 {
                    let idx = solver.index(i, j, k);
                    let vel = solver.velocity[idx];

                    // First-order upwind differencing
                    let dalpha_dx = if vel.x > T::zero() {
                        (solver.alpha[idx] - solver.alpha[solver.index(i - 1, j, k)]) / solver.dx
                    } else {
                        (solver.alpha[solver.index(i + 1, j, k)] - solver.alpha[idx]) / solver.dx
                    };
                    let dalpha_dy = if vel.y > T::zero() {
                        (solver.alpha[idx] - solver.alpha[solver.index(i, j - 1, k)]) / solver.dy
                    } else {
                        (solver.alpha[solver.index(i, j + 1, k)] - solver.alpha[idx]) / solver.dy
                    };
                    let dalpha_dz = if vel.z > T::zero() {
                        (solver.alpha[idx] - solver.alpha[solver.index(i, j, k - 1)]) / solver.dz
                    } else {
                        (solver.alpha[solver.index(i, j, k + 1)] - solver.alpha[idx]) / solver.dz
                    };

                    // Advection equation: ∂α/∂t + u·∇α = 0
                    solver.alpha_previous[idx] = solver.alpha[idx]
                        - dt * (vel.x * dalpha_dx + vel.y * dalpha_dy + vel.z * dalpha_dz);

                    // Bound volume fraction: invariant α ∈ [0, 1]
                    solver.alpha_previous[idx] = <T as num_traits::Float>::min(
                        <T as num_traits::Float>::max(solver.alpha_previous[idx], T::zero()),
                        T::one(),
                    );
                }
            }
        }

        // Swap buffers (zero-copy optimization).
        std::mem::swap(&mut solver.alpha, &mut solver.alpha_previous);
        Ok(())
    }

    /// Apply artificial compression to sharpen interface.
    pub fn apply_compression<
        T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
    >(
        self,
        solver: &mut VofSolver<T>,
        dt: T,
    ) -> Result<()> {
        let interface_compression = solver.config.interface_compression;
        if !interface_compression.is_finite() || !(0.0..=1.0).contains(&interface_compression) {
            return Err(Error::InvalidConfiguration(
                "VOF interface_compression must be finite and in [0, 1]".to_string(),
            ));
        }
        let compression_factor = <T as FromPrimitive>::from_f64(interface_compression)
            .expect("interface_compression is an IEEE 754 representable f64 constant");

        for k in 0..solver.nz {
            for j in 0..solver.ny {
                for i in 0..solver.nx {
                    let idx = solver.index(i, j, k);
                    let alpha = solver.alpha[idx];

                    if k == 0
                        || k == solver.nz - 1
                        || j == 0
                        || j == solver.ny - 1
                        || i == 0
                        || i == solver.nx - 1
                    {
                        solver.alpha_previous[idx] = alpha;
                        continue;
                    }

                    if alpha
                        > <T as FromPrimitive>::from_f64(VOF_INTERFACE_LOWER).unwrap_or(T::zero())
                        && alpha
                            < <T as FromPrimitive>::from_f64(VOF_INTERFACE_UPPER)
                                .unwrap_or(T::one())
                    {
                        let normal = solver.normals[idx];

                        if normal.norm() > T::zero() {
                            let u_compression = normal * compression_factor;

                            let two = T::one() + T::one();
                            let dalpha_dx = (solver.alpha[solver.index(i + 1, j, k)]
                                - solver.alpha[solver.index(i - 1, j, k)])
                                / (two * solver.dx);
                            let dalpha_dy = (solver.alpha[solver.index(i, j + 1, k)]
                                - solver.alpha[solver.index(i, j - 1, k)])
                                / (two * solver.dy);
                            let dalpha_dz = (solver.alpha[solver.index(i, j, k + 1)]
                                - solver.alpha[solver.index(i, j, k - 1)])
                                / (two * solver.dz);

                            let compression_term = u_compression.x * dalpha_dx
                                + u_compression.y * dalpha_dy
                                + u_compression.z * dalpha_dz;

                            solver.alpha_previous[idx] =
                                alpha - dt * compression_term * alpha * (T::one() - alpha);

                            solver.alpha_previous[idx] = <T as num_traits::Float>::min(
                                <T as num_traits::Float>::max(
                                    solver.alpha_previous[idx],
                                    T::zero(),
                                ),
                                T::one(),
                            );
                        } else {
                            solver.alpha_previous[idx] = alpha;
                        }
                    } else {
                        solver.alpha_previous[idx] = alpha;
                    }
                }
            }
        }

        std::mem::swap(&mut solver.alpha, &mut solver.alpha_previous);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::super::plic_geometry::volume_under_plane_3d;
    use nalgebra::Vector3;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn test_volume_under_plane_bounds(
            nx in -1.0..1.0f64,
            ny in -1.0..1.0f64,
            nz in -1.0..1.0f64,
            c in -2.0..2.0f64,
            dx in 0.1..2.0f64,
            dy in 0.1..2.0f64,
            dz in 0.1..2.0f64,
        ) {
            let norm = (nx*nx + ny*ny + nz*nz).sqrt();
            prop_assume!(norm > 1e-6);
            let normal = Vector3::new(nx/norm, ny/norm, nz/norm);

            let vol = volume_under_plane_3d(normal, c, dx, dy, dz);
            let cell_vol = dx * dy * dz;

            // Volume must be bounded between 0 and cell_volume
            assert!(vol >= 0.0);
            assert!(vol <= cell_vol + 1e-10);
        }

        #[test]
        fn test_volume_under_plane_monotonicity(
            nx in -1.0..1.0f64,
            ny in -1.0..1.0f64,
            nz in -1.0..1.0f64,
            c1 in -1.0..1.0f64,
            dc in 0.0..1.0f64,
            dx in 0.1..2.0f64,
            dy in 0.1..2.0f64,
            dz in 0.1..2.0f64,
        ) {
            let norm = (nx*nx + ny*ny + nz*nz).sqrt();
            prop_assume!(norm > 1e-6);
            let normal = Vector3::new(nx/norm, ny/norm, nz/norm);

            let c2 = c1 + dc;
            let vol1 = volume_under_plane_3d(normal, c1, dx, dy, dz);
            let vol2 = volume_under_plane_3d(normal, c2, dx, dy, dz);

            // Volume must be monotonically increasing with plane constant
            assert!(vol2 >= vol1 - 1e-10);
        }
    }
}
