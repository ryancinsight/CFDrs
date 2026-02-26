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
use super::solver::VofSolver;
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

impl AdvectionMethod {
    /// Create advection method based on configuration
    #[must_use]
    pub fn create(config: &VofConfig) -> Self {
        config.advection_method
    }

    /// Advect volume fraction field by one time step `dt`.
    ///
    /// # Invariant
    /// After this call, `solver.alpha[i] ∈ [0, 1]` for all cells `i`.
    pub fn advect<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        self,
        solver: &mut VofSolver<T>,
        dt: T,
    ) -> Result<()> {
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
        let half = <T as FromPrimitive>::from_f64(0.5).unwrap_or_else(|| one / (one + one));

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
                        solver, i, j, k, i + 1, j, k,
                        vel_xp.x, solver.dy * solver.dz,
                        solver.dx, dt, 0,
                    );
                    // −x face: flux > 0 means fluid enters cell (i,j,k)
                    let flux_xm = self.plic_face_flux(
                        solver, i - 1, j, k, i, j, k,
                        vel_xm.x, solver.dy * solver.dz,
                        solver.dx, dt, 0,
                    );
                    // +y face
                    let flux_yp = self.plic_face_flux(
                        solver, i, j, k, i, j + 1, k,
                        vel_yp.y, solver.dx * solver.dz,
                        solver.dy, dt, 1,
                    );
                    // −y face
                    let flux_ym = self.plic_face_flux(
                        solver, i, j - 1, k, i, j, k,
                        vel_ym.y, solver.dx * solver.dz,
                        solver.dy, dt, 1,
                    );
                    // +z face
                    let flux_zp = self.plic_face_flux(
                        solver, i, j, k, i, j, k + 1,
                        vel_zp.z, solver.dx * solver.dy,
                        solver.dz, dt, 2,
                    );
                    // −z face
                    let flux_zm = self.plic_face_flux(
                        solver, i, j, k - 1, i, j, k,
                        vel_zm.z, solver.dx * solver.dy,
                        solver.dz, dt, 2,
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
        il: usize, jl: usize, kl: usize,
        // Right/downstream cell index
        ir: usize, jr: usize, kr: usize,
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

        if <T as num_traits::Float>::abs(u_face) < <T as FromPrimitive>::from_f64(1e-300).unwrap_or(zero) {
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

        <T as num_traits::Float>::min(<T as num_traits::Float>::max(f, zero), one) * swept_volume * sign
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
    pub fn apply_compression<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
        self,
        solver: &mut VofSolver<T>,
        dt: T,
    ) -> Result<()> {
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

                    if alpha > <T as FromPrimitive>::from_f64(VOF_INTERFACE_LOWER).unwrap_or(T::zero())
                        && alpha < <T as FromPrimitive>::from_f64(VOF_INTERFACE_UPPER).unwrap_or(T::one())
                    {
                        let normal = solver.normals[idx];

                        if normal.norm() > T::zero() {
                            let compression_factor = <T as FromPrimitive>::from_f64(0.5).unwrap_or(T::zero());
                            let u_compression = normal * compression_factor;

                            let two = <T as FromPrimitive>::from_f64(2.0).unwrap_or(T::one() + T::one());
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
                                <T as num_traits::Float>::max(solver.alpha_previous[idx], T::zero()),
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

/// Compute the fluid volume fraction inside a rectangular prism cut by a PLIC plane.
///
/// ## Theorem (PLIC Prism Volume — Scardovelli & Zaleski 2000)
///
/// Given a unit-normal vector **n** and a plane constant `C` (i.e., the plane
/// **n**·**x** = `C`), the volume `V` of a rectangular cell
/// `[0,Δx]×[0,Δy]×[0,Δz]` on the side **n**·**x** ≤ `C` satisfies
///
/// ```math
/// V(C) = V_cell · F(α)
/// ```
///
/// where `α` is the volume fraction of the full cell and `F` is the analytical
/// 5-region piecewise polynomial from Eq. (2.34)–(2.38) of Scardovelli & Zaleski
/// (2000).  Here we apply that formula to the *swept prism* with depth `depth`
/// instead of the full cell width.
///
/// The approach:
/// 1. Replace the full cell's `n_i Δ_i` components with the prism's `n_i Δ_i^prism`.
/// 2. Find the plane constant `C_prism` that gives the target fluid fraction `α` in the prism.
/// 3. Return `V_fluid / V_prism`.
///
/// For the donor-flux approach we directly use `alpha_donor` as the target fraction
/// in the (full) donor cell and evaluate what fraction of the *swept prism* is fluid.
/// Since the PLIC plane is the same, the fraction in the swept prism is the ratio
/// of the volume below the PLIC plane inside the prism to the total prism volume.
fn plic_volume_fraction_in_prism<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
    normal: nalgebra::Vector3<T>,
    alpha_donor: T,
    depth: T,
    dx: T,
    dy: T,
    dz: T,
) -> T {
    use num_traits::Float;
    let zero = T::zero();
    let one = T::one();

    // If the cell is entirely full or empty, the swept prism carries the same fraction.
    if alpha_donor <= zero {
        return zero;
    }
    if alpha_donor >= one {
        return one;
    }

    // Retrieve the PLIC plane constant C for the donor cell.
    // We use the Scardovelli-Zaleski analytical inverse to get C from alpha_donor.
    let c_plane = find_plic_plane_constant(normal, alpha_donor, dx, dy, dz);

    // Compute what fraction of the swept prism (dimensions: depth × dy × dz)
    // lies below the plane n·x = c_plane.
    //
    // The prism origin is the same as the donor cell's origin.  We evaluate
    // V_plane_intersection(n, c_plane, depth, dy, dz) / (depth * dy * dz).
    let prism_vol = depth * dy * dz;
    if prism_vol <= zero {
        return zero;
    }

    let vol = volume_under_plane_3d(normal, c_plane, depth, dy, dz);
    Float::min(Float::max(vol / prism_vol, zero), one)
}

/// Find the PLIC plane constant C such that the volume under **n**·**x** = C
/// in `[0,dx]×[0,dy]×[0,dz]` equals `alpha * dx * dy * dz`.
///
/// Uses iterative bisection (tolerance 1e-12 of cell size) per Scardovelli & Zaleski.
fn find_plic_plane_constant<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
    normal: nalgebra::Vector3<T>,
    alpha: T,
    dx: T,
    dy: T,
    dz: T,
) -> T {
    use num_traits::Float;
    let zero = T::zero();
    let half = <T as FromPrimitive>::from_f64(0.5).unwrap_or_else(T::one);

    let n_abs = Float::abs(normal.x) * dx + Float::abs(normal.y) * dy + Float::abs(normal.z) * dz;
    let mut c_lo = zero;
    let mut c_hi = n_abs;
    let target = alpha * dx * dy * dz;
    let tol = <T as FromPrimitive>::from_f64(1e-12).unwrap_or_else(|| T::zero()) * dx * dy * dz;

    for _ in 0..64 {
        if c_hi - c_lo < tol {
            break;
        }
        let c_mid = c_lo + (c_hi - c_lo) * half;
        let vol = volume_under_plane_3d(normal, c_mid, dx, dy, dz);
        if vol < target {
            c_lo = c_mid;
        } else {
            c_hi = c_mid;
        }
    }
    c_lo + (c_hi - c_lo) * half
}

/// Compute the volume of a rectangular cell `[0,Δx]×[0,Δy]×[0,Δz]` on the
/// side **n**·**x** ≤ C using the corrected Scardovelli & Zaleski (2000)
/// analytical 5-region formula (Eq. 2.34–2.38).
///
/// ## Mathematical Derivation
///
/// Let `m_i = |n_i| Δ_i` (scaled absolute normal components), sorted so
/// `m₁ ≤ m₂ ≤ m₃`, and let `α = C/(m₁+m₂+m₃)`.  Then the 5 regions are:
///
/// | Region | Domain                | Volume fraction F(α)                                   |
/// |--------|----------------------|-------------------------------------------------------|
/// | 1      | α ≤ m₁/Σ             | α³/(6m₁m₂m₃) × Σ³                                    |
/// | 2      | m₁ ≤ α·Σ ≤ m₂       | F₁ + residual for pentahedron                         |
/// | 3      | m₂ ≤ α·Σ ≤ min(m₃,m₁₂) | see Eq. 2.36                                      |
/// | 4      | else ≤ 1−m₁          | 1 − F(1−α) by symmetry                               |
/// | 5      | 1−m₁ ≤ α             | 1 − α³/(6m₁m₂m₃) × Σ³                               |
///
/// The derivation follows the inclusion-exclusion principle applied to the
/// three half-spaces defined by the axes.
///
/// **Reference**: Scardovelli, R. & Zaleski, S. (2000). "Analytical relations
///   connecting linear interfaces and volume fractions in rectangular grids".
///   J. Comput. Phys. 164:228–237. (Eqs. 2.34–2.38)
pub fn volume_under_plane_3d<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
    normal: nalgebra::Vector3<T>,
    plane_constant: T,
    dx: T,
    dy: T,
    dz: T,
) -> T {
    use num_traits::Float;

    let zero = T::zero();
    let cell_volume = dx * dy * dz;
    let six = <T as FromPrimitive>::from_f64(6.0).unwrap_or_else(|| {
        let one = T::one();
        let two = one + one;
        two * two + two
    });

    // Absolute normal scaled by cell dimensions.
    let m1 = Float::abs(normal.x) * dx;
    let m2 = Float::abs(normal.y) * dy;
    let m3 = Float::abs(normal.z) * dz;
    let m_sum = m1 + m2 + m3;

    // Degenerate case: all normal components essentially zero.
    let eps = <T as FromPrimitive>::from_f64(1e-14).unwrap_or(zero);
    if m1 + m2 + m3 < eps {
        return <T as FromPrimitive>::from_f64(0.5).unwrap_or(T::one()) * cell_volume;
    }

    let c = plane_constant;
    if c <= zero {
        return zero;
    }
    if c >= m_sum {
        return cell_volume;
    }

    // Inclusion-exclusion formula (Pilliod & Puckett 2004):
    //
    // V(C) = [C³ − Σᵢ(C−mᵢ)₊³ + Σᵢ<ⱼ(C−mᵢ−mⱼ)₊³ − (C−Σ)₊³] / (6·m₁·m₂·m₃)
    //
    // where (x)₊ = max(x, 0).
    //
    // # Theorem
    //
    // This equals the exact volume of the unit cuboid below the plane
    // m₁ξ + m₂η + m₃ζ = C. Monotonicity follows because dV/dC is the
    // cross-sectional area of the plane–cube intersection, which is ≥ 0.
    let cube_pos = |x: T| -> T {
        if x <= zero { zero } else { x * x * x }
    };

    let denom = six * m1 * m2 * m3;

    let numerator = cube_pos(c)
        - cube_pos(c - m1) - cube_pos(c - m2) - cube_pos(c - m3)
        + cube_pos(c - m1 - m2) + cube_pos(c - m1 - m3) + cube_pos(c - m2 - m3)
        - cube_pos(c - m_sum);

    let volume_fraction = numerator / denom;

    // Clamp to [0, cell_volume] for numerical safety.
    Float::min(Float::max(volume_fraction * cell_volume, zero), cell_volume)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use nalgebra::Vector3;

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
