//! PLIC (Piecewise Linear Interface Calculation) geometry primitives.
//!
//! These functions implement the analytical volume computation for a rectangular
//! cell cut by a plane **n**·**x** = C, following the inclusion-exclusion formula
//! of Scardovelli & Zaleski (2000) and the bisection-based plane constant
//! inversion of Pilliod & Puckett (2004).
//!
//! # References
//!
//! - Scardovelli, R. & Zaleski, S. (2000). "Analytical relations connecting
//!   linear interfaces and volume fractions in rectangular grids".
//!   J. Comput. Phys. 164:228–237.
//! - Pilliod, J.E. & Puckett, E.G. (2004). "Second-order accurate
//!   volume-of-fluid algorithms for tracking material interfaces".
//!   J. Comput. Phys. 199:465–502.

use nalgebra::RealField;
use num_traits::FromPrimitive;

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
pub(crate) fn plic_volume_fraction_in_prism<
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
>(
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
pub(crate) fn find_plic_plane_constant<
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
>(
    normal: nalgebra::Vector3<T>,
    alpha: T,
    dx: T,
    dy: T,
    dz: T,
) -> T {
    use num_traits::Float;
    let zero = T::zero();
    let half = T::one() / (T::one() + T::one());

    let n_abs = Float::abs(normal.x) * dx + Float::abs(normal.y) * dy + Float::abs(normal.z) * dz;
    let mut c_lo = zero;
    let mut c_hi = n_abs;
    let target = alpha * dx * dy * dz;
    let tol = <T as FromPrimitive>::from_f64(1e-12)
        .expect("1e-12 is an IEEE 754 representable f64 constant")
        * dx
        * dy
        * dz;

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
pub fn volume_under_plane_3d<
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy,
>(
    normal: nalgebra::Vector3<T>,
    plane_constant: T,
    dx: T,
    dy: T,
    dz: T,
) -> T {
    use num_traits::Float;

    let zero = T::zero();
    let cell_volume = dx * dy * dz;
    let six = T::one() + T::one() + T::one() + T::one() + T::one() + T::one();

    // Absolute normal scaled by cell dimensions.
    let m1 = Float::abs(normal.x) * dx;
    let m2 = Float::abs(normal.y) * dy;
    let m3 = Float::abs(normal.z) * dz;
    let m_sum = m1 + m2 + m3;

    // Degenerate case: all normal components essentially zero.
    let eps = <T as FromPrimitive>::from_f64(1e-14).unwrap_or(zero);
    if m1 + m2 + m3 < eps {
        return (T::one() / (T::one() + T::one())) * cell_volume;
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
        if x <= zero {
            zero
        } else {
            x * x * x
        }
    };

    let denom = six * m1 * m2 * m3;

    let numerator = cube_pos(c) - cube_pos(c - m1) - cube_pos(c - m2) - cube_pos(c - m3)
        + cube_pos(c - m1 - m2)
        + cube_pos(c - m1 - m3)
        + cube_pos(c - m2 - m3)
        - cube_pos(c - m_sum);

    let volume_fraction = numerator / denom;

    // Clamp to [0, cell_volume] for numerical safety.
    Float::min(Float::max(volume_fraction * cell_volume, zero), cell_volume)
}

#[cfg(test)]
mod tests {
    use super::*;
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
