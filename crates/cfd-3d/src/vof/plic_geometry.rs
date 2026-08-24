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

use crate::scalar;
use cfd_core::CfdScalar;
use eunomia::FloatElement;
use leto::geometry::Vector3;

/// Compute the fluid volume fraction inside a rectangular prism cut by a PLIC plane.
///
/// ## Theorem (PLIC Prism Volume — Scardovelli & Zaleski 2000)
///
/// Given a unit-normal vector **n** and a plane constant `C` (i.e., the plane
/// **n**·**x** = `C`), the volume `V` of a rectangular cell
/// `[0,Δx]×[0,Δy]×[0,Δz]` on the side Σ|nᵢ|xᵢ ≤ `C` is given by the analytical
/// 5-region piecewise polynomial of Eq. (2.34)–(2.38) of Scardovelli & Zaleski
/// (2000), implemented in [`volume_under_plane_3d`].
///
/// The donor-flux evaluation in [`plic_volume_fraction_in_prism`] treats the
/// gas as the corner cut (the interface normal points into the fluid) and
/// intersects that cut with the swept slab adjacent to the outflow face.

/// Fluid-volume fraction of the swept PLIC prism in a donor cell.
///
/// # Convention
///
/// The interface normal points **into the fluid** (negative α-gradient), so
/// the gas occupies the origin corner cut `{x : Σ|nᵢ|xᵢ ≤ C}` of the donor
/// cell. `C` is bisected so this cut holds exactly `(1 − α)` of the cell.
///
/// The swept prism is the slab of thickness `depth` adjacent to the donor's
/// **outflow** face along `axis`: `[Δₐ − depth, Δₐ]` when `flow_sign > 0`
/// (outflow through the high-coordinate face) and `[0, depth]` otherwise.
///
/// Prior to 2026-08-20 this function matched the corner cut to `α` itself
/// (mirroring fluid and gas) and evaluated the slab at the donor-cell origin
/// regardless of flow direction, so partially filled donor cells donated
/// wrong fluxes in both directions.
pub(crate) fn plic_volume_fraction_in_prism<T: CfdScalar>(
    normal: Vector3<T>,
    alpha_donor: T,
    depth: T,
    dx: T,
    dy: T,
    dz: T,
    axis: usize,
    flow_sign: T,
) -> T {
    let zero = scalar::zero::<T>();
    let one = scalar::one::<T>();

    // If the cell is entirely full or empty, the swept prism carries the same fraction.
    if alpha_donor <= zero {
        return zero;
    }
    if alpha_donor >= one {
        return one;
    }

    let prism_vol = depth * dy * dz;
    if prism_vol <= zero {
        return zero;
    }

    // Full-cell scaled absolute components and the axis-specific values.
    let (delta_a, n_a_signed) = match axis {
        0 => (dx, normal.x),
        1 => (dy, normal.y),
        _ => (dz, normal.z),
    };
    let m_full = [
        scalar::abs(normal.x) * dx,
        scalar::abs(normal.y) * dy,
        scalar::abs(normal.z) * dz,
    ];
    let n_a = scalar::abs(n_a_signed);

    // The interface normal points into the fluid, so the GAS occupies the
    // origin corner cut {Σ|nᵢ|xᵢ ≤ C} of the donor cell. Bisect C so the cut
    // holds exactly (1 − α) of the cell. Prior to 2026-08-20 the bisection
    // targeted α itself, mirroring fluid and gas.
    let gas_target = one - alpha_donor;
    let mut c_lo = zero;
    let mut c_hi = m_full[0] + m_full[1] + m_full[2];
    for _ in 0..64 {
        let c_mid = c_lo + (c_hi - c_lo) * <T as FloatElement>::from_f64(0.5);
        if corner_cut_volume_fraction(m_full, c_mid) < gas_target {
            c_lo = c_mid;
        } else {
            c_hi = c_mid;
        }
    }
    let c_gas = c_lo + (c_hi - c_lo) * <T as FloatElement>::from_f64(0.5);

    // Distance from the donor-cell origin to the near edge of the swept slab
    // along the axis. Outflow through the high-coordinate face (flow_sign > 0)
    // sweeps [Δₐ − depth, Δₐ]; otherwise [0, depth].
    let slab_offset = if flow_sign > zero {
        delta_a - depth
    } else {
        zero
    };

    // Convert the corner constant into the signed plane position d* (with
    // n·x ≤ d* describing the gas), shift by the slab offset using the SIGNED
    // axis component, then convert back with the slab's own mirror offset:
    //
    //   d*          = C_gas + Σ_{nᵢ<0} nᵢΔᵢ
    //   C_slab      = d* − n_a·t − Σ_{nᵢ<0} nᵢΔᵢ^slab
    //               = C_gas + min(n_a, 0)·(Δ_a − depth) − n_a·slab_offset
    //
    // with n_a the SIGNED axis component throughout. Evaluating this
    // expression for a mirrored normal flips which slab is fluid — the
    // orientation-blind absolute-value form cannot.
    let negative_axis = scalar::min(n_a_signed, zero);
    let c_gas_slab = c_gas + negative_axis * (delta_a - depth) - n_a_signed * slab_offset;

    // Gas fraction inside the slab prism (extent `depth` along the axis, full
    // cell dimensions perpendicular).
    let m_slab = match axis {
        0 => [n_a * depth, m_full[1], m_full[2]],
        1 => [m_full[0], n_a * depth, m_full[2]],
        _ => [m_full[0], m_full[1], n_a * depth],
    };
    let gas_frac = corner_cut_volume_fraction(m_slab, c_gas_slab);

    let fluid_frac = one - gas_frac;
    scalar::min(scalar::max(fluid_frac, zero), one)
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
pub fn volume_under_plane_3d<T: CfdScalar>(
    normal: Vector3<T>,
    plane_constant: T,
    dx: T,
    dy: T,
    dz: T,
) -> T {
    let cell_volume = dx * dy * dz;

    // Absolute normal scaled by cell dimensions.
    let m = [
        scalar::abs(normal.x) * dx,
        scalar::abs(normal.y) * dy,
        scalar::abs(normal.z) * dz,
    ];

    corner_cut_volume_fraction(m, plane_constant) * cell_volume
}

/// Fraction (in `[0, 1]`) of a rectangular box cut off at its origin corner
/// by `{x : m₀x₀ + m₁x₁ + m₂x₂ ≤ C}`, where `mᵢ = |nᵢ|Δᵢ ≥ 0`.
///
/// Uses the inclusion–exclusion formula (Pilliod & Puckett 2004)
///
/// ```text
/// F(C) = [C₊³ − Σᵢ(C−mᵢ)₊³ + Σᵢ<ⱼ(C−mᵢ−mⱼ)₊³ − (C−Σ)₊³] / (6·m₀·m₁·m₂)
/// ```
///
/// reduced to the 2-D and 1-D limit formulas whenever one or two `mᵢ`
/// vanish. Prior to 2026-08-20 the caller divided by `6·m₀·m₁·m₂`
/// unconditionally, producing infinities for axis-aligned normals.
///
/// Monotonicity follows because dF/dC is the cross-sectional area of the
/// plane–box intersection, which is ≥ 0.
fn corner_cut_volume_fraction<T: CfdScalar>(m: [T; 3], c: T) -> T {
    let zero = scalar::zero::<T>();
    let eps = <T as FloatElement>::from_f64(1e-14);

    let mut active: [T; 3] = [zero; 3];
    let mut k = 0_usize;
    for &component in &m {
        if component > eps {
            active[k] = component;
            k += 1;
        }
    }
    active[..k].sort_by(|lhs, rhs| lhs.partial_cmp(rhs).unwrap_or(std::cmp::Ordering::Equal));

    let pos_pow = |x: T, exp: u32| -> T {
        if x <= zero {
            zero
        } else {
            FloatElement::powi(x, exp as i32)
        }
    };

    match k {
        0 => <T as FloatElement>::from_f64(0.5),
        1 => scalar::min(scalar::max(c / active[0], zero), scalar::one::<T>()),
        2 => {
            let (m0, m1) = (active[0], active[1]);
            let two = <T as FloatElement>::from_f64(2.0);
            let numerator =
                pos_pow(c, 2) - pos_pow(c - m0, 2) - pos_pow(c - m1, 2) + pos_pow(c - m0 - m1, 2);
            scalar::min(
                scalar::max(numerator / (two * m0 * m1), zero),
                scalar::one::<T>(),
            )
        }
        _ => {
            let (m0, m1, m2) = (active[0], active[1], active[2]);
            let six = <T as FloatElement>::from_f64(6.0);
            let numerator =
                pos_pow(c, 3) - pos_pow(c - m0, 3) - pos_pow(c - m1, 3) - pos_pow(c - m2, 3)
                    + pos_pow(c - m0 - m1, 3)
                    + pos_pow(c - m0 - m2, 3)
                    + pos_pow(c - m1 - m2, 3)
                    - pos_pow(c - m0 - m1 - m2, 3);
            scalar::min(
                scalar::max(numerator / (six * m0 * m1 * m2), zero),
                scalar::one::<T>(),
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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

    #[test]
    fn swept_slab_follows_flow_direction() {
        // Half-full cell, normal (1, 0, 0) points into the fluid ⇒ fluid at
        // high x. Outflow through the +x face sweeps an all-fluid slab;
        // outflow through the −x face sweeps an all-gas slab.
        let normal = Vector3::new(1.0, 0.0, 0.0);
        let alpha = 0.5_f64;
        let (dx, dy, dz) = (1.0_f64, 1.0_f64, 1.0_f64);
        let depth = 0.2_f64;

        let f_out_high: f64 =
            plic_volume_fraction_in_prism(normal, alpha, depth, dx, dy, dz, 0, 1.0);
        let f_out_low: f64 =
            plic_volume_fraction_in_prism(normal, alpha, depth, dx, dy, dz, 0, -1.0);

        assert!(
            (f_out_high - 1.0).abs() < 1e-9,
            "high-side outflow f={f_out_high}"
        );
        assert!(
            (f_out_low - 0.0).abs() < 1e-9,
            "low-side outflow f={f_out_low}"
        );

        // Mirrored normal puts the fluid at low x; fractions swap.
        let mirrored = Vector3::new(-1.0, 0.0, 0.0);
        let f_mirrored_high: f64 =
            plic_volume_fraction_in_prism(mirrored, alpha, depth, dx, dy, dz, 0, 1.0);
        let f_mirrored_low: f64 =
            plic_volume_fraction_in_prism(mirrored, alpha, depth, dx, dy, dz, 0, -1.0);

        assert!((f_mirrored_high - 0.0).abs() < 1e-9);
        assert!((f_mirrored_low - 1.0).abs() < 1e-9);
    }

    #[test]
    fn swept_slab_partial_coverage_is_exact() {
        // α = 0.75 with normal (1, 0, 0): gas fills x ≤ 0.25. A slab of
        // depth 0.5 on the low side ([0, 0.5]) is half gas, half fluid.
        let normal = Vector3::new(1.0, 0.0, 0.0);
        let f_low: f64 = plic_volume_fraction_in_prism(normal, 0.75, 0.5, 1.0, 1.0, 1.0, 0, -1.0);
        assert!((f_low - 0.5).abs() < 1e-9, "partial slab f={f_low}");

        // The high-side slab [0.5, 1.0] is entirely fluid.
        let f_high: f64 = plic_volume_fraction_in_prism(normal, 0.75, 0.5, 1.0, 1.0, 1.0, 0, 1.0);
        assert!((f_high - 1.0).abs() < 1e-9);
    }

    #[test]
    fn swept_slab_handles_all_axes_and_degenerate_fill() {
        // y-axis analogue of the directional test.
        let normal_y = Vector3::new(0.0, 1.0, 0.0);
        let f_y_high: f64 =
            plic_volume_fraction_in_prism(normal_y, 0.5, 0.2, 1.0, 1.0, 1.0, 1, 1.0);
        let f_y_low: f64 =
            plic_volume_fraction_in_prism(normal_y, 0.5, 0.2, 1.0, 1.0, 1.0, 1, -1.0);
        assert!((f_y_high - 1.0).abs() < 1e-9);
        assert!((f_y_low - 0.0).abs() < 1e-9);

        // z-axis analogue.
        let normal_z = Vector3::new(0.0, 0.0, 1.0);
        let f_z_high: f64 =
            plic_volume_fraction_in_prism(normal_z, 0.5, 0.2, 1.0, 1.0, 1.0, 2, 1.0);
        let f_z_low: f64 =
            plic_volume_fraction_in_prism(normal_z, 0.5, 0.2, 1.0, 1.0, 1.0, 2, -1.0);
        assert!((f_z_high - 1.0).abs() < 1e-9);
        assert!((f_z_low - 0.0).abs() < 1e-9);

        // Full/empty cells carry their fraction regardless of direction.
        let full: f64 = plic_volume_fraction_in_prism(normal_x(), 1.0, 0.2, 1.0, 1.0, 1.0, 0, -1.0);
        let empty: f64 = plic_volume_fraction_in_prism(normal_x(), 0.0, 0.2, 1.0, 1.0, 1.0, 0, 1.0);
        assert!((full - 1.0).abs() < 1e-12);
        assert!((empty - 0.0).abs() < 1e-12);
    }

    #[test]
    fn axis_aligned_normal_volume_is_exact() {
        // Axis-aligned normals previously divided by 6·m₀·m₁·m₂ = 0 and
        // relied on NaN-propagating clamps. The reduced 1-D/2-D corner-cut
        // formulas must return exact slab fractions instead.
        let normal = Vector3::new(1.0, 0.0, 0.0);
        let vol: f64 = volume_under_plane_3d(normal, 0.25, 1.0, 1.0, 1.0);
        assert!((vol - 0.25).abs() < 1e-12, "slab volume {vol}");

        let full: f64 = volume_under_plane_3d(normal, 2.0, 1.0, 1.0, 1.0);
        let empty: f64 = volume_under_plane_3d(normal, -0.1, 1.0, 1.0, 1.0);
        assert!((full - 1.0).abs() < 1e-12);
        assert!(empty.abs() < 1e-12);

        // Two active components exercise the 2-D limit formula:
        // {x + y ≤ 0.5} in the unit square has area 0.125.
        let diagonal = Vector3::new(1.0, 1.0, 0.0);
        let area: f64 = volume_under_plane_3d(diagonal, 0.5, 1.0, 1.0, 1.0);
        assert!((area - 0.125).abs() < 1e-12, "wedge area {area}");
    }

    fn normal_x() -> Vector3<f64> {
        Vector3::new(1.0, 0.0, 0.0)
    }
}
