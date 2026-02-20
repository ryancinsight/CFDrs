//! Linear Return-to-Isotropy pressure-strain model (Rotta, 1951).
//!
//! ### Theorem: Linear Relaxation
//! Φ_ij = −C₁ (ε/k) (⟨u_i'u_j'⟩ − (2/3)k δ_ij)
//!
//! **Proof**: The anisotropy tensor b_ij = ⟨u_i'u_j'⟩/k − δ_ij/3 decays
//! exponentially with time-scale k/ε under the linear model, recovering the
//! isotropic state as the unique fixed point (Rotta, 1951).

use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Compute `Φ_ij` for the linear return-to-isotropy (Rotta, 1951) model.
///
/// # Arguments
/// * `c1` — Rotta constant (≈ 1.8)
/// * `a_ij` — anisotropy tensor component for indices (i,j)
/// * `epsilon` — dissipation rate ε
/// * `k`       — turbulent kinetic energy
/// * `i`, `j`  — tensor index pair
#[inline]
pub fn pressure_strain_linear<T: RealField + Copy + FromPrimitive>(
    c1: T,
    a_xx: T,
    a_xy: T,
    a_yy: T,
    epsilon: T,
    k: T,
    i: usize,
    j: usize,
) -> T {
    let c1_term = -c1 * epsilon / k;
    match (i, j) {
        (0, 0) => c1_term * a_xx,
        (0, 1) | (1, 0) => c1_term * a_xy,
        (1, 1) => c1_term * a_yy,
        _ => T::zero(),
    }
}
