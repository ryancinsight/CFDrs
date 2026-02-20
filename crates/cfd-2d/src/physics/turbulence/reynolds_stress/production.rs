//! Reynolds stress production term.
//!
//! ### Theorem: Exact Production (Pope, 2000; Launder et al., 1975)
//!
//! ```text
//! P_ij = −⟨u_i'u_k'⟩ ∂U_j/∂x_k − ⟨u_j'u_k'⟩ ∂U_i/∂x_k
//! ```
//!
//! This is **not** the Boussinesq approximation `P_ij = −2 ν_t S_ij`.
//! The full tensor contraction captures anisotropy generation within the
//! second-moment closure (Launder, Reece & Rodi, 1975).

use super::tensor::ReynoldsStressTensor;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Compute the production tensor component `P_ij` at grid point `(x, y)`.
///
/// # Arguments
/// * `rs` — current Reynolds stress tensor field
/// * `velocity_gradient` — `[[dU/dx, dU/dy], [dV/dx, dV/dy]]` at `(x, y)`
/// * `i, j` — tensor index pair (symmetric: `(i,j) == (j,i)`)
/// * `x, y` — grid point indices
#[inline]
pub fn production_term<T: RealField + Copy + FromPrimitive>(
    rs: &ReynoldsStressTensor<T>,
    velocity_gradient: &[[T; 2]; 2],
    i: usize,
    j: usize,
    x: usize,
    y: usize,
) -> T {
    let du_dx = velocity_gradient[0][0];
    let du_dy = velocity_gradient[0][1];
    let dv_dx = velocity_gradient[1][0];
    let dv_dy = velocity_gradient[1][1];

    let xx = rs.xx[(x, y)];
    let xy = rs.xy[(x, y)];
    let yy = rs.yy[(x, y)];

    let two = T::from_f64(2.0).expect("2.0 must be representable");

    match (i, j) {
        (0, 0) => -two * xx * du_dx - two * xy * du_dy,
        (0, 1) | (1, 0) => -xx * dv_dx - xy * dv_dy - xy * du_dx - yy * du_dy,
        (1, 1) => -two * xy * dv_dx - two * yy * dv_dy,
        _ => T::zero(),
    }
}
