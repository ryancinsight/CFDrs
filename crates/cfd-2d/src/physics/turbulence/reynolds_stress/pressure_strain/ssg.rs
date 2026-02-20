//! Full SSG (Speziale-Sarkar-Gatski, 1991) pressure-strain model.
//!
//! ### Theorem: SSG Full Pressure-Strain
//! Φ_ij = − (C₁ ε + C₁* P) b_ij
//!         + C₂ ε (b_ik b_kj − 1/3 II_b δ_ij)
//!         + (C₃ − C₃* √II_b) k S_ij
//!         + C₄ k (b_ik S_kj + b_jk S_ki − 2/3 b_kl S_kl δ_ij)
//!         + C₅ k (b_ik W_kj − W_ik b_kj)
//!
//! **Reference**: Speziale, Sarkar & Gatski, J. Fluid Mech. 227:245-272 (1991).

use nalgebra::RealField;
use num_traits::FromPrimitive;

fn c<T: RealField + Copy + FromPrimitive>(v: f64) -> T {
    T::from_f64(v).expect("SSG constant must be representable")
}

/// Compute `Φ_ij` for the full SSG pressure-strain model.
#[allow(clippy::too_many_arguments)]
#[inline]
pub fn pressure_strain_ssg<T: RealField + Copy + FromPrimitive>(
    c1: T,
    c1_star: T,
    c2: T,
    c3: T,
    c3_star: T,
    c4: T,
    c5: T,
    a_xx: T,
    a_xy: T,
    a_yy: T,
    k: T,
    epsilon: T,
    s11: T,
    s12: T,
    s22: T,
    w12: T,
    w21: T,
    i: usize,
    j: usize,
) -> T {
    let a_zz = -(a_xx + a_yy);
    let ii_b = a_xx * a_xx + a_yy * a_yy + a_zz * a_zz + c::<T>(2.0) * a_xy * a_xy;
    let sqrt_ii_b = ii_b.sqrt();

    // P = -2k b:S
    let production = c::<T>(-2.0) * k * (a_xx * s11 + a_yy * s22 + c::<T>(2.0) * a_xy * s12);

    let c1_coeff = -(c1 * epsilon + c1_star * production);
    let c2_coeff = c2 * epsilon;
    let c3_coeff = (c3 - c3_star * sqrt_ii_b) * k;
    let c4_coeff = c4 * k;
    let c5_coeff = c5 * k;

    let one_third = c::<T>(1.0 / 3.0);
    let two_thirds = c::<T>(2.0 / 3.0);

    match (i, j) {
        (0, 0) => {
            let t1 = c1_coeff * a_xx;
            let b_sq_xx = a_xx * a_xx + a_xy * a_xy;
            let t2 = c2_coeff * (b_sq_xx - one_third * ii_b);
            let t3 = c3_coeff * s11;
            let m_xx = c::<T>(2.0) * (a_xx * s11 + a_xy * s12);
            let b_colon_s = a_xx * s11 + a_yy * s22 + c::<T>(2.0) * a_xy * s12;
            let t4 = c4_coeff * (m_xx - two_thirds * b_colon_s);
            let t5 = c5_coeff * c::<T>(2.0) * a_xy * w21;
            t1 + t2 + t3 + t4 + t5
        }
        (0, 1) | (1, 0) => {
            let t1 = c1_coeff * a_xy;
            let b_sq_xy = a_xx * a_xy + a_xy * a_yy;
            let t2 = c2_coeff * b_sq_xy;
            let t3 = c3_coeff * s12;
            let m_xy = a_xx * s12 + a_xy * s22 + a_xy * s11 + a_yy * s12;
            let b_colon_s = a_xx * s11 + a_yy * s22 + c::<T>(2.0) * a_xy * s12;
            let t4 = c4_coeff * (m_xy - two_thirds * b_colon_s * T::zero()); // δ_xy = 0
            let t5 = c5_coeff * (a_xx * w12 + a_yy * w21);
            t1 + t2 + t3 + t4 + t5
        }
        (1, 1) => {
            let t1 = c1_coeff * a_yy;
            let b_sq_yy = a_xy * a_xy + a_yy * a_yy;
            let t2 = c2_coeff * (b_sq_yy - one_third * ii_b);
            let t3 = c3_coeff * s22;
            let m_yy = c::<T>(2.0) * (a_xy * s12 + a_yy * s22);
            let b_colon_s = a_xx * s11 + a_yy * s22 + c::<T>(2.0) * a_xy * s12;
            let t4 = c4_coeff * (m_yy - two_thirds * b_colon_s);
            let t5 = c5_coeff * c::<T>(2.0) * a_xy * w12;
            t1 + t2 + t3 + t4 + t5
        }
        _ => T::zero(),
    }
}
