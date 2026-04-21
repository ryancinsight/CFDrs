//! Healthy-cell protection composite for SDT report metrics.
//!
//! # Theorem - Bounded geometric protection composite
//!
//! If `wbc_targeted_cavitation` and `rbc_venturi_protection` are both in
//! `[0, 1]`, then `healthy_cell_protection_index` is also in `[0, 1]`.
//!
//! **Proof sketch**
//! `1 - wbc_targeted_cavitation` is the WBC sparing fraction, so it lies in
//! `[0, 1]` after clamping. The product of two numbers in `[0, 1]` remains in
//! `[0, 1]`, and the principal square root maps `[0, 1]` to `[0, 1]`.
//!
//! # Theorem - Monotone protection response
//!
//! For fixed RBC protection, the composite decreases as WBC-targeted
//! cavitation increases. For fixed WBC targeting, the composite increases as
//! RBC venturi protection increases.
//!
//! **Proof sketch**
//! The composite is the geometric mean of WBC sparing and RBC protection.
//! Increasing WBC-targeted cavitation reduces WBC sparing, while increasing
//! RBC venturi protection increases the product. The square root is monotone
//! over non-negative reals, so the ordering is preserved.

/// Conservative healthy-cell protection composite.
///
/// The index is the geometric mean of WBC sparing and RBC venturi protection:
/// `sqrt((1 - wbc_targeted_cavitation) × rbc_venturi_protection)`.
#[must_use]
pub fn healthy_cell_protection_index(
    wbc_targeted_cavitation: f64,
    rbc_venturi_protection: f64,
) -> f64 {
    let wbc_sparing = (1.0 - wbc_targeted_cavitation).clamp(0.0, 1.0);
    let rbc_protection = rbc_venturi_protection.clamp(0.0, 1.0);
    (wbc_sparing * rbc_protection).sqrt()
}

#[cfg(test)]
mod tests {
    use super::healthy_cell_protection_index;
    use proptest::prelude::*;

    #[test]
    fn healthy_cell_protection_index_matches_boundary_cases() {
        assert_eq!(healthy_cell_protection_index(1.0, 0.8), 0.0);
        assert_eq!(healthy_cell_protection_index(0.2, 0.0), 0.0);
        assert_eq!(healthy_cell_protection_index(0.0, 1.0), 1.0);
    }

    proptest! {
        #[test]
        fn healthy_cell_protection_index_stays_bounded(
            wbc_targeted_cavitation in 0.0f64..1.0,
            rbc_venturi_protection in 0.0f64..1.0,
        ) {
            let index = healthy_cell_protection_index(
                wbc_targeted_cavitation,
                rbc_venturi_protection,
            );
            prop_assert!(index >= 0.0);
            prop_assert!(index <= 1.0);
        }

        #[test]
        fn healthy_cell_protection_index_is_monotone_in_rbc_protection(
            wbc_targeted_cavitation in 0.0f64..1.0,
            low_rbc in 0.0f64..1.0,
            high_rbc in 0.0f64..1.0,
        ) {
            let (rbc_lo, rbc_hi) = if low_rbc <= high_rbc {
                (low_rbc, high_rbc)
            } else {
                (high_rbc, low_rbc)
            };
            let lo = healthy_cell_protection_index(wbc_targeted_cavitation, rbc_lo);
            let hi = healthy_cell_protection_index(wbc_targeted_cavitation, rbc_hi);
            prop_assert!(lo <= hi + 1.0e-15);
        }

        #[test]
        fn healthy_cell_protection_index_is_monotone_in_wbc_sparing(
            low_wbc_targeted_cavitation in 0.0f64..1.0,
            high_wbc_targeted_cavitation in 0.0f64..1.0,
            rbc_venturi_protection in 0.0f64..1.0,
        ) {
            let (wbc_lo, wbc_hi) = if low_wbc_targeted_cavitation <= high_wbc_targeted_cavitation {
                (low_wbc_targeted_cavitation, high_wbc_targeted_cavitation)
            } else {
                (high_wbc_targeted_cavitation, low_wbc_targeted_cavitation)
            };
            let lo = healthy_cell_protection_index(wbc_hi, rbc_venturi_protection);
            let hi = healthy_cell_protection_index(wbc_lo, rbc_venturi_protection);
            prop_assert!(lo <= hi + 1.0e-15);
        }
    }
}
