//! Richardson extrapolation for grid-independent solutions
//!
//! Implements Richardson extrapolation following ASME V&V 20-2009 guidelines.
//!
//! This is the single consolidated Richardson implementation (CFDRS-GA-012):
//! the duplicate in `manufactured::richardson::core` was retired, and its
//! numerical-stability guards (signed convergence-ratio rejection, order
//! bounds, `r^p ≈ 1` checks) were absorbed here behind typed `cfd-core`
//! errors.

use crate::scalar;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};

/// Richardson extrapolation calculator
///
/// Estimates grid-independent solutions using systematic grid refinement
#[derive(Debug, Clone)]
pub struct RichardsonExtrapolation<T: RealField + Copy> {
    /// Assumed or computed order of accuracy
    pub order: T,
    /// Grid refinement ratio (r = `h_coarse` / `h_fine`)
    pub refinement_ratio: T,
}

impl<T: RealField + Copy + FloatElement> RichardsonExtrapolation<T> {
    /// Create a new Richardson extrapolation with known order
    pub fn with_order(order: T, refinement_ratio: T) -> Result<Self> {
        if order <= scalar::zero::<T>() {
            return Err(Error::InvalidInput(
                "Order of accuracy must be positive".to_string(),
            ));
        }

        if refinement_ratio <= scalar::one::<T>() {
            return Err(Error::InvalidInput(
                "Refinement ratio must be greater than 1".to_string(),
            ));
        }

        Ok(Self {
            order,
            refinement_ratio,
        })
    }

    /// Create with standard second-order accuracy
    pub fn second_order(refinement_ratio: T) -> Result<Self> {
        let two = <T as FloatElement>::from_f64(2.0);
        Self::with_order(two, refinement_ratio)
    }

    /// Extrapolate to zero grid spacing using two solutions
    ///
    /// # Arguments
    /// * `f_fine` - Solution on fine grid
    /// * `f_coarse` - Solution on the grid with spacing `r * h_fine`
    ///
    /// # Returns
    /// Extrapolated solution at h→0
    ///
    /// # Errors
    /// [`Error::Numerical`] when `r^p ≈ 1`, where the extrapolation
    /// denominator vanishes and the estimate is numerically meaningless.
    pub fn extrapolate(&self, f_fine: T, f_coarse: T) -> Result<T> {
        let r_p = scalar::powf(self.refinement_ratio, self.order);
        let denominator = r_p - scalar::one::<T>();
        if scalar::abs(denominator) < <T as FloatElement>::from_f64(1e-8) {
            return Err(Error::Numerical(
                cfd_core::error::NumericalErrorKind::DivisionByZero,
            ));
        }
        Ok((r_p * f_fine - f_coarse) / denominator)
    }

    /// Compute grid convergence index (GCI) following Roache (1998)
    ///
    /// GCI provides an error band for the grid-converged solution. Note this
    /// is the absolute error band; the fractional form reported in the
    /// literature divides by `f_fine` at the call site.
    ///
    /// # Errors
    /// [`Error::Numerical`] when `r^p ≈ 1` (vanishing GCI denominator).
    pub fn grid_convergence_index(&self, f_fine: T, f_coarse: T, safety_factor: T) -> Result<T> {
        let r_p = scalar::powf(self.refinement_ratio, self.order);
        let denominator = r_p - scalar::one::<T>();
        if scalar::abs(denominator) < <T as FloatElement>::from_f64(1e-8) {
            return Err(Error::Numerical(
                cfd_core::error::NumericalErrorKind::DivisionByZero,
            ));
        }
        let epsilon = scalar::abs(f_fine - f_coarse);

        Ok(safety_factor * epsilon / denominator)
    }

    /// Estimate order of accuracy from three grid solutions
    ///
    /// Uses the generalized Richardson extrapolation formula
    pub fn estimate_order(f_coarse: T, f_medium: T, f_fine: T, refinement_ratio: T) -> Result<T>
    where
        T: RealField + Copy + FloatElement,
    {
        let epsilon_21 = f_medium - f_fine;
        let epsilon_32 = f_coarse - f_medium;

        let epsilon_tolerance = <T as FloatElement>::from_f64(
            cfd_core::physics::constants::numerical::solver::EPSILON_TOLERANCE,
        );

        if scalar::abs(epsilon_21) < epsilon_tolerance {
            return Err(Error::InvalidInput(
                "Solutions too close to estimate order".to_string(),
            ));
        }

        let ratio = epsilon_32 / epsilon_21;
        // Reject oscillatory (sign-alternating) or non-finite convergence:
        // the order formula is only meaningful for monotone error sequences.
        if ratio <= scalar::zero::<T>() || !NumericElement::is_finite(ratio) {
            return Err(Error::InvalidInput("Invalid convergence ratio".to_string()));
        }

        let order = scalar::ln(ratio) / scalar::ln(refinement_ratio);

        // Guard the usable CFD order range (absorbed from the retired
        // manufactured duplicate): NaN or out-of-bounds orders indicate
        // numerically meaningless estimates, not real convergence behavior.
        if !NumericElement::is_finite(order)
            || order < <T as FloatElement>::from_f64(0.1)
            || order > <T as FloatElement>::from_f64(15.0)
        {
            return Err(Error::InvalidInput(format!(
                "Richardson extrapolation numerically unstable: order {} out of bounds",
                <T as NumericElement>::to_f64(order)
            )));
        }

        Ok(order)
    }

    /// Check if solutions are in asymptotic range
    ///
    /// Returns true if the convergence ratio is consistent with the expected order
    pub fn is_asymptotic(&self, f_coarse: T, f_medium: T, f_fine: T) -> bool {
        let epsilon_21 = scalar::abs(f_medium - f_fine);
        let epsilon_32 = scalar::abs(f_coarse - f_medium);

        let epsilon_tolerance = <T as FloatElement>::from_f64(
            cfd_core::physics::constants::numerical::solver::EPSILON_TOLERANCE,
        );

        if epsilon_21 < epsilon_tolerance {
            return false;
        }

        let observed_ratio = epsilon_32 / epsilon_21;
        let expected_ratio = scalar::powf(self.refinement_ratio, self.order);

        // Check if within 10% of expected ratio
        let relative_diff = scalar::abs((observed_ratio - expected_ratio) / expected_ratio);
        let ten_percent = <T as FloatElement>::from_f64(0.1);
        relative_diff < ten_percent
    }
}

/// Perform Richardson extrapolation with automatic order estimation
///
/// Uses three grid levels to estimate order and extrapolate
pub fn richardson_extrapolate<T>(solutions: &[T], grid_sizes: &[T]) -> Result<(T, T)>
where
    T: RealField + Copy + FloatElement,
{
    if solutions.len() < 2 || solutions.len() != grid_sizes.len() {
        return Err(Error::InvalidInput(
            "Need at least 2 solutions with corresponding grid sizes".to_string(),
        ));
    }

    // Sort by grid size (finest first)
    let mut paired: Vec<_> = grid_sizes
        .iter()
        .zip(solutions.iter())
        .map(|(h, f)| (*h, *f))
        .collect();
    paired.sort_by(|a, b| a.0.partial_cmp(&b.0).expect("expected value"));

    let f_fine = paired[0].1;
    let f_coarse = paired[1].1;
    let h_fine = paired[0].0;
    let h_coarse = paired[1].0;

    let r21 = h_coarse / h_fine; // Refinement ratio between fine and medium grids

    // Estimate order if we have 3 or more solutions
    let order = if solutions.len() >= 3 {
        let h_medium = paired[1].0;
        let h_coarse_actual = paired[2].0;
        let r32 = h_coarse_actual / h_medium; // Refinement ratio between medium and coarse

        // Check if refinement ratios are uniform (within 1% tolerance)
        let one_percent = <T as FloatElement>::from_f64(0.01);

        if (scalar::abs(r21 - r32) / r21) > one_percent {
            return Err(Error::InvalidInput(format!(
                "Non-uniform grid refinement ratios ({r21:?} and {r32:?}) detected. \
                 Richardson extrapolation requires constant refinement ratio."
            )));
        }

        let f_medium = paired[1].1;
        let f_coarse_actual = paired[2].1;
        RichardsonExtrapolation::estimate_order(f_coarse_actual, f_medium, f_fine, r21)?
    } else {
        // Assume second order if not enough data
        <T as FloatElement>::from_f64(2.0)
    };

    let extrapolator = RichardsonExtrapolation::with_order(order, r21)?;
    let extrapolated = extrapolator.extrapolate(f_fine, f_coarse)?;

    Ok((extrapolated, order))
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    #[test]
    fn test_richardson_second_order() {
        // Test with exact second-order convergence
        let extrapolator =
            RichardsonExtrapolation::<f64>::second_order(2.0).expect("expected value");

        // Solutions: f(h) = 1 + h²
        let f_fine = 1.0 + 0.01; // h = 0.1
        let f_coarse = 1.0 + 0.04; // h = 0.2

        let extrapolated = extrapolator
            .extrapolate(f_fine, f_coarse)
            .expect("expected value");
        assert_relative_eq!(extrapolated, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_order_estimation() {
        // Test order estimation with known convergence
        let f_coarse = 1.16; // h = 0.4, f = 1 + h²
        let f_medium = 1.04; // h = 0.2
        let f_fine = 1.01; // h = 0.1

        let order = RichardsonExtrapolation::<f64>::estimate_order(f_coarse, f_medium, f_fine, 2.0)
            .expect("expected value");

        assert_relative_eq!(order, 2.0, epsilon = 0.01);
    }

    #[test]
    fn test_gci_calculation() {
        let extrapolator =
            RichardsonExtrapolation::<f64>::second_order(2.0).expect("expected value");
        let f_fine = 1.01;
        let f_coarse = 1.04;
        let safety_factor = 1.25; // Recommended for 3+ grids

        let gci = extrapolator
            .grid_convergence_index(f_fine, f_coarse, safety_factor)
            .expect("expected value");

        // GCI should be small for well-converged solutions
        assert!(gci < 0.02);
    }

    /// Verbatim three-grid worked example from Roache (1998) as published in
    /// the NASA GRC "Examining Spatial (Grid) Convergence" tutorial
    /// (supersonic-inlet pressure recovery, grids at r = 2 with solutions
    /// f = (0.97050, 0.96854, 0.96178) fine -> coarse).
    ///
    /// Published reference values: p = 1.786170 (hand calc) / 1.78618479
    /// (NASA's VERIFY program), f_h->0 = 0.97130 / 0.971300304, fractional
    /// GCI_12 = 0.103083% / 0.103080%, GCI_23 = 0.356249% / 0.356244%,
    /// asymptotic ratio ~= 1.002.
    ///
    /// For this data the algebra is exactly rational: |e32/e21| =
    /// 0.00676/0.00196 = 169/49, hence the asserted anchors are exact —
    /// 2^p = 169/49, f_h->0 = 0.97050 + 0.00196*49/120, GCI_12 =
    /// 1.25*0.00196*49/120/0.97050 — with the published (rounded-p) values
    /// checked at looser tolerance to absorb their rounding.
    #[test]
    fn test_roache_1998_three_grid_worked_example() {
        let f_fine = 0.97050;
        let f_medium = 0.96854;
        let f_coarse = 0.96178;
        let r = 2.0;

        // Observed order from the three-grid formula.
        let p = RichardsonExtrapolation::<f64>::estimate_order(f_coarse, f_medium, f_fine, r)
            .expect("expected value");

        // Exact anchor: 2^p = 169/49 (|e32/e21| = 676/196 = 169/49); the
        // decimal published by NASA's hand calculation is 1.786170.
        let two_to_p_exact = 169.0 / 49.0;
        assert!(
            (p - 1.786_170).abs() < 1e-5,
            "observed order {p} must match published 1.786170"
        );
        assert!(
            (r.powf(p) - two_to_p_exact).abs() < 1e-9,
            "2^p must equal 169/49, got {}",
            r.powf(p)
        );

        // Richardson extrapolation on the two finest grids. The exact-p
        // anchor is 0.97050 + 0.00196*49/120 = 0.971300333...; NASA's VERIFY
        // program printed 0.971300304 from its slightly rounded order.
        let extrapolator = RichardsonExtrapolation::with_order(p, r).expect("expected value");
        let f_h0 = extrapolator
            .extrapolate(f_fine, f_medium)
            .expect("expected value");
        assert!(
            (f_h0 - 0.97050 - 0.00196 * 49.0 / 120.0).abs() < 1e-9,
            "extrapolated value {f_h0} must match exact anchor 0.971300333"
        );
        assert!(
            (f_h0 - 0.971_300_304).abs() < 1e-7,
            "extrapolated value {f_h0} must match published 0.971300304"
        );

        // Fractional GCI on the fine grid with the three-grid safety factor.
        let fs = 1.25;
        let gci_abs = extrapolator
            .grid_convergence_index(f_fine, f_medium, fs)
            .expect("expected value");
        let gci_fractional_percent = gci_abs / f_fine * 100.0;
        assert!(
            (gci_fractional_percent - 0.103_080).abs() < 2e-5,
            "GCI_12 {gci_fractional_percent}% must match published 0.103080%"
        );

        // Roache's asymptotic-range check: GCI_23 ~= r^p * GCI_12.
        let gci_coarse_abs = extrapolator
            .grid_convergence_index(f_medium, f_coarse, fs)
            .expect("expected value");
        let gci_coarse_percent = gci_coarse_abs / f_medium * 100.0;
        assert!(
            (gci_coarse_percent - 0.356_244).abs() < 5e-5,
            "GCI_23 {gci_coarse_percent}% must match published 0.356244%"
        );
        let asymptotic_ratio = gci_coarse_percent / gci_fractional_percent / two_to_p_exact;
        assert!(
            (asymptotic_ratio - 1.0).abs() < 0.01,
            "GCI_23 / (r^p GCI_12) = {asymptotic_ratio} must be ~1 (asymptotic range)"
        );

        // The ratio-band asymptotic check also flags this data as asymptotic.
        assert!(extrapolator.is_asymptotic(f_coarse, f_medium, f_fine));
    }
}
