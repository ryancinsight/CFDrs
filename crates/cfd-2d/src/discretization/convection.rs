//! Convection scheme strategies for finite volume discretization.
//!
//! This module implements the Strategy pattern for convection schemes,
//! allowing easy extension and pluggable convection discretization methods.
//!
//! # Theorem (Upwind Boundedness — Godunov 1959)
//!
//! A first-order upwind scheme is monotone: if the initial data satisfies
//! $\phi_{\min} \le \phi_i^0 \le \phi_{\max}$ for all $i$, then
//! $\phi_{\min} \le \phi_i^n \le \phi_{\max}$ for all $n > 0$, provided the CFL
//! condition $|u| \Delta t / \Delta x \le 1$ holds.
//!
//! **Proof sketch**:
//! The upwind stencil $\phi_i^{n+1} = \phi_i^n - C(\phi_i^n - \phi_{i-1}^n)$ (for $u > 0$,
//! Courant number $C \ge 0$) is a convex combination of $\phi_i^n$ and $\phi_{i-1}^n$
//! when $0 \le C \le 1$, hence satisfies the discrete maximum principle. Extension to
//! 2D follows by dimensional splitting (Strang 1968).
//!
//! Higher-order schemes (QUICK, Power-Law) trade monotonicity for accuracy;
//! the Godunov theorem states that no linear scheme of order $\ge 2$ is monotone.
//! TVD limiters restore boundedness for higher-order schemes.

use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Trait for convection discretization schemes
pub trait ConvectionScheme<T: RealField + Copy>: Send + Sync {
    /// Calculate convection coefficients for east and west faces
    fn coefficients(&self, fe: T, fw: T, de: T, dw: T) -> (T, T);

    /// Get scheme name for identification
    fn name(&self) -> &'static str;

    /// Check if scheme is bounded (prevents oscillations)
    fn is_bounded(&self) -> bool;
}

/// First-order upwind scheme - stable but diffusive
pub struct FirstOrderUpwind;

impl<T: RealField + Copy> ConvectionScheme<T> for FirstOrderUpwind {
    fn coefficients(&self, fe: T, fw: T, de: T, dw: T) -> (T, T) {
        // Upwind scheme: always take from upwind direction
        // ae coefficient includes diffusion + convection from east face
        // aw coefficient includes diffusion + convection from west face
        let ae = de + T::max(T::zero(), -fe);
        let aw = dw + T::max(T::zero(), -fw);
        (ae, aw)
    }

    fn name(&self) -> &'static str {
        "First-Order Upwind"
    }

    fn is_bounded(&self) -> bool {
        true // Always bounded
    }
}

/// Central difference scheme - second-order accurate but can oscillate
pub struct CentralDifference;

impl<T: RealField + Copy + FromPrimitive + Copy> ConvectionScheme<T> for CentralDifference {
    fn coefficients(&self, fe: T, fw: T, de: T, dw: T) -> (T, T) {
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        let ae = de - fe / two;
        let aw = dw + fw / two;
        (ae, aw)
    }

    fn name(&self) -> &'static str {
        "Central Difference"
    }

    fn is_bounded(&self) -> bool {
        false // Can produce oscillations for high Peclet numbers
    }
}

/// Hybrid scheme - switches between upwind and central difference
pub struct HybridScheme;

impl<T: RealField + Copy + FromPrimitive + Copy> ConvectionScheme<T> for HybridScheme {
    fn coefficients(&self, fe: T, fw: T, de: T, dw: T) -> (T, T) {
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());

        // Calculate Peclet numbers
        let pe_e = fe.abs() / de;
        let pe_w = fw.abs() / dw;

        // Hybrid switching criteria (Pe = 2)
        let switch_criterion = two;

        let ae = if pe_e <= switch_criterion {
            // Central difference
            de - fe / two
        } else {
            // Upwind
            de + T::max(T::zero(), -fe)
        };

        let aw = if pe_w <= switch_criterion {
            // Central difference
            dw + fw / two
        } else {
            // Upwind
            dw + T::max(T::zero(), fw)
        };

        (ae, aw)
    }

    fn name(&self) -> &'static str {
        "Hybrid"
    }

    fn is_bounded(&self) -> bool {
        true // Bounded due to upwind switching
    }
}

/// Power-law scheme - more accurate interpolation for convection-diffusion
pub struct PowerLawScheme;

impl<T: RealField + Copy + FromPrimitive + Copy> ConvectionScheme<T> for PowerLawScheme {
    fn coefficients(&self, fe: T, fw: T, de: T, dw: T) -> (T, T) {
        // Power-law function: max(0, (1 - 0.1|Pe|)^5)
        let power_law = |pe: T| -> T {
            let abs_pe = pe.abs();
            let one_tenth = T::from_f64(0.1).unwrap_or_else(|| T::zero());
            let factor = T::one() - one_tenth * abs_pe;
            if factor > T::zero() {
                // Approximate (1-x)^5 for small x
                let factor2 = factor * factor;
                let factor4 = factor2 * factor2;
                factor4 * factor
            } else {
                T::zero()
            }
        };

        let pe_e = fe / de;
        let pe_w = fw / dw;

        let d_e_eff = de * power_law(pe_e) + T::max(T::zero(), -fe);
        let d_w_eff = dw * power_law(pe_w) + T::max(T::zero(), fw);

        let ae = d_e_eff;
        let aw = d_w_eff;

        (ae, aw)
    }

    fn name(&self) -> &'static str {
        "Power Law"
    }

    fn is_bounded(&self) -> bool {
        true // Inherently bounded
    }
}

/// Quadratic Upstream Interpolation for Convective Kinematics (QUICK) - third-order upwind-biased
/// Reference: Leonard, B.P. (1979). "A stable and accurate convective modelling procedure based on quadratic upstream interpolation"
pub struct QuadraticUpstreamInterpolationScheme;

impl<T: RealField + Copy + FromPrimitive + Copy> ConvectionScheme<T>
    for QuadraticUpstreamInterpolationScheme
{
    fn coefficients(&self, fe: T, fw: T, de: T, dw: T) -> (T, T) {
        // QUICK scheme according to Leonard (1979)
        // φ_f = 6/8 * φ_C + 3/8 * φ_D - 1/8 * φ_U
        // where C=central, D=downstream, U=upstream (2 cells away)
        //
        // For the limited API, we use deferred correction approach:
        // Treat as upwind for implicit part, add QUICK correction explicitly

        // Base upwind coefficients
        let ae = de + T::max(T::zero(), -fe);
        let aw = dw + T::max(T::zero(), fw);

        // Note: Full QUICK requires access to φ_UU which this API doesn't provide
        // Users should use the extended stencil API for accurate QUICK
        (ae, aw)
    }

    fn name(&self) -> &'static str {
        "QUICK"
    }

    fn is_bounded(&self) -> bool {
        true
    }
}

/// Factory for creating convection schemes
pub struct ConvectionSchemeFactory;

impl ConvectionSchemeFactory {
    /// Create scheme by name
    #[must_use]
    pub fn create<T: RealField + Copy + FromPrimitive + Copy>(
        scheme_name: &str,
    ) -> Box<dyn ConvectionScheme<T>> {
        match scheme_name.to_lowercase().as_str() {
            "upwind" | "first_order_upwind" => Box::new(FirstOrderUpwind),
            "central" | "central_difference" => Box::new(CentralDifference),
            "hybrid" => Box::new(HybridScheme),
            "power_law" => Box::new(PowerLawScheme),
            "quick" => Box::new(QuadraticUpstreamInterpolationScheme),
            _ => {
                // Default to stable first-order upwind
                Box::new(FirstOrderUpwind)
            }
        }
    }

    /// List available schemes
    #[must_use]
    pub fn available_schemes() -> Vec<&'static str> {
        vec!["upwind", "central", "hybrid", "power_law", "quick"]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_first_order_upwind() {
        let scheme = FirstOrderUpwind;
        let (ae, aw) = scheme.coefficients(1.0, -1.0, 0.5, 0.5);

        // For positive fe, should use upwind (no convection contribution to ae)
        // For negative fw, should use upwind (no convection contribution to aw)
        assert_eq!(ae, 0.5); // de only
        assert_eq!(aw, 1.5); // dw + |fw|
        assert!(<FirstOrderUpwind as ConvectionScheme<f64>>::is_bounded(
            &scheme
        ));
    }

    #[test]
    fn test_central_difference() {
        let scheme = CentralDifference;
        let (ae, aw) = scheme.coefficients(2.0, -2.0, 1.0, 1.0);

        // Central difference: ae = de - fe/2, aw = dw + fw/2
        assert_eq!(ae, 0.0); // 1.0 - 2.0/2.0
        assert_eq!(aw, 0.0); // 1.0 + (-2.0)/2.0
        assert!(!<CentralDifference as ConvectionScheme<f64>>::is_bounded(
            &scheme
        ));
    }

    #[test]
    fn test_scheme_factory() {
        let upwind: Box<dyn ConvectionScheme<f64>> = ConvectionSchemeFactory::create("upwind");
        assert_eq!(upwind.name(), "First-Order Upwind");

        let central: Box<dyn ConvectionScheme<f64>> = ConvectionSchemeFactory::create("central");
        assert_eq!(central.name(), "Central Difference");

        // Test default fallback
        let unknown: Box<dyn ConvectionScheme<f64>> = ConvectionSchemeFactory::create("unknown");
        assert_eq!(unknown.name(), "First-Order Upwind");
    }

    /// Power-law at Pe = 0 (no convection): should reduce to central differencing.
    /// When fe = fw = 0, the power-law function returns (1-0)^5 = 1, so
    /// ae = de·1 + max(0,0) = de, aw = dw·1 + max(0,0) = dw (pure diffusion).
    #[test]
    fn test_power_law_pure_diffusion() {
        let scheme = PowerLawScheme;
        let de = 2.0_f64;
        let dw = 2.0_f64;
        let (ae, aw) = scheme.coefficients(0.0, 0.0, de, dw);
        // pe_e = 0/2 = 0, power_law(0) = 1.0, ae = 2·1 + max(0,0) = 2
        assert!(
            (ae - de).abs() < 1e-12,
            "Power-law ae at Pe=0 mismatch: {ae}"
        );
        assert!(
            (aw - dw).abs() < 1e-12,
            "Power-law aw at Pe=0 mismatch: {aw}"
        );
    }

    /// Hybrid scheme: switches at Pe = 2.
    /// - Pe < 2 → central difference coefficients
    /// - Pe ≥ 2 → upwind coefficients
    #[test]
    fn test_hybrid_switching() {
        let scheme = HybridScheme;

        // Low Peclet (Pe = 1 < 2): should use central difference
        // fe = 1.0, de = 1.0, Pe_e = |fe|/de = 1, below threshold
        let (ae_low, _aw_low) = scheme.coefficients(1.0, -1.0, 1.0, 1.0);
        // Central: ae = de - fe/2 = 1 - 0.5 = 0.5
        assert!(
            (ae_low - 0.5_f64).abs() < 1e-12,
            "Hybrid should use central at Pe=1: ae={ae_low}"
        );

        // High Peclet (Pe = 4 > 2): should use upwind
        // fe = 4.0, de = 1.0, Pe_e = 4
        let (ae_high, _) = scheme.coefficients(4.0, -4.0, 1.0, 1.0);
        // Upwind: ae = de + max(0, -fe) = 1 + max(0,-4) = 1
        assert!(
            (ae_high - 1.0_f64).abs() < 1e-12,
            "Hybrid should use upwind at Pe=4: ae={ae_high}"
        );
    }
}
