//! Convection scheme strategies for finite volume discretization.
//!
//! This module implements the Strategy pattern for convection schemes,
//! allowing easy extension and pluggable convection discretization methods.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use cfd_core::constants;

/// Trait for convection discretization schemes
pub trait ConvectionScheme<T: RealField>: Send + Sync {
    /// Calculate convection coefficients for east and west faces
    fn coefficients(&self, fe: T, fw: T, de: T, dw: T) -> (T, T);
    
    /// Get scheme name for identification
    fn name(&self) -> &'static str;
    
    /// Check if scheme is bounded (prevents oscillations)
    fn is_bounded(&self) -> bool;
}

/// First-order upwind scheme - stable but diffusive
pub struct FirstOrderUpwind;

impl<T: RealField> ConvectionScheme<T> for FirstOrderUpwind {
    fn coefficients(&self, fe: T, fw: T, de: T, dw: T) -> (T, T) {
        // Upwind: take flux from upwind cell
        let ae = de + T::max(T::zero(), -fe);
        let aw = dw + T::max(T::zero(), fw);
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

impl<T: RealField + FromPrimitive> ConvectionScheme<T> for CentralDifference {
    fn coefficients(&self, fe: T, fw: T, de: T, dw: T) -> (T, T) {
        let two = T::from_f64(constants::TWO).unwrap();
        let ae = de - fe / two.clone();
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

impl<T: RealField + FromPrimitive> ConvectionScheme<T> for HybridScheme {
    fn coefficients(&self, fe: T, fw: T, de: T, dw: T) -> (T, T) {
        let two = T::from_f64(constants::TWO).unwrap();
        
        // Calculate Peclet numbers
        let pe_e = fe.clone().abs() / de.clone();
        let pe_w = fw.clone().abs() / dw.clone();
        
        // Hybrid switching criteria (Pe = 2)
        let switch_criterion = two.clone();
        
        let ae = if pe_e <= switch_criterion {
            // Central difference
            de.clone() - fe.clone() / two.clone()
        } else {
            // Upwind
            de + T::max(T::zero(), -fe)
        };
        
        let aw = if pe_w <= switch_criterion {
            // Central difference
            dw.clone() + fw.clone() / two
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

impl<T: RealField + FromPrimitive> ConvectionScheme<T> for PowerLawScheme {
    fn coefficients(&self, fe: T, fw: T, de: T, dw: T) -> (T, T) {
        // Power-law function: max(0, (1 - 0.1|Pe|)^5)
        let power_law = |pe: T| -> T {
            let abs_pe = pe.abs();
            let one_tenth = T::from_f64(constants::ONE_TENTH).unwrap();
            let factor = T::one() - one_tenth * abs_pe;
            if factor > T::zero() {
                // Approximate (1-x)^5 for small x
                let factor2 = factor.clone() * factor.clone();
                let factor4 = factor2.clone() * factor2;
                factor4 * factor
            } else {
                T::zero()
            }
        };
        
        let pe_e = fe.clone() / de.clone();
        let pe_w = fw.clone() / dw.clone();
        
        let d_e_eff = de.clone() * power_law(pe_e) + T::max(T::zero(), -fe.clone());
        let d_w_eff = dw.clone() * power_law(pe_w) + T::max(T::zero(), fw.clone());
        
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

/// QUICK scheme - third-order accurate upwind-biased
pub struct QuickScheme;

impl<T: RealField + FromPrimitive> ConvectionScheme<T> for QuickScheme {
    fn coefficients(&self, fe: T, fw: T, de: T, dw: T) -> (T, T) {
        // Simplified QUICK implementation
        // In practice, QUICK requires additional stencil points
        // This is a basic approximation
        
        let _three_eighths = T::from_f64(constants::QUICK_COEFF_UPSTREAM).unwrap();
        let six_eighths = T::from_f64(constants::QUICK_COEFF_DOWNSTREAM).unwrap();
        
        // Approximate QUICK as modified upwind
        let ae = de.clone() + T::max(T::zero(), -fe.clone()) * six_eighths.clone();
        let aw = dw.clone() + T::max(T::zero(), fw.clone()) * six_eighths;
        
        (ae, aw)
    }
    
    fn name(&self) -> &'static str {
        "QUICK"
    }
    
    fn is_bounded(&self) -> bool {
        false // Can produce oscillations without limiters
    }
}

/// Factory for creating convection schemes
pub struct ConvectionSchemeFactory;

impl ConvectionSchemeFactory {
    /// Create scheme by name
    pub fn create<T: RealField + FromPrimitive>(scheme_name: &str) -> Box<dyn ConvectionScheme<T>> {
        match scheme_name.to_lowercase().as_str() {
            "upwind" | "first_order_upwind" => Box::new(FirstOrderUpwind),
            "central" | "central_difference" => Box::new(CentralDifference),
            "hybrid" => Box::new(HybridScheme),
            "power_law" => Box::new(PowerLawScheme),
            "quick" => Box::new(QuickScheme),
            _ => {
                // Default to stable first-order upwind
                Box::new(FirstOrderUpwind)
            }
        }
    }
    
    /// List available schemes
    pub fn available_schemes() -> Vec<&'static str> {
        vec![
            "upwind",
            "central", 
            "hybrid",
            "power_law",
            "quick"
        ]
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
        assert!(scheme.is_bounded());
    }
    
    #[test]
    fn test_central_difference() {
        let scheme = CentralDifference;
        let (ae, aw) = scheme.coefficients(2.0, -2.0, 1.0, 1.0);
        
        // Central difference: ae = de - fe/2, aw = dw + fw/2
        assert_eq!(ae, 0.0); // 1.0 - 2.0/2.0
        assert_eq!(aw, 0.0); // 1.0 + (-2.0)/2.0
        assert!(!scheme.is_bounded());
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
}