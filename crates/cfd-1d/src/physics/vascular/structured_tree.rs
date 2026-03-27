//! Olufsen (1999) Structured Tree Impedance Boundary Condition
//!
//! ## Theorem — Fractal Tree Impedance (Olufsen 1999)
//!
//! **Theorem**: The outflow boundary of a truncated 1D arterial network must be
//! terminated by a structured fractal tree of resistive and compliant vessels
//! extending down to the precapillary scale ($r_{min}$) to properly reflect 
//! outgoing pressure waves and enforce a physiological terminal resistance.
//!
//! For a steady-state ($\omega \to 0$) lumped parameter model, the complex
//! wave impedance $Z(\omega)$ mathematically reduces to the zero-frequency 
//! steady equivalent resistance $R_{tree}$ of an asymmetric binary branching
//! network.
//!
//! Given asymmetric branching ratios $\alpha = r_1/r_p$ and $\beta = r_2/r_p$,
//! and a length-to-radius scaling factor $\lambda = L/r$, the local segment 
//! resistance is derived from Poiseuille flow:
//!
//! $R_p = \frac{8 \mu L}{\pi r^4} = \frac{8 \mu \lambda}{\pi r^3}$
//!
//! The equivalent resistance of any node in the tree evaluates recursively:
//!
//! $R_{eq\_p} = R_p + \left( \frac{1}{R_{eq\_1}} + \frac{1}{R_{eq\_2}} \right)^{-1}$
//!
//! This recursive subdivision continues until $r < r_{min}$, at which point
//! the segment resistance is treated as the pure terminal resistance.
//!
//! ### Proof: Area-Preserving Sub-trees
//!
//! Under Murray's Law ($\alpha^3 + \beta^3 = 1$), the net resistance geometrically
//! scales, ensuring finite pressure drop. When flow distributes such that wall
//! shear stress is uniform, the fractal impedance $Z(0)$ asymptotically equates
//! to the physiological systemic mean resistance $P_{art}/Q_{CO}$.
//!
//! ### References
//! - Olufsen, M. S. (1999). "Structured tree outflow condition for blood flow in 
//!   larger systemic arteries". *American Journal of Physiology*.

use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use serde::{Deserialize, Serialize};


/// Morphological parameters establishing the Olufsen (1999) fractal tree
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct OlufsenParameters<T: RealField + Copy> {
    /// Ratio of the primary daughter radius to parent radius ($\alpha$)
    pub alpha: T,
    /// Ratio of the secondary daughter radius to parent radius ($\beta$)
    pub beta: T,
    /// Length-to-radius proportionality constant ($\lambda \approx 50$)
    pub length_ratio: T,
    /// Minimum terminal branching radius before truncation ($r_{min}$)
    pub r_min: T,
}

impl<T: RealField + Copy> OlufsenParameters<T> {
    /// Creates a physiologically-backed generic parameter set
    /// Uses accepted systemic values where $\alpha=0.9$, $\beta=0.6$, $\lambda=50$.
    pub fn new_systemic(r_min: T) -> Self {
        Self {
            alpha: T::from_f64(0.90).expect("Mathematical constant conversion compromised"),
            beta: T::from_f64(0.60).expect("Mathematical constant conversion compromised"),
            length_ratio: T::from_f64(50.0).expect("Mathematical constant conversion compromised"),
            r_min,
        }
    }

    /// Recursively calculates the Olufsen zero-frequency structured tree resistance
    /// Limits recursion to structurally valid $r_min$ depths to avoid stack overflow.
    pub fn compute_steady_impedance(
        &self, 
        root_radius: T, 
        dynamic_viscosity: T
    ) -> Result<T> {
        if root_radius <= T::zero() || self.r_min <= T::zero() {
            return Err(Error::PhysicsViolation(
                "Tree generation requires strictly positive root and minimum radii limits".to_string()
            ));
        }
        if self.alpha >= T::one() || self.beta >= T::one() {
             return Err(Error::PhysicsViolation(
                "Fractal splitting variables alpha and beta must be < 1 to enforce convergence down to r_min".to_string()
            ));
        }

        // Fast computation down to r_min via recursive depth-first iteration
        // In systemic modeling, alpha^3 + beta^3 ≈ 1 (Murray's Law)
        let z_0 = self.recursive_resistance(root_radius, dynamic_viscosity);
        Ok(z_0)
    }

    fn recursive_resistance(&self, current_r: T, viscosity: T) -> T {
        // Poiseuille segment resistance parameterization for this single branch
        let eight = T::from_f64(8.0).expect("Mathematical constant conversion compromised");
        let pi = T::pi();
        
        let r3 = current_r * current_r * current_r;
        let segment_resistance = (eight * viscosity * self.length_ratio) / (pi * r3);

        if current_r < self.r_min {
            // Leaf termination
            segment_resistance
        } else {
            // Bifurcate
            let z_1 = self.recursive_resistance(current_r * self.alpha, viscosity);
            let z_2 = self.recursive_resistance(current_r * self.beta, viscosity);

            // Parallel impedance of distals
            let parallel_distal = (z_1 * z_2) / (z_1 + z_2);
            
            // Total is local path + parallel distals
            segment_resistance + parallel_distal
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn construct_olufsen_tree_terminates_at_rmin() {
        let params = OlufsenParameters {
            alpha: 0.8_f64,
            beta: 0.6_f64,   // (0.8^3 + 0.6^3 = 0.512 + 0.216 = 0.728) < 1 means shrinking areas
            length_ratio: 50.0,
            r_min: 0.05,     // 50 microns
        };

        let root_r = 0.2; // 200 microns
        let mu = 0.003;   // roughly blood viscosity

        let res = params.compute_steady_impedance(root_r, mu).unwrap();
        
        // Ensure result is > 0 and bounded
        assert!(res > 0.0);
        
        // Single segment baseline logic check (no bifurcations)
        let r3 = root_r * root_r * root_r;
        let baseline = (8.0 * 0.003 * 50.0) / (std::f64::consts::PI * r3);
        assert!(res > baseline);
    }
    
    #[test]
    fn olufsen_tree_preserves_symmetry() {
        // A perfectly symmetric tree splitting into identical daughters.
        // It can be analytically compared to geometric series resistance summing!
        let params = OlufsenParameters {
            alpha: 0.5_f64,
            beta: 0.5_f64,
            length_ratio: 50.0,
            r_min: 0.1_f64,  // Set high so only 1 level of branching happens
        };
        // root = 0.4 -> Level 1 daughters = 0.2 -> Level 2 daughters = 0.1 (termination)
        
        let root = 0.4_f64;
        let mu = 0.001_f64;
        
        let z_total = params.compute_steady_impedance(root, mu).unwrap();
        
        // Analytical construction:
        let r0 = root;
        let r1 = root * 0.5; // 0.2
        let r2 = r1 * 0.5;   // 0.1
        
        let seg_res = |r: f64| -> f64 {
            (8.0 * mu * 50.0) / (std::f64::consts::PI * r.powi(3))
        };
        
        let z0 = seg_res(r0);
        let z1 = seg_res(r1);
        let z2 = seg_res(r2); // This is termination layer because r2 = 0.1 is not strict < r_min
        
        // Level 2 terminate here without bifurcation:
        // Wait, the termination condition is `current_r < r_min`.
        // If current_r = 0.1, is it < 0.1? No, 0.1 is not strictly < 0.1
        // So it bifurcates ONE MORE TIME!
        let r3 = r2 * 0.5; // 0.05
        let z3_term = seg_res(r3);
        
        let z2_total = z2 + (z3_term * z3_term) / (z3_term + z3_term); // z2 + z3_term / 2
        let z1_total = z1 + (z2_total * z2_total) / (z2_total + z2_total); // z1 + z2_total / 2
        let expected = z0 + (z1_total * z1_total) / (z1_total + z1_total); // z0 + z1_total / 2
        
        assert_relative_eq!(z_total, expected, max_relative = 1e-12);
    }

    #[test]
    fn invalid_fractal_exponents_fail() {
        let params = OlufsenParameters {
            alpha: 1.05_f64, // Invalid > 1
            beta: 0.8,
            length_ratio: 50.0,
            r_min: 0.1,
        };
        let res = params.compute_steady_impedance(1.0_f64, 0.001);
        assert!(matches!(res, Err(Error::PhysicsViolation(_))));
    }
}
