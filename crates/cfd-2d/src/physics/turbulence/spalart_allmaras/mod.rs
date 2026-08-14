//! # Spalart-Allmaras One-Equation Turbulence Model
//!
//! ## Mathematical Foundation
//!
//! The Spalart-Allmaras (SA) model is a one-equation turbulence model that solves
//! a transport equation for a modified turbulent kinematic viscosity $\tilde{\nu}$.
//!
//! **Reference:** Spalart, P. R., & Allmaras, S. R. (1994). "A one-equation turbulence
//! model for aerodynamic flows." *AIAA Journal*, 32(8), 1598-1605.
//!
//! ### Governing Transport Equation
//!
//! The SA model solves a single transport equation for the modified turbulent viscosity:
//!
//! ```math
//! \frac{D\tilde{\nu}}{Dt} = C_{b1} \tilde{S} \tilde{\nu} - C_{w1} f_w \left( \frac{\tilde{\nu}}{\tilde{d}} \right)^2 + \frac{1}{\sigma} \left[ \nabla \cdot \left( (\nu + \tilde{\nu}) \nabla \tilde{\nu} \right) + C_{b2} |\nabla \tilde{\nu}|^2 \right]
//! ```
//!
//! where:
//! - $\tilde{S}$ is the modified vorticity magnitude
//! - $\tilde{d}$ is the distance to the nearest wall
//! - $`f_w`$ is the wall destruction function
//!
//! ### Turbulent Viscosity Calculation
//!
//! The eddy viscosity is computed from the modified viscosity using damping functions:
//!
//! ```math
//! \nu_t = \tilde{\nu} f_{v1}
//! ```
//!
//! where the damping function $f_{v1}$ is:
//!
//! ```math
//! f_{v1} = \frac{\chi^3}{\chi^3 + C_{v1}^3}, \quad \chi = \frac{\tilde{\nu}}{\nu}
//! ```
//!
//! ### Modified Vorticity Magnitude
//!
//! The modified vorticity $\tilde{S}$ includes rotation and strain rate effects:
//!
//! ```math
//! \tilde{S} = \Omega + \frac{\tilde{\nu}}{\kappa^2 \tilde{d}^2} f_{v2}
//! ```
//!
//! where:
//! - $\Omega = \sqrt{2 \Omega_{ij} \Omega_{ij}}$ is the vorticity magnitude
//! - $\kappa = 0.41$ is the von Kármán constant
//! - $f_{v2}$ is the secondary damping function:
//!
//! ```math
//! f_{v2} = 1 - \frac{\chi}{1 + \chi f_{v1}}
//! ```
//!
//! ## Model Functions and Terms
//!
//! ### Production Term
//!
//! The production of turbulent kinetic energy is modeled as:
//!
//! ```math
//! P = C_{b1} \tilde{S} \tilde{\nu}
//! ```
//!
//! ### Wall Destruction Function
//!
//! The wall destruction function $`f_w`$ provides the correct near-wall behavior:
//!
//! ```math
//! f_w = g \left[ \frac{1 + C_{w3}^6}{g^6 + C_{w3}^6} \right]^{1/6}
//! ```
//!
//! where:
//!
//! ```math
//! g = r + C_{w2} (r^6 - r), \quad r = \frac{\tilde{\nu}}{\tilde{S} \kappa^2 \tilde{d}^2}
//! ```
//!
//! ### Destruction Term
//!
//! The destruction term ensures decay near walls:
//!
//! ```math
//! D = C_{w1} f_w \left( \frac{\tilde{\nu}}{\tilde{d}} \right)^2
//! ```
//!
//! ## Model Constants
//!
//! ### Standard SA Constants (Spalart & Allmaras, 1994)
//!
//! | Constant | Value | Description |
//! |----------|-------|-------------|
//! | $C_{b1}$ | 0.1355 | Production coefficient |
//! | $C_{b2}$ | 0.622 | Cross-diffusion coefficient |
//! | $C_{v1}$ | 7.1 | Primary damping coefficient |
//! | $C_{v2}$ | 0.3 | Secondary damping coefficient |
//! | $C_{w1}$ | $\frac{C_{b1}}{\kappa^2} + \frac{1 + C_{b2}}{\sigma}$ | Destruction coefficient |
//! | $C_{w2}$ | 0.3 | Wall destruction coefficient |
//! | $C_{w3}$ | 2.0 | Wall destruction coefficient |
//! | $\sigma$ | $\frac{2}{3}$ | Turbulent Prandtl number |
//! | $\kappa$ | 0.41 | von Kármán constant |
//!
//! ### Trip Terms (Optional)
//!
//! For transition modeling, additional trip terms can be included:
//!
//! ```math
//! F_{t1} = C_{t1} g_t \exp\left( -C_{t2} \frac{\omega_t^2}{\Delta U^2} - C_{t3} \sqrt{\frac{\omega_t}{\Delta U}} \tilde{d}^2 \right)
//! ```
//!
//! where $C_{t1} = 1$, $C_{t2} = 2$, $C_{t3} = 1.2$, $C_{t4} = 0.5$.
//!
//! ## Boundary Conditions
//!
//! ### Wall Boundary Conditions
//!
//! At solid walls, the modified turbulent viscosity is set to zero:
//!
//! ```math
//! \tilde{\nu}|_{wall} = 0
//! ```
//!
//! ### Inlet/Outlet Conditions
//!
//! Inlet conditions are specified based on experimental data or freestream values.
//! Outlet conditions use zero-gradient (Neumann) boundary conditions.
//!
//! ### Far-Field Conditions
//!
//! In the freestream, $\tilde{\nu}$ approaches a small positive value to maintain
//! numerical stability while ensuring negligible turbulent viscosity.
//!
//! ## Numerical Implementation
//!
//! ### Discretization
//!
//! The transport equation is discretized using finite differences:
//!
//! ```math
//! \frac{\tilde{\nu}^{n+1} - \tilde{\nu}^n}{\Delta t} = P - D + \frac{1}{\sigma} \nabla \cdot [(\nu + \tilde{\nu}) \nabla \tilde{\nu}] + \frac{C_{b2}}{\sigma} |\nabla \tilde{\nu}|^2
//! ```
//!
//! ### Stability Considerations
//!
//! 1. **Time step limitation**: $\Delta t \leq \frac{\Delta x^2}{2(\nu + \tilde{\nu})_{max}}$
//! 2. **Minimum values**: $\tilde{\nu}_{min} = 0$ (enforced at walls)
//! 3. **Wall distance calculation**: Accurate $\tilde{d}$ computation is crucial
//! 4. **Cross-diffusion limiting**: Prevents numerical oscillations
//!
//! ## Validation and Accuracy
//!
//! ### Theoretical Validation
//!
//! The SA model has been extensively validated against:
//! - Flat plate boundary layers
//! - Airfoil flows (NACA 0012, RAE 2822)
//! - Backward-facing step flows
//! - Jet flows and mixing layers
//! - High-lift configurations
//!
//! ### Key Strengths
//!
//! 1. **Computational efficiency**: Single transport equation
//! 2. **Wall treatment**: No wall functions required
//! 3. **Robustness**: Stable for wide range of flow conditions
//! 4. **Accuracy**: Good prediction of separation and reattachment
//! 5. **Convergence**: Fast convergence properties
//!
//! ### Limitations
//!
//! 1. **Freestream dependence**: Requires careful specification of freestream $\tilde{\nu}$
//! 2. **3D flows**: May require rotation/curvature corrections
//! 3. **Compressibility**: Needs modifications for high-speed flows
//! 4. **Transition**: No built-in transition modeling (requires extensions)
//!
//! ## Implementation Notes
//!
//! This implementation provides:
//! - Complete SA transport equation with all terms
//! - Accurate wall distance calculation for 2D domains
//! - Proper boundary condition enforcement
//! - Numerical stability safeguards
//! - Comprehensive validation tests
//!
//! The SA model is particularly well-suited for:
//! - Aerospace applications (airfoils, wings, nacelles)
//! - External aerodynamics with complex geometries
//! - Flows with massive separation
//! - Industrial applications requiring robust turbulence modeling
//!
//! ## Extensions and Variants
//!
//! ### SA with Rotation/Curvature Correction (SARC)
//!
//! Adds terms to account for rotational and curvature effects:
//!
//! ```math
//! \frac{D\tilde{\nu}}{Dt} = \dots + C_{rot} \frac{\tilde{\nu}}{\tilde{d}^2} f_{r1} |\nabla \times \vec{U}|
//! ```
//!
//! ### Negative SA (SA-)
//!
//! Modified version with improved freestream behavior:
//!
//! ```math
//! \nu_t = \tilde{\nu} f_{v1}, \quad f_{v1} = \frac{\chi^3}{\chi^3 + C_{v1}^3} \cdot \frac{1 + \chi/C_{v1}}{1 + \chi}
//! ```
//!
//! ### SA with Quadratic Constitutive Relation (SA-QCR)
//!
//! Includes non-linear eddy viscosity terms for improved anisotropy prediction.
//!
//! ## References
//!
//! - Spalart, P. R., & Allmaras, S. R. (1992). "A one-equation turbulence model for aerodynamic flows." In *30th Aerospace Sciences Meeting* (p. 439).
//! - Spalart, P. R., & Allmaras, S. R. (1994). "A one-equation turbulence model for aerodynamic flows." *AIAA Journal*, 32(8), 1598-1605.
//! - Allmaras, S. R., & Johnson, F. T. (2012). "Modifications and clarifications for the implementation of the Spalart-Allmaras turbulence model." In *7th International Conference on Computational Fluid Dynamics* (p. 189).
//!
//! # Theorem
//! The turbulence model must satisfy the realizability conditions for the Reynolds stress tensor.
//!
//! **Proof sketch**:
//! For any turbulent flow, the Reynolds stress tensor $\tau_{ij} = -\rho \`overline{u_i^\prime` `u_j^\prime`}$
//! must be positive semi-definite. This requires that the turbulent kinetic energy $k \ge 0$
//! and the normal stresses $\`overline{u_i^\prime` `u_i^\prime`} \ge 0$. The implemented model
//! enforces these constraints either through exact transport equations or bounded eddy-viscosity
//! formulations, ensuring physical realizability and numerical stability.

mod model;
mod solver;
mod wall_distance;

pub use model::SpalartAllmaras;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::turbulence::constants::{SA_CB1, SA_CW1};
    use crate::physics::turbulence::traits::TurbulenceModel;
    use eunomia::assert_relative_eq;

    #[test]
    fn test_spalart_allmaras_creation() {
        let model = SpalartAllmaras::<f64>::new(10, 10);
        assert_eq!(model.nx, 10);
        assert_eq!(model.ny, 10);
    }

    #[test]
    fn test_eddy_viscosity_zero_nu_tilde() {
        let model = SpalartAllmaras::<f64>::new(10, 10);
        let nu_t = model.eddy_viscosity(0.0, 1e-5);
        assert_eq!(nu_t, 0.0);
    }

    #[test]
    fn test_eddy_viscosity_large_nu_tilde() {
        let model = SpalartAllmaras::<f64>::new(10, 10);
        let nu = 1e-5;
        let nu_tilde = 1e-3;
        let nu_t = model.eddy_viscosity(nu_tilde, nu);

        // For χ >> 1, fv1 → 1, so νt → ν̃
        assert!(nu_t > 0.9 * nu_tilde);
        assert!(nu_t <= nu_tilde);
    }

    #[test]
    fn test_vorticity_magnitude() {
        let model = SpalartAllmaras::<f64>::new(10, 10);

        // Test case: uniform rotation
        let velocity_gradient = [[0.0, 1.0], [-1.0, 0.0]];
        let vorticity = model.vorticity_magnitude(&velocity_gradient);

        // Ω = |∂v/∂x - ∂u/∂y| = |-1 - 1| = 2
        assert_relative_eq!(vorticity, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_production_term() {
        let model = SpalartAllmaras::<f64>::new(10, 10);
        let nu_tilde = 1e-4;
        let s_tilde = 100.0;

        let production = model.production(nu_tilde, s_tilde);

        // P = Cb1 * S̃ * ν̃
        let expected = SA_CB1 * s_tilde * nu_tilde;
        assert_relative_eq!(production, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_destruction_term() {
        let model = SpalartAllmaras::<f64>::new(10, 10);
        let nu_tilde = 1e-4;
        let wall_distance = 0.01;
        let fw = 1.0;

        let destruction = model.destruction(nu_tilde, wall_distance, fw);

        // D = Cw1 * fw * (ν̃/d)²
        let ratio = nu_tilde / wall_distance;
        let expected = SA_CW1 * fw * ratio * ratio;
        assert_relative_eq!(destruction, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_dissipation_term_mapping() {
        let model = SpalartAllmaras::<f64>::new(10, 10);
        let nu_tilde = 1e-4;
        let wall_distance = 0.01;
        let dissipation = model.dissipation_term(nu_tilde, wall_distance);
        let ratio = nu_tilde / wall_distance;
        let expected = SA_CW1 * ratio * ratio;
        assert_relative_eq!(dissipation, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_cbrt_function() {
        // Test cube root computation via helper
        use super::wall_distance::cbrt;
        assert_relative_eq!(cbrt(8.0), 2.0, epsilon = 1e-8);
        assert_relative_eq!(cbrt(27.0), 3.0, epsilon = 1e-8);
        assert_relative_eq!(cbrt(1.0), 1.0, epsilon = 1e-8);
        assert_eq!(cbrt(0.0), 0.0);
    }

    #[test]
    fn test_wall_distance_field() {
        let model = SpalartAllmaras::<f64>::new(5, 5);
        let distances = model.wall_distance_field(0.1, 0.1);

        assert_eq!(distances.len(), 25);

        // Corner should have minimum distance
        // Center should have maximum distance
        let center_idx = 2 * 5 + 2;
        let corner_idx = 0;

        assert!(distances[center_idx] > distances[corner_idx]);
        assert_relative_eq!(distances[corner_idx], 0.05, epsilon = 1e-10);
    }

    #[test]
    fn test_boundary_conditions() {
        let model = SpalartAllmaras::<f64>::new(5, 5);
        let mut nu_tilde = vec![1.0; 25];

        model.apply_boundary_conditions(&mut nu_tilde);

        // Check walls are zero
        for i in 0..5 {
            assert_eq!(nu_tilde[i], 0.0); // Bottom
            assert_eq!(nu_tilde[4 * 5 + i], 0.0); // Top
        }
        for j in 0..5 {
            assert_eq!(nu_tilde[j * 5], 0.0); // Left
            assert_eq!(nu_tilde[j * 5 + 4], 0.0); // Right
        }
    }

    #[test]
    fn test_trip_term_zero() {
        let model = SpalartAllmaras::<f64>::new(10, 10);
        let trip = model.trip_term(1e-4, 0.01);

        // For fully turbulent flows, trip term should be zero
        assert_eq!(trip, 0.0);
    }
}
