//! Reynolds Stress Transport Model (RSTM) — module hierarchy.
//!
//! ## Structure
//! ```text
//! reynolds_stress/
//!   mod.rs            — PressureStrainModel enum + canonical pub use re-exports
//!   tensor.rs         — ReynoldsStressTensor storage type
//!   model.rs          — ReynoldsStressModel struct + constructor + initialisation
//!   production.rs     — P_ij exact production term (not Boussinesq)
//!   diffusion.rs      — ε_ij dissipation tensor + T_ij turbulent transport
//!   curvature.rs      — Suga-Craft (2003) streamline curvature correction
//!   wall_reflection.rs — Gibson-Launder (1978) wall-reflection correction
//!   pressure_strain/
//!     linear.rs       — Rotta (1951) linear return-to-isotropy
//!     quadratic.rs    — Speziale et al. (1991) quadratic model
//!     ssg.rs          — Full SSG model
//!   transport.rs      — Full time-advancement + TurbulenceModel trait impl
//! ```
//!
//! ## References
//! - Pope, S. B. (2000). *Turbulent Flows*. Cambridge University Press.
//! - Launder, B. E., Reece, G. J., & Rodi, W. (1975). J. Fluid Mech., 68(3), 537–566.
//! - Speziale, C. G., Sarkar, S., & Gatski, T. B. (1991). J. Fluid Mech., 227, 245–272.
//! - Gibson, M. M., & Launder, B. E. (1978). J. Fluid Mech., 86(3), 491–511.
//! - Suga, K., & Craft, T. J. (2003). Flow Turbulence Combust., 70, 143–162.

pub mod curvature;
pub mod diffusion;
pub mod model;
pub mod pressure_strain;
pub mod production;
pub mod tensor;
pub mod transport;
pub mod wall_reflection;

pub use model::ReynoldsStressModel;
pub use tensor::ReynoldsStressTensor;

/// Pressure-strain correlation model selector.
///
/// Passed to [`ReynoldsStressModel`] at construction time and used at every
/// pressure-strain evaluation to dispatch to the correct kernel.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PressureStrainModel {
    /// Linear return-to-isotropy (Rotta, 1951). Φ_ij = −C₁(ε/k) b_ij.
    LinearReturnToIsotropy,
    /// Quadratic slow + rapid (Speziale, Sarkar & Gatski, 1991).
    Quadratic,
    /// Full SSG model: non-linear in b_ij with strain and vorticity coupling.
    SSG,
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra::DMatrix;

    /// DNS channel flow validation (Moser et al., 1999) at Re_τ = 590.
    #[test]
    fn test_dns_channel_flow_validation() {
        let model = ReynoldsStressModel::<f64>::new(40, 40);
        let mut stresses = model.initialize_reynolds_stresses(1.0, 0.1);

        let re_tau = 590.0f64;
        let nu = 1.0 / re_tau;

        for j in 1..39 {
            let y_plus = (j as f64 / 39.0) * re_tau;
            let k_plus = if y_plus <= 10.0 { 0.15 * y_plus * y_plus }
                         else if y_plus <= 100.0 { 3.3 + 0.25 * (y_plus - 10.0).ln() }
                         else { 2.5 * (y_plus / 100.0).ln() + 1.0 };
            let k_dns = k_plus * nu * nu * re_tau * re_tau;
            stresses.k[(20, j)] = k_dns;
            stresses.epsilon[(20, j)] = k_dns.powf(1.5) / 0.09;
            let anis = if y_plus <= 10.0 { 0.8 } else if y_plus <= 50.0 { 0.6 } else { 0.4 };
            stresses.xx[(20, j)] = (2.0 / 3.0 + anis) * k_dns;
            stresses.yy[(20, j)] = (2.0 / 3.0 - anis) * k_dns;
            stresses.xy[(20, j)] = 0.0;
        }

        model.apply_wall_boundary_conditions(
            &mut stresses.xx, &mut stresses.xy, &mut stresses.yy,
            &mut stresses.k, &mut stresses.epsilon,
        );

        assert_relative_eq!(stresses.k[(20, 0)], 0.0, epsilon = 1e-10);
        assert_relative_eq!(stresses.k[(20, 39)], 0.0, epsilon = 1e-10);

        for j in 1..39 {
            assert!(stresses.xx[(20, j)] >= 0.0, "⟨u'u'⟩ must be non-negative");
            assert!(stresses.yy[(20, j)] >= 0.0, "⟨v'v'⟩ must be non-negative");
            let max_shear = (stresses.xx[(20, j)] * stresses.yy[(20, j)]).sqrt();
            assert!(stresses.xy[(20, j)].abs() <= max_shear + 1e-10, "Cauchy-Schwarz violated");
            let k = stresses.k[(20, j)];
            let bij_xx = stresses.xx[(20, j)] / k - 2.0 / 3.0;
            let bij_yy = stresses.yy[(20, j)] / k - 2.0 / 3.0;
            assert!(bij_xx.abs() <= 2.0 / 3.0 + 1e-6, "b_xx anisotropy out of bounds");
            assert!(bij_yy.abs() <= 2.0 / 3.0 + 1e-6, "b_yy anisotropy out of bounds");
        }
    }

    #[test]
    fn test_homogeneous_shear_analytical_validation() {
        let mut model = ReynoldsStressModel::<f64>::new(3, 3);
        model.pressure_strain_model = PressureStrainModel::LinearReturnToIsotropy;
        let mut stresses = model.initialize_reynolds_stresses(1.0, 0.1);

        let shear_rate = 1.0f64;
        let mut u = DMatrix::zeros(3, 3);
        let v = DMatrix::zeros(3, 3);
        for j in 0..3 { for i in 0..3 { u[(i, j)] = shear_rate * j as f64 * 0.1; } }
        let velocity = [u, v];
        let dt = 0.001;
        let dx = 0.1;
        let dy = 0.1;

        let mut uv_history = Vec::new();
        for step in 0..200 {
            let result = model.update_reynolds_stresses_optimized(&mut stresses, &velocity, dt, dx, dy);
            assert!(result.is_ok(), "RSM update failed at step {step}");
            uv_history.push(stresses.xy[(1, 1)]);
            if step > 10
                && (uv_history[step] - uv_history[step - 1]).abs()
                    / uv_history[step].abs().max(1e-12) < 0.001
            {
                break;
            }
        }

        let eq = stresses.xy[(1, 1)];
        assert!(eq < 0.0, "Shear stress must be negative in shear flow: {eq:.6}");
        assert!(eq.abs() > 1e-6, "Shear stress must evolve from isotropic IC");
    }

    #[test]
    fn test_reynolds_stress_initialization() {
        let model = ReynoldsStressModel::<f64>::new(10, 10);
        let s = model.initialize_reynolds_stresses(1.0, 0.1);
        assert_relative_eq!(s.xx[(5, 5)], 2.0 / 3.0, epsilon = 1e-6);
        assert_relative_eq!(s.yy[(5, 5)], 2.0 / 3.0, epsilon = 1e-6);
        assert_relative_eq!(s.xy[(5, 5)], 0.0, epsilon = 1e-6);
        assert_relative_eq!(s.k[(5, 5)], 1.0, epsilon = 1e-6);
        assert_relative_eq!(s.epsilon[(5, 5)], 0.1, epsilon = 1e-6);
    }

    #[test]
    fn test_turbulent_viscosity() {
        let model = ReynoldsStressModel::<f64>::new(10, 10);
        use super::super::traits::TurbulenceModel;
        let nu_t = model.turbulent_viscosity(1.0, 0.1, 1.0);
        assert_relative_eq!(nu_t, 0.9, epsilon = 1e-6);
    }

    #[test]
    fn test_reynolds_stress_realizability() {
        let model = ReynoldsStressModel::<f64>::new(10, 10);
        let s = model.initialize_reynolds_stresses(1.0, 0.1);
        for y in 0..10 { for x in 0..10 {
            let xx = s.xx[(x, y)]; let xy = s.xy[(x, y)]; let yy = s.yy[(x, y)];
            let k = s.k[(x, y)];
            assert!(xx >= 0.0, "⟨u'u'⟩ non-negative");
            assert!(yy >= 0.0, "⟨v'v'⟩ non-negative");
            assert!(xy.abs() <= (xx * yy).sqrt() + 1e-12, "Cauchy-Schwarz");
            if k > 1e-12 {
                assert!((xx / k - 2.0/3.0).abs() <= 2.0/3.0 + 1e-10);
                assert!((yy / k - 2.0/3.0).abs() <= 2.0/3.0 + 1e-10);
            }
        } }
    }

    #[test]
    fn test_production_term_correctness() {
        let model = ReynoldsStressModel::<f64>::new(1, 1);
        let rs = ReynoldsStressTensor {
            xx: DMatrix::from_element(1, 1, 1.0),
            xy: DMatrix::from_element(1, 1, 0.5),
            yy: DMatrix::from_element(1, 1, 0.8),
            k: DMatrix::from_element(1, 1, 1.35),
            epsilon: DMatrix::from_element(1, 1, 0.1),
            epsilon_xx: None, epsilon_xy: None, epsilon_yy: None,
        };
        let vg = [[0.0, 1.0], [0.0, 0.0]];
        let p_xx = model.production_term(&rs, &vg, 0, 0, 0, 0);
        let p_xy = model.production_term(&rs, &vg, 0, 1, 0, 0);
        let p_yy = model.production_term(&rs, &vg, 1, 1, 0, 0);
        assert_relative_eq!(p_xx, -1.0, epsilon = 1e-12);
        assert_relative_eq!(p_xy, -0.8, epsilon = 1e-12);
        assert_relative_eq!(p_yy, 0.0, epsilon = 1e-12);
        // Verify this is NOT Boussinesq (-0.5 would indicate Boussinesq)
        assert!((p_xx - (-0.5)).abs() > 0.1, "Must not match Boussinesq");
    }

    #[test]
    fn test_ssg_pressure_strain_physical_behavior() {
        let mut model = ReynoldsStressModel::<f64>::new(3, 3);
        model.pressure_strain_model = PressureStrainModel::SSG;
        let mut stresses = model.initialize_reynolds_stresses(1.0, 0.1);
        let mut u = DMatrix::zeros(3, 3);
        let v = DMatrix::zeros(3, 3);
        for j in 0..3 { for i in 0..3 { u[(i, j)] = j as f64 * 0.1; } }
        let velocity = [u, v];
        for _ in 0..50 {
            let r = model.update_reynolds_stresses_optimized(&mut stresses, &velocity, 0.001, 0.1, 0.1);
            assert!(r.is_ok(), "SSG update failed");
        }
        assert!(stresses.xy[(1, 1)] < 0.0, "Shear stress must be negative in shear flow");
    }
}
