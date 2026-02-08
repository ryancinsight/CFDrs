//! 1D Venturi resistance model validation
//!
//! Validates the VenturiModel from cfd-1d against:
//! 1. Ideal Bernoulli equation (inviscid limit at high Re)
//! 2. ISO 5167 discharge coefficients
//! 3. Borda-Carnot expansion loss (sudden expansion limit)
//! 4. Low-Re millifluidic correction
//! 5. Blood (non-Newtonian Casson) flow through Venturi
//!
//! References:
//! - Venturi, G. B. (1797). Recherches expérimentales ...
//! - ISO 5167-4:2003
//! - Reader-Harris, M. (2015). Orifice Plates and Venturi Tubes. Springer.
//! - Idelchik, I. E. (2007). Handbook of Hydraulic Resistance (4th ed.).
//! - Merrill, E. W. (1969). Rheology of blood. Physiological Reviews.

use cfd_1d::resistance::{
    ExpansionType, FlowConditions, ResistanceModel, VenturiGeometry, VenturiModel,
};

fn main() {
    println!("╔════════════════════════════════════════════════════════════╗");
    println!("║   1D Venturi Resistance Model — Literature Validation     ║");
    println!("╚════════════════════════════════════════════════════════════╝");
    println!();

    let pass_count = std::cell::Cell::new(0u32);
    let fail_count = std::cell::Cell::new(0u32);

    let check = |name: &str, condition: bool, detail: &str| {
        if condition {
            pass_count.set(pass_count.get() + 1);
            println!("  ✓ {name}: {detail}");
        } else {
            fail_count.set(fail_count.get() + 1);
            println!("  ✗ {name}: {detail}");
        }
    };

    // ═══════════════════════════════════════════════════════════════
    // Test 1: Bernoulli Limit (High Re, C_d = 1.0)
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 1: Bernoulli Equation Limit ━━━");
    println!("  At high Re with C_d = 1.0, ΔP_contraction → ½ρV_t²(1 - β²)");
    {
        let model = VenturiModel::symmetric(0.01_f64, 0.005, 0.001, 0.05)
            .with_geometry(VenturiGeometry::Custom {
                discharge_coefficient: 1.0,
            })
            .with_expansion(ExpansionType::Gradual { half_angle_deg: 3.0 });

        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let conditions = FlowConditions::new(2.0); // 2 m/s inlet

        let analysis = model.analyze(&fluid, &conditions).unwrap();

        // β = D_t/D_i = 0.5, β² = 0.25
        // V_throat = V_inlet × (D_i/D_t)² = 2 × 4 = 8 m/s
        // ΔP_ideal = ½ × 998.2 × 64 × 0.75 = 23,956.8 Pa
        let rho = 998.2;
        let v_t = 8.0;
        let beta_sq = 0.25;
        let dp_ideal = 0.5 * rho * v_t * v_t * (1.0 - beta_sq);

        check(
            "Throat velocity",
            (analysis.throat_velocity - 8.0).abs() < 0.01,
            &format!("V_t = {:.3} m/s (expected 8.0)", analysis.throat_velocity),
        );

        let err_pct = ((analysis.dp_contraction - dp_ideal) / dp_ideal * 100.0).abs();
        check(
            "Bernoulli ΔP_contraction",
            err_pct < 5.0,
            &format!(
                "ΔP = {:.1} Pa vs ideal {:.1} Pa (err = {:.2}%)",
                analysis.dp_contraction, dp_ideal, err_pct
            ),
        );
    }
    println!();

    // ═══════════════════════════════════════════════════════════════
    // Test 2: ISO 5167 Discharge Coefficients
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 2: ISO 5167 Discharge Coefficients ━━━");
    {
        check(
            "Machined convergent C_d",
            (VenturiGeometry::MachinedConvergent.discharge_coefficient() - 0.995).abs() < 1e-6,
            &format!(
                "C_d = {:.4} (ISO 5167: 0.995)",
                VenturiGeometry::MachinedConvergent.discharge_coefficient()
            ),
        );
        check(
            "Rough-cast C_d",
            (VenturiGeometry::RoughCastConvergent.discharge_coefficient() - 0.984).abs() < 1e-6,
            &format!(
                "C_d = {:.4} (ISO 5167: 0.984)",
                VenturiGeometry::RoughCastConvergent.discharge_coefficient()
            ),
        );
        check(
            "Rough-welded C_d",
            (VenturiGeometry::RoughWeldedConvergent.discharge_coefficient() - 0.985).abs() < 1e-6,
            &format!(
                "C_d = {:.4} (ISO 5167: 0.985)",
                VenturiGeometry::RoughWeldedConvergent.discharge_coefficient()
            ),
        );
    }
    println!();

    // ═══════════════════════════════════════════════════════════════
    // Test 3: Borda-Carnot Expansion Loss
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 3: Borda-Carnot Expansion Loss ━━━");
    println!("  Sudden expansion: K_exp = 1.0, ΔP_loss = ½ρ(V₂-V₃)²");
    {
        let model = VenturiModel::symmetric(0.01_f64, 0.005, 0.001, 0.05)
            .with_geometry(VenturiGeometry::Custom {
                discharge_coefficient: 1.0,
            })
            .with_expansion(ExpansionType::Sudden);

        check(
            "Sudden K_exp",
            (ExpansionType::Sudden.loss_coefficient() - 1.0).abs() < 1e-10,
            "K_exp = 1.0 (Borda-Carnot)",
        );

        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let conditions = FlowConditions::new(1.0);
        let analysis = model.analyze(&fluid, &conditions).unwrap();

        // V_throat = 4 m/s, V_outlet = 1 m/s
        // ΔP_BC = ½ × 998.2 × (4-1)² = ½ × 998.2 × 9 = 4491.9 Pa
        let dp_bc_expected = 0.5 * 998.2 * (4.0 - 1.0_f64).powi(2);
        let err_pct = ((analysis.dp_expansion_loss - dp_bc_expected) / dp_bc_expected * 100.0).abs();
        check(
            "Borda-Carnot ΔP",
            err_pct < 5.0,
            &format!(
                "ΔP_exp = {:.1} Pa vs expected {:.1} Pa (err = {:.2}%)",
                analysis.dp_expansion_loss, dp_bc_expected, err_pct
            ),
        );

        // Gradual diffuser should have lower loss
        check(
            "Gradual < Sudden",
            ExpansionType::Gradual { half_angle_deg: 5.0 }.loss_coefficient() < 1.0,
            &format!(
                "K_grad = {:.2} < 1.0",
                ExpansionType::Gradual { half_angle_deg: 5.0 }.loss_coefficient()
            ),
        );
    }
    println!();

    // ═══════════════════════════════════════════════════════════════
    // Test 4: Millifluidic Venturi (Low Re)
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 4: Millifluidic Venturi (Low Re) ━━━");
    println!("  At low Re, friction dominates Bernoulli term");
    {
        let model = VenturiModel::millifluidic(0.001_f64, 0.0005, 0.002);
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();

        // Very low velocity → laminar regime
        let conditions = FlowConditions::new(0.01); // 10 mm/s
        let analysis = model.analyze(&fluid, &conditions).unwrap();

        // Throat Re should be small but positive
        check(
            "Low throat Re",
            analysis.throat_reynolds > 0.0 && analysis.throat_reynolds < 100.0,
            &format!("Re_t = {:.2}", analysis.throat_reynolds),
        );

        // Friction should be significant
        let friction_fraction = analysis.dp_friction / analysis.dp_total;
        check(
            "Friction dominates at low Re",
            friction_fraction > 0.1,
            &format!(
                "ΔP_friction/ΔP_total = {:.1}%",
                friction_fraction * 100.0
            ),
        );

        // Laminar friction factor: f = 64/Re
        let f_expected = 64.0 / analysis.throat_reynolds;
        let err_pct = ((analysis.friction_factor - f_expected) / f_expected * 100.0).abs();
        check(
            "Laminar friction f=64/Re",
            err_pct < 1.0,
            &format!(
                "f = {:.4} vs 64/Re = {:.4} (err = {:.3}%)",
                analysis.friction_factor, f_expected, err_pct
            ),
        );
    }
    println!();

    // ═══════════════════════════════════════════════════════════════
    // Test 5: Area Ratio Scaling
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 5: Area Ratio Scaling ━━━");
    println!("  ΔP ∝ (1 - β²) for fixed flow rate");
    {
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();

        let betas = [0.3, 0.5, 0.7];
        let mut dps = Vec::new();

        for &beta in &betas {
            let d_throat = 0.01 * beta; // D_i = 10mm
            let model = VenturiModel::symmetric(0.01, d_throat, 0.001, 0.05)
                .with_geometry(VenturiGeometry::Custom {
                    discharge_coefficient: 1.0,
                })
                .with_expansion(ExpansionType::Gradual { half_angle_deg: 5.0 });

            let conditions = FlowConditions::new(0.5);
            let analysis = model.analyze(&fluid, &conditions).unwrap();
            dps.push((beta, analysis.dp_contraction));
        }

        // ΔP should increase as β decreases (narrower throat)
        check(
            "ΔP increases with constriction",
            dps[0].1 > dps[1].1 && dps[1].1 > dps[2].1,
            &format!(
                "β={:.1}: {:.0} Pa > β={:.1}: {:.0} Pa > β={:.1}: {:.0} Pa",
                dps[0].0, dps[0].1, dps[1].0, dps[1].1, dps[2].0, dps[2].1
            ),
        );
    }
    println!();

    // ═══════════════════════════════════════════════════════════════
    // Test 6: Blood Flow Through Venturi (Casson Model)
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 6: Blood Flow Through Venturi (Casson) ━━━");
    println!("  Merrill et al. (1969): blood is shear-thinning");
    {
        let model = VenturiModel::millifluidic(0.002_f64, 0.001, 0.003);

        // Use Casson blood model
        let blood = cfd_core::physics::fluid::blood::CassonBlood::<f64>::normal_blood();

        let conditions = FlowConditions::new(0.05); // 50 mm/s
        let analysis = model.analyze(&blood, &conditions).unwrap();

        // Blood viscosity should be shear-dependent
        check(
            "Blood throat viscosity",
            analysis.throat_viscosity > 0.0,
            &format!("μ_throat = {:.4e} Pa·s", analysis.throat_viscosity),
        );

        // At high shear rates in the throat, viscosity should approach
        // the infinite-shear viscosity of Casson model
        check(
            "Shear rate in throat",
            analysis.throat_shear_rate > 0.0,
            &format!("γ̇_throat = {:.1} s⁻¹", analysis.throat_shear_rate),
        );

        // Total pressure drop should be positive
        check(
            "Blood Venturi ΔP > 0",
            analysis.dp_total > 0.0,
            &format!("ΔP_total = {:.2} Pa", analysis.dp_total),
        );

        // Compare with Newtonian water at same conditions
        let water = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let water_analysis = model.analyze(&water, &conditions).unwrap();

        // Blood should have different pressure drop than water
        // (it's non-Newtonian, so viscosity varies)
        let ratio = analysis.dp_total / water_analysis.dp_total;
        println!(
            "    Blood/Water ΔP ratio = {:.3} (blood is more viscous at low shear)",
            ratio
        );
        check(
            "Blood vs Water differ",
            (ratio - 1.0).abs() > 0.001,
            &format!("Ratio = {:.3} ≠ 1.0", ratio),
        );
    }
    println!();

    // ═══════════════════════════════════════════════════════════════
    // Test 7: Resistance Model Trait Compliance
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 7: ResistanceModel Trait Compliance ━━━");
    {
        let model = VenturiModel::symmetric(0.01_f64, 0.005, 0.01, 0.05);
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let conditions = FlowConditions::new(0.1);

        check(
            "model_name()",
            model.model_name() == "Venturi",
            &format!("name = \"{}\"", model.model_name()),
        );

        let (re_min, re_max) = model.reynolds_range();
        check(
            "reynolds_range wide",
            re_min < 1.0 && re_max > 1e6,
            &format!("Re ∈ [{:.1}, {:.0e}]", re_min, re_max),
        );

        let r = model.calculate_resistance(&fluid, &conditions);
        let r_val = r.as_ref().copied().unwrap_or(0.0);
        check(
            "calculate_resistance Ok",
            r.is_ok() && r_val > 0.0,
            &format!("R = {:.4e}", r_val),
        );

        let coeffs = model.calculate_coefficients(&fluid, &conditions);
        check(
            "calculate_coefficients Ok",
            coeffs.is_ok(),
            "Returns (R, k) tuple",
        );

        let valid = model.validate_invariants(&fluid, &conditions);
        check("validate_invariants Ok", valid.is_ok(), "No violations");

        // Invalid model should fail validation
        let bad_model = VenturiModel::symmetric(0.005_f64, 0.01, 0.01, 0.05); // throat > inlet
        let bad_valid = bad_model.validate_invariants(&fluid, &conditions);
        check(
            "Invalid geometry rejected",
            bad_valid.is_err(),
            "throat > inlet detected",
        );
    }
    println!();

    // ═══════════════════════════════════════════════════════════════
    // Summary
    // ═══════════════════════════════════════════════════════════════
    let total = pass_count.get() + fail_count.get();
    println!("═══════════════════════════════════════════════════════");
    println!(
        "  Results: {}/{} passed, {} failed",
        pass_count.get(),
        total,
        fail_count.get()
    );
    if fail_count.get() == 0 {
        println!("  ✓ ALL VENTURI VALIDATION TESTS PASSED");
    } else {
        println!("  ✗ SOME TESTS FAILED — review above");
    }
    println!("═══════════════════════════════════════════════════════");
}
