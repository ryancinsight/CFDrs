//! 1D Serpentine channel resistance model validation
//!
//! Validates the SerpentineModel from cfd-1d against:
//! 1. Dean number calculation (De = Re √(D_h / 2R_c))
//! 2. White (1929) low-De curvature enhancement
//! 3. Ito (1959) moderate-De friction factor correlation
//! 4. Shah-London (1978) rectangular duct fRe product
//! 5. Idelchik (2007) bend minor loss coefficients
//! 6. More bends → more resistance (physical monotonicity)
//! 7. Blood (Casson) vs water in serpentine
//!
//! References:
//! - Dean, W. R. (1928). Proc. R. Soc. Lond. A, 121(787), 402-420.
//! - White, C. M. (1929). Proc. R. Soc. Lond. A, 123(792), 645-663.
//! - Ito, H. (1959). ASME J. Basic Eng., 81(2), 123-134.
//! - Shah, R. K. & London, A. L. (1978). Laminar Flow Forced Convection in Ducts.
//! - Idelchik, I. E. (2007). Handbook of Hydraulic Resistance (4th ed.).
//! - Squires, T. M. & Quake, S. R. (2005). Rev. Mod. Phys. 77, 977.
//! - Stroock, A. D. et al. (2002). Science 295, 647.

use cfd_1d::resistance::{
    BendType, FlowConditions, ResistanceModel, SerpentineCrossSection,
    SerpentineModel,
};

fn main() {
    println!("╔════════════════════════════════════════════════════════════╗");
    println!("║  1D Serpentine Resistance Model — Literature Validation   ║");
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
    // Test 1: Dean Number Calculation
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 1: Dean Number (Dean 1928) ━━━");
    println!("  De = Re × √(D_h / (2 R_c))");
    {
        // D_h = 1 mm, R_c = 5 mm
        // De = Re × √(0.001 / 0.01) = Re × 0.3162
        let model = SerpentineModel {
            straight_length: 0.02_f64,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };

        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();

        // At Re ≈ 100, De ≈ 31.6
        let conditions = FlowConditions::new(0.1); // velocity m/s
        let analysis = model.analyze(&fluid, &conditions).unwrap();

        let expected_de = analysis.reynolds * (0.001_f64 / (2.0 * 0.005)).sqrt();
        let err_pct = ((analysis.dean_number - expected_de) / expected_de * 100.0).abs();
        check(
            "Dean number formula",
            err_pct < 0.01,
            &format!(
                "De = {:.4} vs expected {:.4} (err = {:.4}%)",
                analysis.dean_number, expected_de, err_pct
            ),
        );

        // Dean number should scale linearly with Re
        let conditions_2x = FlowConditions::new(0.2);
        let analysis_2x = model.analyze(&fluid, &conditions_2x).unwrap();
        let ratio = analysis_2x.dean_number / analysis.dean_number;
        let re_ratio = analysis_2x.reynolds / analysis.reynolds;
        check(
            "De ∝ Re linearity",
            ((ratio - re_ratio) / re_ratio).abs() < 0.05,
            &format!(
                "De ratio = {:.3}, Re ratio = {:.3}",
                ratio, re_ratio
            ),
        );
    }
    println!();

    // ═══════════════════════════════════════════════════════════════
    // Test 2: Curvature Enhancement (White 1929, Ito 1959)
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 2: Curvature Enhancement (Ito 1959) ━━━");
    {
        let model = SerpentineModel {
            straight_length: 0.02_f64,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();

        // At very low De (< 11.6), enhancement ≈ 1.0
        let conditions_slow = FlowConditions::new(0.001);
        let analysis_slow = model.analyze(&fluid, &conditions_slow).unwrap();
        check(
            "Low De: enhancement ≈ 1",
            (analysis_slow.curvature_enhancement - 1.0).abs() < 0.1,
            &format!(
                "De = {:.3}, f_c/f_s = {:.4}",
                analysis_slow.dean_number, analysis_slow.curvature_enhancement
            ),
        );

        // At moderate De (> 11.6), enhancement > 1.0
        let conditions_fast = FlowConditions::new(1.0);
        let analysis_fast = model.analyze(&fluid, &conditions_fast).unwrap();
        check(
            "Moderate De: enhancement > 1",
            analysis_fast.curvature_enhancement > 1.0,
            &format!(
                "De = {:.1}, f_c/f_s = {:.4}",
                analysis_fast.dean_number, analysis_fast.curvature_enhancement
            ),
        );

        // Enhancement should increase monotonically with De
        check(
            "Enhancement monotonic",
            analysis_fast.curvature_enhancement > analysis_slow.curvature_enhancement,
            &format!(
                "{:.4} > {:.4}",
                analysis_fast.curvature_enhancement, analysis_slow.curvature_enhancement
            ),
        );

        // f_curved should be higher than f_straight
        check(
            "f_curved > f_straight",
            analysis_fast.friction_factor_curved > analysis_fast.friction_factor_straight,
            &format!(
                "f_c = {:.6} > f_s = {:.6}",
                analysis_fast.friction_factor_curved, analysis_fast.friction_factor_straight
            ),
        );
    }
    println!();

    // ═══════════════════════════════════════════════════════════════
    // Test 3: Shah-London Rectangular Correction
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 3: Shah-London Rectangular Duct (1978) ━━━");
    println!("  f·Re = 96(1 - 1.3553α + 1.9467α² - ...) for rect ducts");
    {
        // Square duct: α = 1, fRe = 56.92 → fRe/64 ≈ 0.889
        let sq = SerpentineCrossSection::Rectangular {
            width: 0.001,
            height: 0.001,
        };
        let factor_sq = sq.shah_london_fre_factor();
        let expected_sq = 96.0
            * (1.0 - 1.3553 + 1.9467 - 1.7012 + 0.9564 - 0.2537)
            / 64.0;
        check(
            "Square α=1 factor",
            (factor_sq - expected_sq).abs() < 0.001,
            &format!("C = {:.4} vs {:.4}", factor_sq, expected_sq),
        );

        // High aspect ratio: α ≈ 0.1 (10:1), fRe → 96 × ≈0.893 / 64
        let rect = SerpentineCrossSection::Rectangular {
            width: 0.010,
            height: 0.001,
        };
        let factor_rect = rect.shah_london_fre_factor();
        check(
            "10:1 aspect factor > 1",
            factor_rect > 1.0,
            &format!("C(α=0.1) = {:.4} (fRe > 64)", factor_rect),
        );

        // Circular should have factor = 1.0
        let circ = SerpentineCrossSection::Circular { diameter: 0.001 };
        check(
            "Circular factor = 1.0",
            (circ.shah_london_fre_factor() - 1.0).abs() < 1e-10,
            &format!("C_circ = {:.6}", circ.shah_london_fre_factor()),
        );

        // Rectangular serpentine vs circular at same D_h
        let model_circ = SerpentineModel {
            straight_length: 0.02_f64,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };

        // Square with same D_h as circular (D_h = 2wh/(w+h) = w for square)
        let model_rect = SerpentineModel {
            straight_length: 0.02_f64,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Rectangular {
                width: 0.001,
                height: 0.001,
            },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };

        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let conditions = FlowConditions::new(0.1);

        let dp_circ = model_circ.analyze(&fluid, &conditions).unwrap().dp_total;
        let dp_rect = model_rect.analyze(&fluid, &conditions).unwrap().dp_total;

        // Square duct has lower fRe (56.9 vs 64), so slightly lower friction
        // but different area → different velocity → different ΔP
        println!(
            "    ΔP_circ = {:.2} Pa, ΔP_rect = {:.2} Pa (ratio = {:.3})",
            dp_circ,
            dp_rect,
            dp_rect / dp_circ
        );
    }
    println!();

    // ═══════════════════════════════════════════════════════════════
    // Test 4: Bend Minor Loss Coefficients (Idelchik 2007)
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 4: Bend Minor Loss Coefficients (Idelchik 2007) ━━━");
    {
        // Sharp 180° bend: K = 2.2 + 250/Re
        let k_sharp_100 = BendType::Sharp.loss_coefficient(100.0_f64);
        let expected_sharp = 2.2 + 250.0 / 100.0;
        check(
            "Sharp K(Re=100)",
            (k_sharp_100 - expected_sharp).abs() < 1e-6,
            &format!("K = {:.4} (expected {:.1})", k_sharp_100, expected_sharp),
        );

        // Smooth (R/D=5): K = 0.3 + 75/Re
        let k_smooth_1000 =
            BendType::Smooth { radius_to_dh_ratio: 5.0 }.loss_coefficient(1000.0_f64);
        let expected_smooth = 0.3 + 75.0 / 1000.0;
        check(
            "Smooth K(Re=1000)",
            (k_smooth_1000 - expected_smooth).abs() < 1e-6,
            &format!("K = {:.4} (expected {:.4})", k_smooth_1000, expected_smooth),
        );

        // Sharp > smooth at any Re
        let k_sharp_500 = BendType::Sharp.loss_coefficient(500.0_f64);
        let k_smooth_500 =
            BendType::Smooth { radius_to_dh_ratio: 3.0 }.loss_coefficient(500.0_f64);
        check(
            "Sharp > Smooth",
            k_sharp_500 > k_smooth_500,
            &format!(
                "K_sharp = {:.3} > K_smooth = {:.3} at Re=500",
                k_sharp_500, k_smooth_500
            ),
        );

        // At high Re, K approaches the constant C₁
        let k_sharp_inf = BendType::Sharp.loss_coefficient(1e6_f64);
        check(
            "High Re → C₁",
            (k_sharp_inf - 2.2).abs() < 0.001,
            &format!("K(Re=10⁶) = {:.4} ≈ C₁ = 2.2", k_sharp_inf),
        );
    }
    println!();

    // ═══════════════════════════════════════════════════════════════
    // Test 5: More Bends → More Resistance
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 5: Monotonicity — More Bends → More Resistance ━━━");
    {
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let conditions = FlowConditions::new(0.1);

        let mut resistances = Vec::new();
        for n in [2, 5, 10, 20] {
            let model = SerpentineModel {
                straight_length: 0.02_f64,
                num_segments: n,
                cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
                bend_radius: 0.005,
                bend_type: BendType::Smooth {
                    radius_to_dh_ratio: 5.0,
                },
            };
            let r = model.calculate_resistance(&fluid, &conditions).unwrap();
            resistances.push((n, r));
        }

        // Each successive should be larger
        let mut monotonic = true;
        for i in 1..resistances.len() {
            if resistances[i].1 <= resistances[i - 1].1 {
                monotonic = false;
                break;
            }
        }
        check(
            "R monotonically increasing",
            monotonic,
            &format!(
                "R(2)={:.2e}, R(5)={:.2e}, R(10)={:.2e}, R(20)={:.2e}",
                resistances[0].1, resistances[1].1, resistances[2].1, resistances[3].1
            ),
        );

        // The incremental resistance per bend should be roughly constant
        // (at fixed Re, K_bend is the same)
        let dr_per_bend: Vec<f64> = (1..resistances.len())
            .map(|i| {
                let dn = (resistances[i].0 - resistances[i - 1].0) as f64;
                (resistances[i].1 - resistances[i - 1].1) / dn
            })
            .collect();

        let mean_dr = dr_per_bend.iter().sum::<f64>() / dr_per_bend.len() as f64;
        let max_dev = dr_per_bend
            .iter()
            .map(|x| ((x - mean_dr) / mean_dr).abs())
            .fold(0.0_f64, f64::max);

        check(
            "ΔR/bend ≈ constant",
            max_dev < 0.5,
            &format!(
                "Mean ΔR/bend = {:.4e}, max deviation = {:.1}%",
                mean_dr,
                max_dev * 100.0
            ),
        );
    }
    println!();

    // ═══════════════════════════════════════════════════════════════
    // Test 6: Laminar Friction Factor f = C/Re
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 6: Laminar Base Friction (Hagen-Poiseuille) ━━━");
    {
        let model = SerpentineModel {
            straight_length: 0.02_f64,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.01, // Large R_c → low curvature
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 10.0,
            },
        };
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();

        // At very low velocity (Re ~ 1), f ≈ 64/Re for circular
        let conditions = FlowConditions::new(0.001); // 1 mm/s
        let analysis = model.analyze(&fluid, &conditions).unwrap();

        let f_expected = 64.0 / analysis.reynolds;
        let err_pct = ((analysis.friction_factor_straight - f_expected) / f_expected * 100.0).abs();
        check(
            "f_straight = 64/Re",
            err_pct < 1.0,
            &format!(
                "f = {:.4} vs 64/Re = {:.4} (Re = {:.3}, err = {:.3}%)",
                analysis.friction_factor_straight, f_expected, analysis.reynolds, err_pct
            ),
        );
    }
    println!();

    // ═══════════════════════════════════════════════════════════════
    // Test 7: Blood Flow in Serpentine (Casson Model)
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 7: Blood vs Water in Serpentine ━━━");
    println!("  Merrill (1969): blood is shear-thinning (Casson fluid)");
    {
        let model = SerpentineModel {
            straight_length: 0.03_f64,
            num_segments: 8,
            cross_section: SerpentineCrossSection::Rectangular {
                width: 0.0005,
                height: 0.0003,
            },
            bend_radius: 0.002,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };

        let water = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let blood = cfd_core::physics::fluid::blood::CassonBlood::<f64>::normal_blood();

        let conditions = FlowConditions::new(0.05); // 50 mm/s

        let water_analysis = model.analyze(&water, &conditions).unwrap();
        let blood_analysis = model.analyze(&blood, &conditions).unwrap();

        // Blood viscosity should be higher than water
        check(
            "Blood μ > water μ",
            blood_analysis.wall_viscosity > water_analysis.wall_viscosity,
            &format!(
                "μ_blood = {:.4e} > μ_water = {:.4e} Pa·s",
                blood_analysis.wall_viscosity, water_analysis.wall_viscosity
            ),
        );

        // Blood ΔP should be higher than water (higher effective viscosity)
        check(
            "Blood ΔP > water ΔP",
            blood_analysis.dp_total > water_analysis.dp_total,
            &format!(
                "ΔP_blood = {:.2} Pa > ΔP_water = {:.2} Pa",
                blood_analysis.dp_total, water_analysis.dp_total
            ),
        );

        let viscosity_ratio = blood_analysis.wall_viscosity / water_analysis.wall_viscosity;
        println!(
            "    Viscosity ratio μ_blood/μ_water = {:.2}",
            viscosity_ratio
        );
        println!(
            "    ΔP ratio = {:.2}",
            blood_analysis.dp_total / water_analysis.dp_total
        );

        // Wall shear rate should be positive
        check(
            "Wall shear rate > 0",
            blood_analysis.wall_shear_rate > 0.0,
            &format!("γ̇ = {:.1} s⁻¹", blood_analysis.wall_shear_rate),
        );
    }
    println!();

    // ═══════════════════════════════════════════════════════════════
    // Test 8: ResistanceModel Trait Compliance
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 8: ResistanceModel Trait Compliance ━━━");
    {
        let model = SerpentineModel {
            straight_length: 0.02_f64,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let conditions = FlowConditions::new(0.1);

        check(
            "model_name()",
            model.model_name() == "Serpentine",
            &format!("name = \"{}\"", model.model_name()),
        );

        let (re_min, re_max) = model.reynolds_range();
        check(
            "reynolds_range wide",
            re_min < 1.0 && re_max > 1e4,
            &format!("Re ∈ [{:.2}, {:.0e}]", re_min, re_max),
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

        // Invalid model should fail
        let bad_model = SerpentineModel {
            straight_length: -0.01_f64,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };
        let bad_valid = bad_model.validate_invariants(&fluid, &conditions);
        check(
            "Negative length rejected",
            bad_valid.is_err(),
            "straight_length < 0 detected",
        );
    }
    println!();

    // ═══════════════════════════════════════════════════════════════
    // Test 9: Millifluidic Design Space (Squires & Quake 2005)
    // ═══════════════════════════════════════════════════════════════
    println!("━━━ Test 9: Millifluidic Design Space (Squires & Quake 2005) ━━━");
    println!("  Typical millifluidic: Re ~ 0.1-100, De ~ 0.01-10");
    {
        let model = SerpentineModel::millifluidic_rectangular(
            0.0005, // 500 μm wide
            0.0003, // 300 μm deep
            0.005_f64,
            10,
            0.002,
        );

        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let conditions = FlowConditions::new(0.01); // 10 mm/s

        let analysis = model.analyze(&fluid, &conditions).unwrap();

        check(
            "Millifluidic Re < 100",
            analysis.reynolds < 100.0,
            &format!("Re = {:.3}", analysis.reynolds),
        );

        check(
            "Millifluidic De < Re",
            analysis.dean_number < analysis.reynolds,
            &format!(
                "De = {:.3} < Re = {:.3}",
                analysis.dean_number, analysis.reynolds
            ),
        );

        check(
            "Friction dominates in millifluidic",
            analysis.dp_friction > analysis.dp_bends,
            &format!(
                "ΔP_friction = {:.3} Pa > ΔP_bends = {:.3} Pa",
                analysis.dp_friction, analysis.dp_bends
            ),
        );

        println!(
            "    D_h = {:.0} μm, # bends = {}, ΔP_total = {:.3} Pa",
            model.cross_section.hydraulic_diameter() * 1e6,
            analysis.num_bends,
            analysis.dp_total
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
        println!("  ✓ ALL SERPENTINE VALIDATION TESTS PASSED");
    } else {
        println!("  ✗ SOME TESTS FAILED — review above");
    }
    println!("═══════════════════════════════════════════════════════");
}
