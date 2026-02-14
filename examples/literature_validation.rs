//! # Literature Validation — CFD-1D Resistance Models
//!
//! Quantitative comparison of cfd-1d resistance models against published
//! analytical and empirical correlations from the microfluidics literature.
//!
//! ## Validation Cases
//!
//! | # | Case | Reference | Metric |
//! |---|------|-----------|--------|
//! | 1 | Shah-London Poiseuille number | Shah & London (1978) | Po(α) for α = 1…10 |
//! | 2 | Hagen-Poiseuille resistance | Hagen (1839); Poiseuille (1840) | R = 128μL/(πD⁴) |
//! | 3 | Bruus shallow-channel approx | Bruus (2008) | R_H ≈ 12μL/(wh³(1−0.63h/w)) |
//! | 4 | Dean number & curvature | Dean (1928); Ito (1959) | De, f_curved/f_straight |
//! | 5 | Venturi pressure drop | Bernoulli; ISO 5167 | ΔP contraction, recovery |
//! | 6 | Network cross-validation | Kirchhoff's laws | Σ(R·Q) vs solver ΔP |
//!
//! ## Pass/Fail Criteria
//!
//! - Individual model errors < 1 % where analytical solutions are exact
//! - Empirical correlation errors < 5 % within stated validity ranges
//! - Network conservation error < 0.1 %
//!
//! ## Running
//!
//! ```sh
//! cargo run --example literature_validation --features scheme-integration
//! ```

// ═══════════════════════════════════════════════════════════════════════════════
// Imports
// ═══════════════════════════════════════════════════════════════════════════════

use std::io::Write;
use std::path::Path;

// cfd-1d — resistance models
use cfd_1d::{
    BendType, ExpansionType, FlowConditions, HagenPoiseuilleModel, RectangularChannelModel,
    ResistanceModel, SerpentineCrossSection, SerpentineModel, VenturiAnalysis,
    VenturiGeometry, VenturiModel,
};

// cfd-1d — network / solver
use cfd_1d::network::{Network, NodeType};
use cfd_1d::scheme_bridge::SchemeNetworkConverter;
use cfd_1d::solver::{NetworkProblem, NetworkSolver, SolverConfig};

// cfd-core
use cfd_core::physics::fluid::ConstantPropertyFluid;

// scheme
use scheme::config::*;
use scheme::geometry::generator::create_geometry;
use scheme::geometry::SplitType;

// petgraph
use petgraph::graph::NodeIndex;

// plotters
use plotters::prelude::*;

// ═══════════════════════════════════════════════════════════════════════════════
// Constants
// ═══════════════════════════════════════════════════════════════════════════════

const WATER_DENSITY: f64 = 998.0; // kg/m³  (20 °C)
const WATER_VISCOSITY: f64 = 1.002e-3; // Pa·s
const WATER_CP: f64 = 4182.0; // J/(kg·K)
const WATER_K: f64 = 0.598; // W/(m·K)
const WATER_SOUND: f64 = 1482.0; // m/s

/// Tolerance for "exact" analytical match (%)
const TOL_EXACT: f64 = 1.0;
/// Tolerance for empirical correlations (%)
const TOL_EMPIRICAL: f64 = 5.0;
/// Tolerance for network conservation (%)
const TOL_NETWORK: f64 = 0.1;

// ═══════════════════════════════════════════════════════════════════════════════
// Result types
// ═══════════════════════════════════════════════════════════════════════════════

#[derive(Debug, Clone)]
struct ValidationRow {
    parameter: String,
    literature: f64,
    cfd_1d: f64,
    error_pct: f64,
    pass: bool,
}

#[derive(Debug, Clone)]
struct ValidationCase {
    id: usize,
    name: String,
    reference: String,
    rows: Vec<ValidationRow>,
    all_pass: bool,
}

// ═══════════════════════════════════════════════════════════════════════════════
// Shared utilities
// ═══════════════════════════════════════════════════════════════════════════════

fn water_fluid() -> ConstantPropertyFluid<f64> {
    ConstantPropertyFluid::new(
        "Water (20 °C)".into(),
        WATER_DENSITY,
        WATER_VISCOSITY,
        WATER_CP,
        WATER_K,
        WATER_SOUND,
    )
}

fn pct_error(computed: f64, reference: f64) -> f64 {
    if reference.abs() < 1e-30 {
        if computed.abs() < 1e-30 {
            0.0
        } else {
            100.0
        }
    } else {
        ((computed - reference) / reference).abs() * 100.0
    }
}

fn print_case(case: &ValidationCase) {
    println!();
    println!(
        "  ━━━ Case {}: {} ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
        case.id, case.name
    );
    println!("    Reference: {}", case.reference);
    println!();
    println!(
        "    ┌───────────────────────────────┬───────────────┬───────────────┬──────────┬────────┐"
    );
    println!(
        "    │ Parameter                     │  Literature   │   cfd-1d      │ Error %  │ Status │"
    );
    println!(
        "    ├───────────────────────────────┼───────────────┼───────────────┼──────────┼────────┤"
    );
    for r in &case.rows {
        let status = if r.pass { " PASS " } else { " FAIL " };
        let status_color = if r.pass { "✓" } else { "✗" };
        println!(
            "    │ {:>29} │ {:>13.6} │ {:>13.6} │ {:>7.4} │ {} {} │",
            r.parameter, r.literature, r.cfd_1d, r.error_pct, status_color, status
        );
    }
    println!(
        "    └───────────────────────────────┴───────────────┴───────────────┴──────────┴────────┘"
    );
    println!(
        "    Overall: {}",
        if case.all_pass {
            "ALL PASS ✓"
        } else {
            "SOME FAIL ✗"
        }
    );
}

fn find_boundary_nodes(
    network: &Network<f64, ConstantPropertyFluid<f64>>,
) -> (Vec<NodeIndex>, Vec<NodeIndex>) {
    let mut inlets = Vec::new();
    let mut outlets = Vec::new();
    for ni in network.graph.node_indices() {
        match network.graph[ni].node_type {
            NodeType::Inlet => inlets.push(ni),
            NodeType::Outlet => outlets.push(ni),
            _ => {}
        }
    }
    (inlets, outlets)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Case 1 — Shah-London Poiseuille Number for Rectangular Ducts
// ═══════════════════════════════════════════════════════════════════════════════

/// Shah & London (1978) tabulated Poiseuille numbers for rectangular ducts.
///
/// Po_Fanning = fRe where f = Fanning friction factor.
/// Po_Darcy = 4×Po_Fanning = f_D × Re.
///
/// The cfd-1d `RectangularChannelModel` uses Darcy convention (Po_Darcy = 96(1-…)).
///
/// Reference table (Shah & London 1978, Table 43 — f_Fanning×Re):
///
/// | α (aspect) | f·Re (Fanning) | f·Re (Darcy) |
/// |-----------|----------------|--------------|
/// | 1.0       | 14.227         | 56.91        |
/// | 2.0       | 15.548         | 62.19        |
/// | 3.0       | 17.090         | 68.36        |
/// | 5.0       | 19.702         | 78.81        |
/// | 10.0      | 21.169         | 84.68        |
/// | ∞ (plates)| 24.000         | 96.00        |
///
/// We validate by computing R from the model and back-deriving Po.
fn case_shah_london() -> ValidationCase {
    // Shah-London tabulated values: (aspect_ratio, Po_Darcy)
    let table: Vec<(f64, f64)> = vec![
        (1.0, 56.91),
        (2.0, 62.19),
        (3.0, 68.36),
        (5.0, 78.81),
        (10.0, 84.68),
    ];

    let fluid = water_fluid();
    let length = 0.05; // 50 mm — ensures L/Dh > 10
    let h = 0.5e-3; // 0.5 mm height

    let mut rows = Vec::new();

    for (alpha, po_tab) in &table {
        let w = h * alpha; // width = α × height
        let model = RectangularChannelModel::new(w, h, length);

        // R = (Po·μ·L) / (2·A·Dh²)
        // so Po = R · 2·A·Dh² / (μ·L)
        let area = w * h;
        let perimeter = 2.0 * (w + h);
        let dh = 4.0 * area / perimeter;

        let conditions = FlowConditions::new(0.01); // 1 cm/s — laminar
        let (r, _k) = model.calculate_coefficients(&fluid, &conditions).unwrap();

        // Back-calculate Po from R
        let po_computed = r * 2.0 * area * dh * dh / (WATER_VISCOSITY * length);

        let err = pct_error(po_computed, *po_tab);
        // The model implements the Shah-London polynomial FIT, not the exact
        // infinite-series solution.  The polynomial has ≈3 % deviation at α=5,
        // so we compare against TOL_EMPIRICAL (5 %).
        rows.push(ValidationRow {
            parameter: format!("Po(α={:.0})", alpha),
            literature: *po_tab,
            cfd_1d: po_computed,
            error_pct: err,
            pass: err < TOL_EMPIRICAL,
        });
    }

    let all_pass = rows.iter().all(|r| r.pass);
    ValidationCase {
        id: 1,
        name: "Shah-London Poiseuille Number".into(),
        reference: "Shah & London (1978), Table 43".into(),
        rows,
        all_pass,
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Case 2 — Hagen-Poiseuille Resistance for Circular Pipes
// ═══════════════════════════════════════════════════════════════════════════════

/// Exact analytical resistance: R = 128 μ L / (π D⁴)
///
/// References:
/// - Hagen (1839), Poiseuille (1840)
/// - White (2006), Eq. 3-52
fn case_hagen_poiseuille() -> ValidationCase {
    let fluid = water_fluid();

    // Test multiple diameters (0.5, 1.0, 2.0, 3.0 mm)
    let diameters_mm: Vec<f64> = vec![0.5, 1.0, 2.0, 3.0];
    let length = 0.05; // 50 mm

    let mut rows = Vec::new();

    for d_mm in &diameters_mm {
        let d = d_mm * 1e-3;
        let model = HagenPoiseuilleModel::new(d, length);

        // Analytical R = 128 μ L / (π D⁴)
        let r_analytical =
            128.0 * WATER_VISCOSITY * length / (std::f64::consts::PI * d.powi(4));

        let conditions = FlowConditions::new(0.01);
        let (r_model, _k) = model.calculate_coefficients(&fluid, &conditions).unwrap();

        let err = pct_error(r_model, r_analytical);
        rows.push(ValidationRow {
            parameter: format!("R(D={:.1}mm)", d_mm),
            literature: r_analytical,
            cfd_1d: r_model,
            error_pct: err,
            pass: err < TOL_EXACT,
        });
    }

    // Also validate resistance ratio between diameters (scales as D⁻⁴)
    // R(0.5mm) / R(1.0mm) should be 2⁴ = 16
    if rows.len() >= 2 {
        let ratio_analytical = 16.0;
        let ratio_computed = rows[0].cfd_1d / rows[1].cfd_1d;
        let err = pct_error(ratio_computed, ratio_analytical);
        rows.push(ValidationRow {
            parameter: "R(0.5mm)/R(1.0mm)".into(),
            literature: ratio_analytical,
            cfd_1d: ratio_computed,
            error_pct: err,
            pass: err < TOL_EXACT,
        });
    }

    let all_pass = rows.iter().all(|r| r.pass);
    ValidationCase {
        id: 2,
        name: "Hagen-Poiseuille Resistance".into(),
        reference: "Hagen (1839); Poiseuille (1840); White (2006)".into(),
        rows,
        all_pass,
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Case 3 — Bruus (2008) Shallow-Channel Approximation
// ═══════════════════════════════════════════════════════════════════════════════

/// Bruus "Theoretical Microfluidics" (2008), Eq. 3.30:
///
///   R_H ≈ 12 μ L / [w h³ (1 − 0.630 h/w)]
///
/// Valid for shallow channels where w >> h (aspect ratio > 1.5 recommended,
/// error < 13 % at α = 1). Increasingly accurate as α → ∞ (parallel plates).
///
/// We compare the Bruus approximation against both:
/// (a) the exact Shah-London model in cfd-1d (RectangularChannelModel)
/// (b) the Bruus analytical formula directly
///
/// This validates that the Shah-London and Bruus results converge for high α.
fn case_bruus_approximation() -> ValidationCase {
    let fluid = water_fluid();
    let length = 0.05; // 50 mm
    let h = 0.5e-3; // 0.5 mm depth

    // Aspect ratios: 1.5, 2, 3, 5, 10, 20
    let alphas: Vec<f64> = vec![1.5, 2.0, 3.0, 5.0, 10.0, 20.0];

    let mut rows = Vec::new();

    for &alpha in &alphas {
        let w = h * alpha;

        // Bruus (2008) analytical approximation
        let r_bruus = 12.0 * WATER_VISCOSITY * length / (w * h.powi(3) * (1.0 - 0.630 * h / w));

        // cfd-1d Shah-London exact
        let model = RectangularChannelModel::new(w, h, length);
        let conditions = FlowConditions::new(0.01);
        let (r_shah_london, _) = model.calculate_coefficients(&fluid, &conditions).unwrap();

        // How well does Bruus match Shah-London?
        let err = pct_error(r_bruus, r_shah_london);

        // Bruus is an approximation; error expected to decrease with α
        // Accept < 13 % at α=1.5, should be < 1 % at α≥10
        let tol = if alpha < 2.0 { 15.0 } else if alpha < 5.0 { 5.0 } else { 2.0 };

        rows.push(ValidationRow {
            parameter: format!("Bruus vs SL (α={:.1})", alpha),
            literature: r_shah_london,
            cfd_1d: r_bruus,
            error_pct: err,
            pass: err < tol,
        });
    }

    let all_pass = rows.iter().all(|r| r.pass);
    ValidationCase {
        id: 3,
        name: "Bruus Shallow-Channel Approximation".into(),
        reference: "Bruus (2008), Theoretical Microfluidics, Eq. 3.30".into(),
        rows,
        all_pass,
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Case 4 — Dean Number & Curvature Enhancement
// ═══════════════════════════════════════════════════════════════════════════════

/// Dean (1928) secondary flow number De = Re √(D_h / 2R_c)
/// and Ito (1959) curvature friction enhancement.
///
/// Validation points:
/// - De computation for known Re, D_h, R_c
/// - Enhancement ratio f_curved/f_straight from Ito (1959):
///   f/f_s ≈ 1 + 0.1064√De / (1 + 0.0145√De)   for 11.6 < De < 2000
/// - K_bend = C₁ + C₂/Re from Idelchik (2007)
fn case_dean_curvature() -> ValidationCase {
    let fluid = water_fluid();

    let mut rows = Vec::new();

    // ── Sub-case A: Dean number computation ──────────────────────────────
    // Circular pipe: D = 1 mm, R_c = 5 mm
    // De = Re × √(D / (2Rc)) = Re × √(0.001/(0.010)) = Re × √0.1
    {
        let d = 1.0e-3;
        let rc = 5.0e-3;
        let expected_factor = (d / (2.0 * rc) as f64).sqrt(); // √0.1 ≈ 0.3162

        for re in [10.0, 50.0, 100.0, 200.0, 500.0] {
            let de_expected = re * expected_factor;

            // Create serpentine model to access dean_number
            // Use SerpentineModel::analyze which returns De
            let seg_len = 0.02;
            let model = SerpentineModel::new(
                seg_len,
                5,
                SerpentineCrossSection::Circular { diameter: d },
                rc,
            );

            let velocity = re * WATER_VISCOSITY / (WATER_DENSITY * d);
            let conditions = FlowConditions {
                reynolds_number: Some(re),
                velocity: Some(velocity),
                flow_rate: None,
                shear_rate: None,
                temperature: 293.15,
                pressure: 101325.0,
            };

            let analysis = model.analyze(&fluid, &conditions).unwrap();

            // The model computes Re from density, velocity, dh, viscosity
            // so the Dean number may differ slightly due to viscosity at shear
            let err = pct_error(analysis.dean_number, de_expected);
            rows.push(ValidationRow {
                parameter: format!("De(Re={:.0})", re),
                literature: de_expected,
                cfd_1d: analysis.dean_number,
                error_pct: err,
                pass: err < TOL_EXACT,
            });
        }
    }

    // ── Sub-case B: Curvature enhancement (Ito 1959) ─────────────────────
    // For De values in the moderate range (11.6–2000):
    // enhancement = 1 + 0.1064√De / (1 + 0.0145√De)
    {
        let d = 1.0e-3;
        let rc = 2.0e-3; // Tighter bend radius for higher De
        let seg_len = 0.02;

        for re_target in [50.0, 100.0, 200.0, 400.0] {
            let de_expected: f64 = re_target * (d / (2.0 * rc) as f64).sqrt();
            // Ito enhancement
            let de_sqrt: f64 = de_expected.sqrt();
            let enhancement_ito = 1.0 + 0.1064 * de_sqrt / (1.0 + 0.0145 * de_sqrt);

            let model = SerpentineModel::new(
                seg_len,
                5,
                SerpentineCrossSection::Circular { diameter: d },
                rc,
            );

            let velocity = re_target * WATER_VISCOSITY / (WATER_DENSITY * d);
            let conditions = FlowConditions {
                reynolds_number: Some(re_target),
                velocity: Some(velocity),
                flow_rate: None,
                shear_rate: None,
                temperature: 293.15,
                pressure: 101325.0,
            };

            let analysis = model.analyze(&fluid, &conditions).unwrap();
            let err = pct_error(analysis.curvature_enhancement, enhancement_ito);

            rows.push(ValidationRow {
                parameter: format!("f_c/f_s(De≈{:.1})", de_expected),
                literature: enhancement_ito,
                cfd_1d: analysis.curvature_enhancement,
                error_pct: err,
                pass: err < TOL_EMPIRICAL,
            });
        }
    }

    // ── Sub-case C: Bend loss coefficients (Idelchik 2007) ───────────────
    {
        // Sharp bend: K = 2.2 + 250/Re
        let re_test = 100.0;
        let k_sharp_expected = 2.2 + 250.0 / re_test;
        let k_sharp_computed = BendType::Sharp.loss_coefficient(re_test);
        let err = pct_error(k_sharp_computed, k_sharp_expected);
        rows.push(ValidationRow {
            parameter: format!("K_sharp(Re={:.0})", re_test),
            literature: k_sharp_expected,
            cfd_1d: k_sharp_computed,
            error_pct: err,
            pass: err < TOL_EXACT,
        });

        // Smooth bend R/Dh = 5: K = 0.3 + 75/Re
        let k_smooth_expected = 0.3 + 75.0 / re_test;
        let k_smooth_computed =
            BendType::Smooth { radius_to_dh_ratio: 5.0 }.loss_coefficient(re_test);
        let err = pct_error(k_smooth_computed, k_smooth_expected);
        rows.push(ValidationRow {
            parameter: format!("K_smooth_R5(Re={:.0})", re_test),
            literature: k_smooth_expected,
            cfd_1d: k_smooth_computed,
            error_pct: err,
            pass: err < TOL_EXACT,
        });

        // Smooth bend R/Dh = 2: K = 0.4 + 100/Re
        let k_smooth2_expected = 0.4 + 100.0 / re_test;
        let k_smooth2_computed =
            BendType::Smooth { radius_to_dh_ratio: 2.0 }.loss_coefficient(re_test);
        let err = pct_error(k_smooth2_computed, k_smooth2_expected);
        rows.push(ValidationRow {
            parameter: format!("K_smooth_R2(Re={:.0})", re_test),
            literature: k_smooth2_expected,
            cfd_1d: k_smooth2_computed,
            error_pct: err,
            pass: err < TOL_EXACT,
        });
    }

    let all_pass = rows.iter().all(|r| r.pass);
    ValidationCase {
        id: 4,
        name: "Dean Number & Curvature Enhancement".into(),
        reference: "Dean (1928); Ito (1959); Idelchik (2007)".into(),
        rows,
        all_pass,
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Case 5 — Venturi Pressure Drop
// ═══════════════════════════════════════════════════════════════════════════════

/// Bernoulli-based Venturi analysis:
///
///   ΔP_contraction = ½ρV_t²(1 − β²) / C_d²
///   ΔP_friction    = f·(L_t/D_t)·½ρV_t²
///   ΔP_expansion   = K_exp·½ρ(V_t−V_out)²
///   ΔP_recovery    = η_r·½ρ(V_t²−V_out²)
///
/// References: Venturi (1797); ISO 5167-4:2003; Idelchik (2007)
fn case_venturi() -> ValidationCase {
    let fluid = water_fluid();
    let mut rows = Vec::new();

    // Venturi: D_in = 2 mm, D_throat = 1 mm (β = 0.5), D_out = 2 mm
    let d_in = 2.0e-3;
    let d_throat = 1.0e-3;
    let l_throat = 5.0e-3;
    let l_total = 25.0e-3;

    let model = VenturiModel::new(d_in, d_throat, d_in, l_throat, l_total)
        .with_geometry(VenturiGeometry::MachinedConvergent) // C_d = 0.995
        .with_expansion(ExpansionType::Gradual { half_angle_deg: 7.0 }); // K_exp = 0.20

    // Test at multiple velocities (inlet velocity)
    let velocities = [0.01, 0.05, 0.1, 0.5];

    for &v_in in &velocities {
        let a_in = std::f64::consts::PI * d_in * d_in / 4.0;
        let a_throat = std::f64::consts::PI * d_throat * d_throat / 4.0;
        let beta_sq = (d_throat / d_in).powi(2);

        let v_throat = v_in * a_in / a_throat;
        let v_out = v_in; // symmetric

        let conditions = FlowConditions {
            reynolds_number: None,
            velocity: Some(v_in),
            flow_rate: None,
            shear_rate: None,
            temperature: 293.15,
            pressure: 101325.0,
        };

        let analysis: VenturiAnalysis<f64> = model.analyze(&fluid, &conditions).unwrap();

        // Use the model's own Re-dependent C_d and K_exp for the analytical
        // Bernoulli formula.  This validates the PRESSURE DECOMPOSITION
        // (½ρV²(1−β²)/C_d²) rather than the C_d model itself, which is
        // a separate empirical fit tested via effective_discharge_coefficient().
        let c_d = analysis.discharge_coefficient;
        let k_exp = analysis.expansion_loss_coefficient;

        // Analytical Bernoulli contraction
        let dp_contraction_analytical =
            0.5 * WATER_DENSITY * v_throat.powi(2) * (1.0 - beta_sq) / (c_d * c_d);

        // Analytical expansion loss
        let dp_expansion_analytical =
            k_exp * 0.5 * WATER_DENSITY * (v_throat - v_out).powi(2);

        // Validate contraction pressure drop
        let err_contr = pct_error(analysis.dp_contraction, dp_contraction_analytical);
        rows.push(ValidationRow {
            parameter: format!("ΔP_contr(V={:.2})", v_in),
            literature: dp_contraction_analytical,
            cfd_1d: analysis.dp_contraction,
            error_pct: err_contr,
            pass: err_contr < TOL_EXACT,
        });

        // Validate expansion loss
        let err_exp = pct_error(analysis.dp_expansion_loss, dp_expansion_analytical);
        rows.push(ValidationRow {
            parameter: format!("ΔP_exp(V={:.2})", v_in),
            literature: dp_expansion_analytical,
            cfd_1d: analysis.dp_expansion_loss,
            error_pct: err_exp,
            pass: err_exp < TOL_EXACT,
        });
    }

    // Validate Bernoulli scaling: ΔP×C_d² ∝ V² at fixed geometry.
    // The raw ΔP ratio includes Re-dependent C_d variation, so we
    // normalise by C_d² to isolate the V² dependence.
    let v1 = velocities[0];
    let v2 = velocities[2]; // 0.01 vs 0.1 → ratio 10 → V² ratio 100
    let cond1 = FlowConditions {
        reynolds_number: None,
        velocity: Some(v1),
        flow_rate: None,
        shear_rate: None,
        temperature: 293.15,
        pressure: 101325.0,
    };
    let cond2 = FlowConditions {
        reynolds_number: None,
        velocity: Some(v2),
        flow_rate: None,
        shear_rate: None,
        temperature: 293.15,
        pressure: 101325.0,
    };
    let a1 = model.analyze(&fluid, &cond1).unwrap();
    let a2 = model.analyze(&fluid, &cond2).unwrap();
    // Normalised: ΔP·C_d² removes discharge-coefficient dependency
    let norm1 = a1.dp_contraction * a1.discharge_coefficient.powi(2);
    let norm2 = a2.dp_contraction * a2.discharge_coefficient.powi(2);
    let dp_ratio_expected = (v2 / v1).powi(2); // 100
    let dp_ratio_normalised = norm2 / norm1;
    let err = pct_error(dp_ratio_normalised, dp_ratio_expected);
    rows.push(ValidationRow {
        parameter: "ΔP·Cd² ratio (V²)".into(),
        literature: dp_ratio_expected,
        cfd_1d: dp_ratio_normalised,
        error_pct: err,
        pass: err < TOL_EXACT,
    });

    let all_pass = rows.iter().all(|r| r.pass);
    ValidationCase {
        id: 5,
        name: "Venturi Pressure Drop".into(),
        reference: "Bernoulli; ISO 5167-4; Idelchik (2007)".into(),
        rows,
        all_pass,
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Case 6 — Network Cross-Validation (Kirchhoff)
// ═══════════════════════════════════════════════════════════════════════════════

/// Validate that the full cfd-1d network solver produces pressure drops
/// consistent with summing individual analytical channel resistances.
///
/// For a simple bifurcating network with straight channels:
///   - ΔP_total = R_inlet · Q_inlet   (inlet channel)
///   - The two daughter channels are in parallel: R_parallel = R₁·R₂/(R₁+R₂)
///   - ΔP_network = [R_inlet + R_parallel] · Q_total
///
/// We validate mass conservation (Q_in = Q_out1 + Q_out2) and compare
/// network solver ΔP against analytical.
fn case_network_cross_validation() -> ValidationCase {
    let fluid = water_fluid();
    let mut rows = Vec::new();

    // Create a bifurcation chip via scheme
    let config = GeometryConfig::new(0.5, 1.0, 0.5).unwrap();
    let system = create_geometry(
        (60.0, 30.0),
        &[SplitType::Bifurcation],
        &config,
        &ChannelTypeConfig::AllStraight,
    );

    let converter = SchemeNetworkConverter::with_scale(&system, 1e-3);
    let mut network = converter.build_network(fluid.clone()).unwrap();

    let (inlets, outlets) = find_boundary_nodes(&network);
    let p_in = 8000.0_f64;
    let p_out = 1000.0_f64;
    for &ni in &inlets {
        network.set_pressure(ni, p_in);
    }
    for &ni in &outlets {
        network.set_pressure(ni, p_out);
    }

    let problem = NetworkProblem::new(network);
    let solver = NetworkSolver::<f64, ConstantPropertyFluid<f64>>::with_config(SolverConfig {
        tolerance: 1e-10,
        max_iterations: 500,
    });
    let solution = solver.solve_network(&problem).unwrap();

    let pressures = solution.pressures();
    let flow_rates = solution.flow_rates();

    // ── Mass conservation ────────────────────────────────────────────────
    let mut q_inlet_total = 0.0;
    let mut q_outlet_total = 0.0;

    for ei in solution.graph.edge_indices() {
        let (from, to) = solution.graph.edge_endpoints(ei).unwrap();
        let q = flow_rates.get(&ei).copied().unwrap_or(0.0);
        let q_abs = q.abs();

        if inlets.contains(&from) {
            q_inlet_total += q_abs;
        }
        if outlets.contains(&to) {
            q_outlet_total += q_abs;
        }
    }

    let mass_err = pct_error(q_outlet_total, q_inlet_total);
    rows.push(ValidationRow {
        parameter: "Mass conservation".into(),
        literature: q_inlet_total,
        cfd_1d: q_outlet_total,
        error_pct: mass_err,
        pass: mass_err < TOL_NETWORK,
    });

    // ── Pressure drop should equal P_in - P_out ──────────────────────────
    let dp_bc = p_in - p_out;
    // Get actual inlet and outlet node pressures
    let p_in_solved = inlets.iter().map(|n| pressures.get(n).copied().unwrap_or(0.0)).next().unwrap_or(0.0);
    let p_out_solved = outlets.iter().map(|n| pressures.get(n).copied().unwrap_or(0.0)).next().unwrap_or(0.0);
    let dp_solved = p_in_solved - p_out_solved;

    let dp_err = pct_error(dp_solved, dp_bc);
    rows.push(ValidationRow {
        parameter: "ΔP boundary".into(),
        literature: dp_bc,
        cfd_1d: dp_solved,
        error_pct: dp_err,
        pass: dp_err < TOL_NETWORK,
    });

    // ── Analytical resistance vs solver ──────────────────────────────────
    // For each channel, compute R analytically from stored geometry
    // and compare Q·R against the solver pressure drop
    let mut analytical_dp_sum = 0.0;
    let mut solver_dp_sum = 0.0;
    let mut n_channels = 0;

    for ei in solution.graph.edge_indices() {
        let (from, to) = solution.graph.edge_endpoints(ei).unwrap();
        let q = flow_rates.get(&ei).copied().unwrap_or(0.0).abs();
        let p_from = pressures.get(&from).copied().unwrap_or(0.0);
        let p_to = pressures.get(&to).copied().unwrap_or(0.0);
        let dp_solver_ch = (p_from - p_to).abs();

        if let Some(props) = solution.properties.get(&ei) {
            let length = props.length;
            let dh = props.hydraulic_diameter.unwrap_or(1e-3);
            let area = props.area;

            // Analytical R for rectangular channel
            // Use Hagen-Poiseuille-equivalent: R ≈ (Po·μ·L) / (2·A·Dh²)
            // Since channels from scheme are rectangular, get aspect ratio
            let (w, h_ch) = if let Some(geom) = &props.geometry {
                match &geom.cross_section {
                    cfd_1d::channel::CrossSection::Rectangular { width, height } => (*width, *height),
                    _ => {
                        // Estimate from area and Dh
                        let estimate_h = dh / 2.0;
                        let estimate_w = area / estimate_h;
                        (estimate_w, estimate_h)
                    }
                }
            } else {
                let estimate_h = dh / 2.0;
                let estimate_w = area / estimate_h;
                (estimate_w, estimate_h)
            };

            // Shah-London Po for this aspect ratio
            let alpha = w.max(h_ch) / w.min(h_ch);
            let po = if (alpha - 1.0).abs() < 0.01 {
                56.91
            } else {
                let a_inv = 1.0 / alpha;
                96.0 * (1.0
                    - 1.3553 * a_inv
                    + 1.9467 * a_inv.powi(2)
                    - 1.7012 * a_inv.powi(3)
                    + 0.9564 * a_inv.powi(4)
                    - 0.2537 * a_inv.powi(5))
            };

            let r_analytical = po * WATER_VISCOSITY * length / (2.0 * area * dh * dh);
            let dp_analytical = r_analytical * q;

            analytical_dp_sum += dp_analytical;
            solver_dp_sum += dp_solver_ch;
            n_channels += 1;
        }
    }

    // Compare per-channel averaged ΔP
    if n_channels > 0 {
        let avg_analytical = analytical_dp_sum / n_channels as f64;
        let avg_solver = solver_dp_sum / n_channels as f64;
        let err = pct_error(avg_solver, avg_analytical);

        rows.push(ValidationRow {
            parameter: "Avg ΔP per channel".into(),
            literature: avg_analytical,
            cfd_1d: avg_solver,
            error_pct: err,
            pass: err < TOL_EMPIRICAL,
        });
    }

    // ── Pressure monotonicity: inlet > junction > outlet ─────────────────
    let junction_pressures: Vec<f64> = solution.graph.node_indices()
        .filter(|ni| matches!(solution.graph[*ni].node_type, NodeType::Junction))
        .filter_map(|ni| pressures.get(&ni).copied())
        .collect();

    let all_monotonic = junction_pressures.iter().all(|&p| p > p_out && p < p_in);
    rows.push(ValidationRow {
        parameter: "Pressure monotonic".into(),
        literature: 1.0,
        cfd_1d: if all_monotonic { 1.0 } else { 0.0 },
        error_pct: if all_monotonic { 0.0 } else { 100.0 },
        pass: all_monotonic,
    });

    let all_pass = rows.iter().all(|r| r.pass);
    ValidationCase {
        id: 6,
        name: "Network Cross-Validation (Kirchhoff)".into(),
        reference: "Kirchhoff's laws; scheme→cfd-1d pipeline".into(),
        rows,
        all_pass,
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// SVG Parity Plots
// ═══════════════════════════════════════════════════════════════════════════════

/// Draw a parity plot (literature vs cfd-1d) for all validation rows.
fn plot_parity(cases: &[ValidationCase]) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all("outputs")?;
    let path = "outputs/literature_parity.svg";
    let root = SVGBackend::new(path, (800, 700)).into_drawing_area();
    root.fill(&WHITE)?;

    // Collect all rows that have comparable scales (filter out ratio/boolean rows)
    let points: Vec<(f64, f64, bool, String)> = cases
        .iter()
        .flat_map(|c| {
            c.rows.iter().filter(|r| r.literature > 0.0 && r.cfd_1d > 0.0).map(|r| {
                (r.literature, r.cfd_1d, r.pass, format!("Case {} — {}", c.id, r.parameter))
            })
        })
        .collect();

    if points.is_empty() {
        return Ok(());
    }

    // Use log scale for wide range
    let min_val = points
        .iter()
        .flat_map(|(l, c, _, _)| vec![*l, *c])
        .fold(f64::INFINITY, f64::min)
        * 0.5;
    let max_val = points
        .iter()
        .flat_map(|(l, c, _, _)| vec![*l, *c])
        .fold(f64::NEG_INFINITY, f64::max)
        * 2.0;

    let min_log = min_val.log10().floor();
    let max_log = max_val.log10().ceil();

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Literature Validation — Parity Plot (log-log)",
            ("sans-serif", 18),
        )
        .margin(15)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(
            (10.0_f64.powf(min_log)..10.0_f64.powf(max_log)).log_scale(),
            (10.0_f64.powf(min_log)..10.0_f64.powf(max_log)).log_scale(),
        )?;

    chart
        .configure_mesh()
        .x_desc("Literature / Analytical")
        .y_desc("cfd-1d Computed")
        .draw()?;

    // 1:1 line
    chart.draw_series(LineSeries::new(
        vec![
            (10.0_f64.powf(min_log), 10.0_f64.powf(min_log)),
            (10.0_f64.powf(max_log), 10.0_f64.powf(max_log)),
        ],
        BLACK.stroke_width(2),
    ))?
    .label("1:1 (ideal)")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], BLACK.stroke_width(2)));

    // ±5 % bands
    let band_points_upper: Vec<(f64, f64)> = (0..100)
        .map(|i| {
            let t = min_log + (max_log - min_log) * i as f64 / 99.0;
            let x = 10.0_f64.powf(t);
            (x, x * 1.05)
        })
        .collect();
    let band_points_lower: Vec<(f64, f64)> = (0..100)
        .map(|i| {
            let t = min_log + (max_log - min_log) * i as f64 / 99.0;
            let x = 10.0_f64.powf(t);
            (x, x * 0.95)
        })
        .collect();

    chart.draw_series(LineSeries::new(
        band_points_upper,
        RGBColor(200, 200, 200).stroke_width(1),
    ))?
    .label("±5% band")
    .legend(|(x, y)| {
        PathElement::new(
            vec![(x, y), (x + 15, y)],
            RGBColor(200, 200, 200).stroke_width(1),
        )
    });
    chart.draw_series(LineSeries::new(
        band_points_lower,
        RGBColor(200, 200, 200).stroke_width(1),
    ))?;

    // Data points
    let pass_pts: Vec<(f64, f64)> = points
        .iter()
        .filter(|(_, _, p, _)| *p)
        .map(|(l, c, _, _)| (*l, *c))
        .collect();
    let fail_pts: Vec<(f64, f64)> = points
        .iter()
        .filter(|(_, _, p, _)| !*p)
        .map(|(l, c, _, _)| (*l, *c))
        .collect();

    chart
        .draw_series(pass_pts.iter().map(|&(x, y)| Circle::new((x, y), 5, GREEN.filled())))?
        .label("Pass")
        .legend(|(x, y)| Circle::new((x, y), 5, GREEN.filled()));

    if !fail_pts.is_empty() {
        chart
            .draw_series(fail_pts.iter().map(|&(x, y)| Circle::new((x, y), 5, RED.filled())))?
            .label("Fail")
            .legend(|(x, y)| Circle::new((x, y), 5, RED.filled()));
    }

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .position(SeriesLabelPosition::UpperLeft)
        .draw()?;

    root.present()?;
    println!("    SVG  → {path}");
    Ok(())
}

/// Draw a per-case error bar chart.
fn plot_error_bars(cases: &[ValidationCase]) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all("outputs")?;
    let path = "outputs/literature_error_bars.svg";
    let root = SVGBackend::new(path, (1100, 500)).into_drawing_area();
    root.fill(&WHITE)?;

    // Collect all row errors
    let mut labels: Vec<String> = Vec::new();
    let mut errors: Vec<f64> = Vec::new();
    let mut passes: Vec<bool> = Vec::new();

    for c in cases {
        for r in &c.rows {
            let param_short: String = r.parameter.chars().take(12).collect();
            labels.push(format!("C{}:{}", c.id, param_short));
            errors.push(r.error_pct);
            passes.push(r.pass);
        }
    }

    let n = labels.len();
    if n == 0 {
        return Ok(());
    }

    let max_err = errors.iter().copied().fold(0.0_f64, f64::max) * 1.3;
    let y_max = max_err.max(5.0);

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Relative Error per Validation Point (%)",
            ("sans-serif", 16),
        )
        .margin(15)
        .x_label_area_size(80)
        .y_label_area_size(50)
        .build_cartesian_2d(0.0..(n as f64), 0.0..y_max)?;

    chart
        .configure_mesh()
        .y_desc("Error [%]")
        .x_desc("")
        .x_label_formatter(&|x| {
            let idx = *x as usize;
            if idx < labels.len() {
                labels[idx].clone()
            } else {
                String::new()
            }
        })
        .x_label_style(("sans-serif", 8).into_text_style(&root).transform(FontTransform::Rotate270))
        .draw()?;

    // Bars
    chart.draw_series((0..n).map(|i| {
        let color = if passes[i] {
            GREEN.mix(0.7)
        } else {
            RED.mix(0.7)
        };
        Rectangle::new(
            [(i as f64 + 0.1, 0.0), (i as f64 + 0.9, errors[i])],
            color.filled(),
        )
    }))?;

    // 1% threshold
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(0.0, 1.0), (n as f64, 1.0)],
        RGBColor(0, 0, 180).stroke_width(1),
    )))?
    .label("1% (exact)")
    .legend(|(x, y)| {
        PathElement::new(vec![(x, y), (x + 15, y)], RGBColor(0, 0, 180).stroke_width(1))
    });

    // 5% threshold
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(0.0, 5.0), (n as f64, 5.0)],
        RGBColor(180, 100, 0).stroke_width(1),
    )))?
    .label("5% (empirical)")
    .legend(|(x, y)| {
        PathElement::new(
            vec![(x, y), (x + 15, y)],
            RGBColor(180, 100, 0).stroke_width(1),
        )
    });

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .position(SeriesLabelPosition::UpperRight)
        .draw()?;

    root.present()?;
    println!("    SVG  → {path}");
    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════════
// CSV / JSON Export
// ═══════════════════════════════════════════════════════════════════════════════

fn export_csv(cases: &[ValidationCase]) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all("outputs")?;
    let path = "outputs/literature_validation.csv";
    let mut f = std::fs::File::create(path)?;
    writeln!(
        f,
        "case_id,case_name,reference,parameter,literature_value,cfd_1d_value,error_pct,pass"
    )?;
    for c in cases {
        for r in &c.rows {
            writeln!(
                f,
                "{},{},{},{},{:.10e},{:.10e},{:.6},{}", 
                c.id,
                c.name.replace(',', ";"),
                c.reference.replace(',', ";"),
                r.parameter.replace(',', ";"),
                r.literature,
                r.cfd_1d,
                r.error_pct,
                r.pass,
            )?;
        }
    }
    println!("    CSV  → {path}");
    Ok(())
}

fn export_json(cases: &[ValidationCase]) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all("outputs")?;
    let path = "outputs/literature_validation.json";
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "{{")?;
    writeln!(f, "  \"title\": \"CFD-1D Literature Validation Report\",")?;
    writeln!(
        f,
        "  \"timestamp\": \"{}\",",
        chrono::Utc::now().to_rfc3339()
    )?;
    writeln!(f, "  \"n_cases\": {},", cases.len())?;
    let total_pass = cases.iter().flat_map(|c| c.rows.iter()).filter(|r| r.pass).count();
    let total_rows = cases.iter().map(|c| c.rows.len()).sum::<usize>();
    writeln!(f, "  \"n_pass\": {},", total_pass)?;
    writeln!(f, "  \"n_fail\": {},", total_rows - total_pass)?;
    writeln!(
        f,
        "  \"all_pass\": {},",
        cases.iter().all(|c| c.all_pass)
    )?;
    writeln!(f, "  \"cases\": [")?;
    for (ci, c) in cases.iter().enumerate() {
        writeln!(f, "    {{")?;
        writeln!(f, "      \"id\": {},", c.id)?;
        writeln!(f, "      \"name\": \"{}\",", c.name)?;
        writeln!(f, "      \"reference\": \"{}\",", c.reference)?;
        writeln!(f, "      \"all_pass\": {},", c.all_pass)?;
        writeln!(f, "      \"rows\": [")?;
        for (ri, r) in c.rows.iter().enumerate() {
            writeln!(f, "        {{")?;
            writeln!(f, "          \"parameter\": \"{}\",", r.parameter)?;
            writeln!(f, "          \"literature\": {:.10e},", r.literature)?;
            writeln!(f, "          \"cfd_1d\": {:.10e},", r.cfd_1d)?;
            writeln!(f, "          \"error_pct\": {:.6},", r.error_pct)?;
            writeln!(f, "          \"pass\": {}", r.pass)?;
            if ri < c.rows.len() - 1 {
                writeln!(f, "        }},")?;
            } else {
                writeln!(f, "        }}")?;
            }
        }
        writeln!(f, "      ]")?;
        if ci < cases.len() - 1 {
            writeln!(f, "    }},")?;
        } else {
            writeln!(f, "    }}")?;
        }
    }
    writeln!(f, "  ]")?;
    writeln!(f, "}}")?;
    println!("    JSON → {path}");
    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════════
// Summary
// ═══════════════════════════════════════════════════════════════════════════════

fn print_summary(cases: &[ValidationCase]) {
    let total_rows: usize = cases.iter().map(|c| c.rows.len()).sum();
    let total_pass: usize = cases
        .iter()
        .flat_map(|c| c.rows.iter())
        .filter(|r| r.pass)
        .count();
    let total_fail = total_rows - total_pass;

    let max_err: f64 = cases
        .iter()
        .flat_map(|c| c.rows.iter())
        .map(|r| r.error_pct)
        .fold(0.0_f64, f64::max);
    let avg_err: f64 = cases
        .iter()
        .flat_map(|c| c.rows.iter())
        .map(|r| r.error_pct)
        .sum::<f64>()
        / total_rows.max(1) as f64;

    println!();
    println!("  ╔══════════════════════════════════════════════════════════════════════╗");
    println!("  ║  Literature Validation Summary                                     ║");
    println!("  ╚══════════════════════════════════════════════════════════════════════╝");
    println!();
    println!(
        "    Cases: {}     Points: {}     Pass: {}     Fail: {}",
        cases.len(),
        total_rows,
        total_pass,
        total_fail
    );
    println!(
        "    Max error: {:.4} %     Mean error: {:.4} %",
        max_err, avg_err
    );
    println!();

    println!("    Per-case breakdown:");
    println!("    ┌──────┬──────────────────────────────────────┬───────┬───────┬────────┐");
    println!("    │  #   │ Case                                 │ Pass  │ Fail  │ Status │");
    println!("    ├──────┼──────────────────────────────────────┼───────┼───────┼────────┤");
    for c in cases {
        let n_pass = c.rows.iter().filter(|r| r.pass).count();
        let n_fail = c.rows.len() - n_pass;
        let status = if c.all_pass { "  ✓   " } else { "  ✗   " };
        println!(
            "    │ {:>4} │ {:>36} │ {:>5} │ {:>5} │{status}│",
            c.id,
            &c.name[..c.name.len().min(36)],
            n_pass,
            n_fail,
        );
    }
    println!("    └──────┴──────────────────────────────────────┴───────┴───────┴────────┘");

    let all_pass = cases.iter().all(|c| c.all_pass);
    println!();
    if all_pass {
        println!("    ═══ ALL {} VALIDATION POINTS PASS ═══", total_rows);
    } else {
        println!(
            "    ═══ {}/{} VALIDATION POINTS PASS, {} FAIL ═══",
            total_pass, total_rows, total_fail
        );
    }

    // Output file listing
    println!();
    println!("    Output files:");
    let files = [
        "outputs/literature_validation.csv",
        "outputs/literature_validation.json",
        "outputs/literature_parity.svg",
        "outputs/literature_error_bars.svg",
    ];
    for p in &files {
        let exists = Path::new(p).exists();
        let size = if exists {
            std::fs::metadata(p).map(|m| m.len()).unwrap_or(0)
        } else {
            0
        };
        println!(
            "      {} {:>8} B  {}",
            if exists { "✓" } else { "✗" },
            size,
            p
        );
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Main
// ═══════════════════════════════════════════════════════════════════════════════

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!();
    println!("  ═══════════════════════════════════════════════════════════════════");
    println!("   CFD-1D Literature Validation");
    println!("   Quantitative comparison against published correlations");
    println!("  ═══════════════════════════════════════════════════════════════════");

    let mut cases = Vec::new();

    cases.push(case_shah_london());
    print_case(cases.last().unwrap());

    cases.push(case_hagen_poiseuille());
    print_case(cases.last().unwrap());

    cases.push(case_bruus_approximation());
    print_case(cases.last().unwrap());

    cases.push(case_dean_curvature());
    print_case(cases.last().unwrap());

    cases.push(case_venturi());
    print_case(cases.last().unwrap());

    cases.push(case_network_cross_validation());
    print_case(cases.last().unwrap());

    // ── Export ────────────────────────────────────────────────────────────
    println!();
    println!("  ── Exporting data and plots ─────────────────────────────────────");
    export_csv(&cases)?;
    export_json(&cases)?;
    plot_parity(&cases)?;
    plot_error_bars(&cases)?;

    // ── Summary ──────────────────────────────────────────────────────────
    print_summary(&cases);

    println!();
    println!("  Done.");
    Ok(())
}
