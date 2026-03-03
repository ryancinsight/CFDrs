//! End-to-end Milestone 12 pipeline for the top 5 RbcProtectedSdt designs.
//!
//! Drives the 5 highest-ranked designs through the full multi-fidelity ladder
//! and writes report-ready artefacts to `report/milestone12/` (gitignored):
//!
//! | File | Contents |
//! |------|----------|
//! | `top5_designs.json` | Full ranked list with SdtMetrics |
//! | `{id}_1d_metrics.json` | 1D network-solver SdtMetrics |
//! | `{id}_schematic.svg` | Chip layout SVG |
//! | `{id}_2d_venturi.json` | 2D FVM throat confirmation (venturi only) |
//! | `{id}_3d_result.json` | 3D FEM summary (venturi only) |
//! | `{id}_validation.json` | 1D ↔ 2D ↔ 3D consistency table |
//! | `pipeline_summary.md` | Human-readable summary table |
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example optim_pipeline --no-default-features
//! ```

use cfd_2d::solvers::ns_fvm::BloodModel;
use cfd_2d::solvers::venturi_flow::{VenturiGeometry, VenturiSolver2D};
use cfd_3d::venturi::{VenturiConfig3D, VenturiSolver3D};
use cfd_core::physics::fluid::blood::CarreauYasudaBlood;
use cfd_core::physics::fluid::CassonBlood;
use cfd_mesh::VenturiMeshBuilder;
use cfd_optim::{save_schematic_svg, save_top5_json, OptimMode, SdtOptimizer, SdtWeights};
use serde::Serialize;
use std::path::Path;

// ── Physical constants ────────────────────────────────────────────────────────

/// Blood density [kg/m³]
const RHO: f64 = 1_060.0;
/// Atmospheric pressure [Pa]
const P_ATM_PA: f64 = 101_325.0;
/// Blood vapour pressure at 37 °C [Pa]
const P_VAPOR_PA: f64 = 6_280.0;

// ── Output structures ─────────────────────────────────────────────────────────

/// 2D FVM venturi throat confirmation result.
#[derive(Debug, Serialize)]
struct Venturi2DResult {
    /// Inlet velocity used by the 2D solver [m/s]
    u_inlet_m_s: f64,
    /// 2D FVM throat velocity [m/s]
    u_throat_m_s: f64,
    /// 2D FVM pressure drop inlet → throat [Pa]  (positive = drop)
    dp_throat_pa: f64,
    /// 2D FVM pressure recovery [Pa]
    dp_recovery_pa: f64,
    /// Cavitation number σ anchored to the candidate's absolute inlet pressure
    sigma_2d: f64,
    /// Bernoulli 1D analytical pressure drop for comparison [Pa]
    dp_bernoulli_1d_pa: f64,
}

/// 3D FEM venturi throat confirmation result.
#[derive(Debug, Serialize)]
struct Venturi3DResult {
    /// 3D FEM inlet velocity [m/s]
    u_inlet_m_s: f64,
    /// 3D FEM throat velocity [m/s]
    u_throat_m_s: f64,
    /// 3D FEM pressure drop inlet → throat [Pa]
    dp_throat_pa: f64,
    /// 3D FEM pressure recovery [Pa]
    dp_recovery_pa: f64,
    /// Relative mass balance error (|Q_in - Q_out| / Q_in)
    mass_error: f64,
    /// Mesh resolution used (axial × radial)
    resolution: (usize, usize),
}

/// 1D ↔ 2D ↔ 3D consistency validation row for one design.
#[derive(Debug, Serialize)]
struct ValidationRow {
    /// Design ID
    id: String,
    /// Topology short name
    topology: String,
    /// Bernoulli 1D dp_throat [Pa]
    dp_1d_bernoulli_pa: f64,
    /// 2D FVM dp_throat [Pa]; NaN if not computed
    dp_2d_fvm_pa: f64,
    /// 3D FEM dp_throat [Pa]; NaN if not computed
    dp_3d_fem_pa: f64,
    /// Relative agreement |dp_1d - dp_2d| / dp_1d [%]; NaN if unavailable
    agreement_1d_2d_pct: f64,
    /// Relative agreement |dp_2d - dp_3d| / dp_2d [%]; NaN if unavailable
    agreement_2d_3d_pct: f64,
    /// 3D mass conservation error [%]; NaN if unavailable
    mass_error_3d_pct: f64,
    /// Cavitation number σ from 1D optimizer
    sigma_1d: f64,
    /// Cavitation number σ from 2D FVM; NaN if unavailable
    sigma_2d: f64,
    /// RbcProtectedSdt score [0, 1]
    score: f64,
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Relative percent deviation between two values, returning NaN when denom ≈ 0.
fn pct_diff(a: f64, b: f64) -> f64 {
    let denom = a.abs();
    if denom < 1e-12 {
        f64::NAN
    } else {
        (a - b).abs() / denom * 100.0
    }
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out_dir = Path::new("report/milestone12");
    std::fs::create_dir_all(out_dir)?;

    println!("=== Milestone 12 Multi-Fidelity Pipeline ===\n");
    println!("Output directory: {}", out_dir.display());

    // ── Part 1: Obtain top-5 RbcProtectedSdt designs ─────────────────────────
    println!("\n[1/6] Running RbcProtectedSdt parametric optimizer …");
    let optimizer = SdtOptimizer::new(OptimMode::RbcProtectedSdt, SdtWeights::default());
    let top5 = optimizer.top_k(5)?;

    println!("      Top-5 designs:");
    for d in &top5 {
        println!(
            "      {}  score={:.4}  σ={:.3}  HI={:.2e}",
            d.candidate.id,
            d.score,
            d.metrics.cavitation_number,
            d.metrics.hemolysis_index_per_pass
        );
    }

    save_top5_json(&top5, &out_dir.join("top5_designs.json"))?;
    println!("      Saved: top5_designs.json");

    // ── Part 2: Per-design 1D metrics and schematic SVG ───────────────────────
    println!("\n[2/6] Saving 1D metrics and schematic SVGs …");
    for d in &top5 {
        let id = &d.candidate.id;
        let json = serde_json::to_string_pretty(&d.metrics)?;
        std::fs::write(out_dir.join(format!("{id}_1d_metrics.json")), json)?;

        match save_schematic_svg(&d.candidate, &out_dir.join(format!("{id}_schematic.svg"))) {
            Ok(()) => {}
            Err(e) => eprintln!("      [WARN] schematic SVG failed for {id}: {e}"),
        }
        println!("      {id}: 1D metrics + SVG saved");
    }

    // ── Parts 3-5: 2D/3D confirmation per venturi design ─────────────────────
    println!("\n[3/6] Running 2D FVM + 3D FEM venturi confirmation …");
    let mut validation_rows: Vec<ValidationRow> = Vec::new();

    for d in &top5 {
        let c = &d.candidate;
        let m = &d.metrics;
        let id = &c.id;

        if !c.topology.has_venturi() {
            println!("      {id}: no venturi — skipping 2D/3D");
            validation_rows.push(ValidationRow {
                id: id.clone(),
                topology: c.topology.short().to_string(),
                dp_1d_bernoulli_pa: f64::NAN,
                dp_2d_fvm_pa: f64::NAN,
                dp_3d_fem_pa: f64::NAN,
                agreement_1d_2d_pct: f64::NAN,
                agreement_2d_3d_pct: f64::NAN,
                mass_error_3d_pct: f64::NAN,
                sigma_1d: m.cavitation_number,
                sigma_2d: f64::NAN,
                score: d.score,
            });
            continue;
        }

        // ── Derived geometry quantities ───────────────────────────────────────
        let q = c.per_venturi_flow();
        let a_in = c.inlet_area_m2();
        let a_th = c.throat_area_m2();
        let v_in_1d = q / a_in.max(1e-30);
        let v_th_1d = q / a_th.max(1e-30);
        let dp_bernoulli = 0.5 * RHO * (v_th_1d * v_th_1d - v_in_1d * v_in_1d);

        // ── Part 3: 2D FVM confirmation ───────────────────────────────────────
        let geom2d = VenturiGeometry::<f64>::new(
            c.inlet_diameter_m,  // w_inlet
            c.throat_diameter_m, // w_throat
            3e-3,                // l_inlet  (3 mm fixed approach)
            2e-3,                // l_converge (2 mm linear taper)
            c.throat_length_m,   // l_throat
            4e-3,                // l_diverge (4 mm recovery)
            c.channel_height_m,  // height (rectangular cross-section)
        );

        // Use Casson blood model for physiological accuracy.
        let blood_2d = BloodModel::Casson(CassonBlood::<f64>::normal_blood());

        // Adaptive grid: stretch y-cells toward centre to resolve throat.
        let cr = c.inlet_diameter_m / c.throat_diameter_m.max(1e-12);
        let ny_2d = (4.0 * cr).round().clamp(40.0, 200.0) as usize;
        let beta_2d =
            (1.0 - 4.0 * c.throat_diameter_m / c.inlet_diameter_m.max(1e-12)).clamp(0.0, 0.9);
        let mut solver2d =
            VenturiSolver2D::new_stretched(geom2d, blood_2d, RHO, 60, ny_2d, beta_2d);

        // Map volumetric flow to 2D rectangular cross-section velocity.
        let area_2d = c.inlet_diameter_m * c.channel_height_m;
        let u_inlet_2d = q / area_2d.max(1e-30);

        let (dp_2d, sigma_2d, result2d_opt) = match solver2d.solve(u_inlet_2d) {
            Ok(sol) => {
                let p_abs_inlet = P_ATM_PA + c.inlet_gauge_pa;
                // dp_throat from VenturiFlowSolution is negative (inlet > throat)
                let p_abs_throat = p_abs_inlet + sol.dp_throat;
                let dyn_p = 0.5 * RHO * sol.u_throat * sol.u_throat;
                let sigma = if dyn_p > 1e-12 {
                    (p_abs_throat - P_VAPOR_PA) / dyn_p
                } else {
                    f64::INFINITY
                };
                let dp = -sol.dp_throat; // positive pressure drop
                (
                    dp,
                    sigma,
                    Some(Venturi2DResult {
                        u_inlet_m_s: u_inlet_2d,
                        u_throat_m_s: sol.u_throat,
                        dp_throat_pa: dp,
                        dp_recovery_pa: -sol.dp_recovery,
                        sigma_2d: sigma,
                        dp_bernoulli_1d_pa: dp_bernoulli,
                    }),
                )
            }
            Err(e) => {
                eprintln!("      [WARN] 2D FVM failed for {id}: {e}");
                (f64::NAN, f64::NAN, None)
            }
        };

        if let Some(ref r2d) = result2d_opt {
            let json = serde_json::to_string_pretty(r2d)?;
            std::fs::write(out_dir.join(format!("{id}_2d_venturi.json")), json)?;
            println!(
                "      {id}: 2D  dp={:.1} Pa  σ={:.3}  (Bernoulli: {:.1} Pa)",
                r2d.dp_throat_pa, r2d.sigma_2d, dp_bernoulli
            );
        }

        // ── Part 4: 3D FEM confirmation ───────────────────────────────────────
        let d_in = c.inlet_diameter_m;
        let d_th = c.throat_diameter_m;
        let res3d = (60_usize, 10_usize); // adequate for M12 confirmation

        let builder3d = VenturiMeshBuilder::<f64>::new(
            d_in,
            d_th,
            5.0 * d_in, // l_inlet  (5 inlet-diam approach)
            3.0 * d_in, // l_convergent
            c.throat_length_m,
            7.0 * d_in, // l_divergent
            5.0 * d_in, // l_outlet
        )
        .with_resolution(res3d.0, res3d.1)
        .with_circular(false);

        let config3d = VenturiConfig3D::<f64> {
            inlet_flow_rate: q,
            resolution: res3d,
            circular: false,
            rect_height: Some(c.channel_height_m),
            ..Default::default()
        };

        let solver3d = VenturiSolver3D::new(builder3d, config3d);
        let fluid3d = CarreauYasudaBlood::<f64>::normal_blood();

        let (dp_3d, mass_err_3d, result3d_opt) = match solver3d.solve(fluid3d) {
            Ok(sol) => {
                let dp = sol.dp_throat.abs();
                let mass = sol.mass_error.abs();
                (
                    dp,
                    mass,
                    Some(Venturi3DResult {
                        u_inlet_m_s: sol.u_inlet,
                        u_throat_m_s: sol.u_throat,
                        dp_throat_pa: dp,
                        dp_recovery_pa: sol.dp_recovery.abs(),
                        mass_error: mass,
                        resolution: res3d,
                    }),
                )
            }
            Err(e) => {
                eprintln!("      [WARN] 3D FEM failed for {id}: {e}");
                (f64::NAN, f64::NAN, None)
            }
        };

        if let Some(ref r3d) = result3d_opt {
            let json = serde_json::to_string_pretty(r3d)?;
            std::fs::write(out_dir.join(format!("{id}_3d_result.json")), json)?;
            println!(
                "      {id}: 3D  dp={:.1} Pa  u_th={:.4} m/s  mass_err={:.2e}",
                r3d.dp_throat_pa, r3d.u_throat_m_s, r3d.mass_error
            );
        }

        // ── Part 5: 1D ↔ 2D ↔ 3D consistency ────────────────────────────────
        let row = ValidationRow {
            id: id.clone(),
            topology: c.topology.short().to_string(),
            dp_1d_bernoulli_pa: dp_bernoulli,
            dp_2d_fvm_pa: dp_2d,
            dp_3d_fem_pa: dp_3d,
            agreement_1d_2d_pct: pct_diff(dp_bernoulli, dp_2d),
            agreement_2d_3d_pct: pct_diff(dp_2d, dp_3d),
            mass_error_3d_pct: mass_err_3d * 100.0,
            sigma_1d: m.cavitation_number,
            sigma_2d,
            score: d.score,
        };
        let json = serde_json::to_string_pretty(&row)?;
        std::fs::write(out_dir.join(format!("{id}_validation.json")), json)?;
        validation_rows.push(row);
    }

    // ── Part 6: Write pipeline summary report ─────────────────────────────────
    println!("\n[4/6] Building summary report …");
    let summary_path = out_dir.join("pipeline_summary.md");
    write_summary(&validation_rows, &top5, &summary_path)?;
    println!("      Saved: pipeline_summary.md");

    // ── Final console table ───────────────────────────────────────────────────
    println!("\n=== 1D ↔ 2D ↔ 3D Consistency Summary ===\n");
    println!(
        "{:<42} {:>12} {:>12} {:>12} {:>10} {:>10}",
        "Design ID", "dp_1D [Pa]", "dp_2D [Pa]", "dp_3D [Pa]", "1D-2D [%]", "2D-3D [%]"
    );
    println!("{}", "-".repeat(100));
    for row in &validation_rows {
        let fmt = |v: f64| {
            if v.is_nan() {
                "  —".to_string()
            } else {
                format!("{v:>12.1}")
            }
        };
        let fmtp = |v: f64| {
            if v.is_nan() {
                "  —".to_string()
            } else {
                format!("{v:>10.1}")
            }
        };
        println!(
            "{:<42} {} {} {} {} {}",
            row.id,
            fmt(row.dp_1d_bernoulli_pa),
            fmt(row.dp_2d_fvm_pa),
            fmt(row.dp_3d_fem_pa),
            fmtp(row.agreement_1d_2d_pct),
            fmtp(row.agreement_2d_3d_pct),
        );
    }
    println!("\nAll artefacts written to: {}", out_dir.display());

    Ok(())
}

// ── Summary report writer ─────────────────────────────────────────────────────

fn write_summary(
    rows: &[ValidationRow],
    designs: &[cfd_optim::RankedDesign],
    path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::fmt::Write as FmtWrite;
    let mut md = String::new();

    writeln!(md, "# Milestone 12 Multi-Fidelity Pipeline Summary\n")?;
    writeln!(
        md,
        "Generated by `cargo run -p cfd-optim --example optim_pipeline`\n"
    )?;

    writeln!(md, "## Top-5 Designs (RbcProtectedSdt)\n")?;
    writeln!(
        md,
        "| Rank | Design ID | Score | σ (1D) | HI | FDA | ΔP (Pa) |"
    )?;
    writeln!(
        md,
        "|------|-----------|-------|--------|-----|-----|---------|"
    )?;
    for d in designs {
        let m = &d.metrics;
        writeln!(
            md,
            "| {} | `{}` | {:.4} | {:.3} | {:.2e} | {} | {:.0} |",
            d.rank,
            d.candidate.id,
            d.score,
            m.cavitation_number,
            m.hemolysis_index_per_pass,
            if m.fda_main_compliant { "✓" } else { "✗" },
            m.total_pressure_drop_pa,
        )?;
    }

    writeln!(md, "\n## 1D ↔ 2D ↔ 3D Pressure-Drop Consistency\n")?;
    writeln!(
        md,
        "| Design ID | dp₁D [Pa] | dp₂D [Pa] | dp₃D [Pa] | 1D–2D [%] | 2D–3D [%] | Mass err [%] |"
    )?;
    writeln!(
        md,
        "|-----------|-----------|-----------|-----------|-----------|-----------|-------------|"
    )?;
    for row in rows {
        let f = |v: f64| {
            if v.is_nan() {
                "—".to_string()
            } else {
                format!("{v:.1}")
            }
        };
        writeln!(
            md,
            "| `{}` | {} | {} | {} | {} | {} | {} |",
            row.id,
            f(row.dp_1d_bernoulli_pa),
            f(row.dp_2d_fvm_pa),
            f(row.dp_3d_fem_pa),
            f(row.agreement_1d_2d_pct),
            f(row.agreement_2d_3d_pct),
            f(row.mass_error_3d_pct),
        )?;
    }

    writeln!(md, "\n## Cavitation Number (σ) Across Fidelities\n")?;
    writeln!(md, "| Design ID | σ₁D | σ₂D | Δσ [%] |")?;
    writeln!(md, "|-----------|-----|-----|--------|")?;
    for row in rows {
        let f = |v: f64| {
            if v.is_nan() || !v.is_finite() {
                "—".to_string()
            } else {
                format!("{v:.3}")
            }
        };
        writeln!(
            md,
            "| `{}` | {} | {} | {} |",
            row.id,
            f(row.sigma_1d),
            f(row.sigma_2d),
            f(pct_diff(row.sigma_1d, row.sigma_2d)),
        )?;
    }

    writeln!(md, "\n## Validation Criteria\n")?;
    writeln!(md, "- **Pressure agreement**: |dp₁D − dp₂D| / dp₁D < 20% (coarse 2D mesh, 2D vs 3D geometry)\n")?;
    writeln!(
        md,
        "- **3D mass conservation**: |Q_in − Q_out| / Q_in < 2%\n"
    )?;
    writeln!(md, "- **σ confirmation**: 2D σ within 30% of 1D σ (viscous effects not captured by Bernoulli)\n")?;

    std::fs::write(path, md)?;
    Ok(())
}
