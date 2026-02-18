//! Top-5 millifluidic design candidates for Sonodynamic Therapy (SDT).
//!
//! Runs the `cfd-optim` parametric sweep over the full 5-topology, ~800-candidate
//! design space and prints two ranked tables:
//!
//! 1. **SDT Cavitation mode** — venturi-throat designs ranked by hydrodynamic
//!    cavitation potential (σ < 1), all fitting within a 96-well plate and
//!    meeting FDA blood-shear requirements in sustained-flow channels.
//!
//! 2. **Uniform Exposure mode** — serpentine / bifurcation designs ranked by
//!    flow uniformity and residence time across the 6 × 6 centre treatment zone.
//!
//! Run with:
//! ```bash
//! cargo run -p cfd-optim --example sdt_top5
//! ```

use cfd_optim::{OptimMode, SdtOptimizer, SdtWeights};

fn main() {
    let weights = SdtWeights::default();

    // ── Mode 1: SDT Cavitation ────────────────────────────────────────────
    println!("\n{}", "=".repeat(110));
    println!("  cfd-optim  |  SDT Cavitation Mode  |  96-well plate, blood, FDA <= 150 Pa main channel");
    println!("{}", "=".repeat(110));
    println!("  Objective: Maximise hydrodynamic cavitation (sigma < 1) at venturi throats");
    println!("           + low haemolysis index (Giersiepen 1990) + well coverage");
    println!("  Constraint: total dP <= inlet gauge pressure; main channel shear <= 150 Pa");
    println!("{}", "-".repeat(110));

    match SdtOptimizer::new(OptimMode::SdtCavitation, weights).top_5() {
        Ok(designs) => {
            println!(
                "{:>4}  {:<42}  {:>7}  {:>9}  {:>9}  {:>7}  {:>7}  {:>7}  {:>7}",
                "#",
                "Candidate ID",
                "Score",
                "σ",
                "τ_throat",
                "FDA",
                "HI",
                "Cov%",
                "ΔP kPa"
            );
            println!("{}", "-".repeat(110));
            for d in &designs {
                let m = &d.metrics;
                let sigma_str = if m.cavitation_number.is_finite() {
                    format!("{:>9.3}", m.cavitation_number)
                } else {
                    format!("{:>9}", "—")
                };
                println!(
                    "{:>4}  {:<42}  {:>7.4}  {}  {:>9.0}  {:>7}  {:>7.2e}  {:>6.0}%  {:>7.1}",
                    d.rank,
                    truncate(&d.candidate.id, 42),
                    d.score,
                    sigma_str,
                    m.throat_shear_pa,
                    if m.fda_main_compliant { "  OK  " } else { "  !!  " },
                    m.hemolysis_index_per_pass,
                    m.well_coverage_fraction * 100.0,
                    m.total_pressure_drop_pa * 1e-3,
                );
            }

            println!("\n  ── Detailed breakdown ──");
            for d in &designs {
                println!();
                print_detailed(d);
            }
        }
        Err(e) => {
            eprintln!("  ERROR: {e}");
            eprintln!("  Hint: check that the default parameter sweep contains at least");
            eprintln!("  5 candidates with sigma < 1 AND main-channel shear < 150 Pa.");
            eprintln!("  Try increasing inlet gauge pressure or reducing throat diameter.");
        }
    }

    // ── Mode 2: Uniform Exposure ─────────────────────────────────────────
    println!("\n{}", "=".repeat(110));
    println!("  cfd-optim  |  Uniform Exposure Mode  |  Maximise light / ultrasound coverage of 6x6 well centre zone");
    println!("{}", "=".repeat(110));
    println!("  Objective: Maximise flow uniformity + well coverage fraction");
    println!("           + normalised residence time in the treatment zone");
    println!("  Constraint: same as above (dP feasibility + FDA main channel)");
    println!("{}", "-".repeat(110));

    match SdtOptimizer::new(OptimMode::UniformExposure, weights).top_5() {
        Ok(designs) => {
            println!(
                "{:>4}  {:<42}  {:>7}  {:>9}  {:>9}  {:>7}  {:>7}  {:>7}  {:>9}",
                "#",
                "Candidate ID",
                "Score",
                "Uniform",
                "τ_main Pa",
                "FDA",
                "Cov%",
                "t_res s",
                "Path mm"
            );
            println!("{}", "-".repeat(110));
            for d in &designs {
                let m = &d.metrics;
                println!(
                    "{:>4}  {:<42}  {:>7.4}  {:>9.4}  {:>9.1}  {:>7}  {:>6.0}%  {:>7.2}  {:>9.1}",
                    d.rank,
                    truncate(&d.candidate.id, 42),
                    d.score,
                    m.flow_uniformity,
                    m.max_main_channel_shear_pa,
                    if m.fda_main_compliant { "  OK  " } else { "  !!  " },
                    m.well_coverage_fraction * 100.0,
                    m.mean_residence_time_s,
                    m.total_path_length_mm,
                );
            }

            println!("\n  ── Detailed breakdown ──");
            for d in &designs {
                println!();
                print_detailed(d);
            }
        }
        Err(e) => {
            eprintln!("  ERROR: {e}");
        }
    }

    println!("\n{}", "=".repeat(110));
    println!("  Done.  Physical units: Pa, m3/s, m, s  |  Plate: ANSI/SLAS 1-2004 96-well");
    println!("{}", "=".repeat(110));
}

// ── Helpers ──────────────────────────────────────────────────────────────────

fn truncate(s: &str, n: usize) -> &str {
    if s.len() <= n {
        s
    } else {
        &s[..n]
    }
}

fn print_detailed(d: &cfd_optim::RankedDesign) {
    let c = &d.candidate;
    let m = &d.metrics;

    println!("  Rank #{}: {}", d.rank, c.id);
    println!("    Topology   : {}", c.topology.name());
    println!(
        "    Flow rate  : {:.2} mL/min  ({:.3e} m³/s)",
        c.flow_rate_m3_s * 6e7,
        c.flow_rate_m3_s
    );
    println!(
        "    Inlet gauge: {:.0} kPa  (abs: {:.1} kPa)",
        c.inlet_gauge_pa * 1e-3,
        c.inlet_pressure_pa() * 1e-3
    );
    if c.topology.has_venturi() {
        println!(
            "    Throat Ø   : {:.0} μm  (inlet Ø: {:.0} μm,  L_throat: {:.0} μm)",
            c.throat_diameter_m * 1e6,
            c.inlet_diameter_m * 1e6,
            c.throat_length_m * 1e6
        );
        println!(
            "    Cavitation : σ = {:.4}  |  potential = {:.4}  (cavitation if σ < 1)",
            m.cavitation_number, m.cavitation_potential
        );
        println!(
            "    Throat shear: {:.0} Pa  (throat shear rate: {:.0} 1/s)  FDA: {}",
            m.throat_shear_pa,
            m.throat_shear_rate_inv_s,
            if m.throat_exceeds_fda { "exceeds 150 Pa (brief transit)" } else { "within 150 Pa" }
        );
    }
    println!(
        "    Chan W × H : {:.0} × {:.0} μm  |  n_segs: {}  seg_len: {:.0} mm",
        c.channel_width_m * 1e6,
        c.channel_height_m * 1e6,
        c.serpentine_segments,
        c.segment_length_m * 1e3,
    );
    println!(
        "    Main shear : {:.1} Pa  (FDA ≤ 150 Pa: {})",
        m.max_main_channel_shear_pa,
        if m.fda_main_compliant { "PASS" } else { "FAIL" }
    );
    println!(
        "    Haemolysis : HI/pass = {:.3e}  (limit: {:.3e}  →  {})",
        m.hemolysis_index_per_pass,
        cfd_optim::constraints::HI_PASS_LIMIT,
        if m.hemolysis_index_per_pass <= cfd_optim::constraints::HI_PASS_LIMIT {
            "PASS"
        } else {
            "FAIL"
        }
    );
    println!(
        "    Coverage   : {:.0}%  |  Uniformity: {:.4}  |  Residence: {:.3} s",
        m.well_coverage_fraction * 100.0,
        m.flow_uniformity,
        m.mean_residence_time_s
    );
    println!(
        "    Pressure   : ΔP = {:.1} kPa  (gauge = {:.1} kPa  →  {})",
        m.total_pressure_drop_pa * 1e-3,
        c.inlet_gauge_pa * 1e-3,
        if m.pressure_feasible { "FEASIBLE" } else { "OVER LIMIT" }
    );
    println!(
        "    Path length: {:.1} mm  |  Score: {:.4}",
        m.total_path_length_mm, d.score
    );
}
