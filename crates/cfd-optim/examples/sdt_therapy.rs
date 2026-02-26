//! End-to-end SDT therapy pipeline.
//!
//! Runs two complementary searches over all 14 millifluidic topology families
//! and ranks designs by the combined [`OptimMode::SdtTherapy`] objective:
//!
//! - **35%** three-population separation efficiency (WBC + cancer → center,
//!   RBC → periphery)
//! - **30%** haemolysis minimisation (1 − HI / HI_limit)
//! - **20%** cavitation potential (σ < 1 → bubbles form)
//! - **15%** flow uniformity
//!
//! Hard constraints: total ΔP ≤ inlet gauge pressure AND main-channel shear ≤ 150 Pa.
//!
//! # Output
//!
//! - Ranked console tables (parametric top-5 and GA top-5)
//! - JSON result files in `cfd-optim/outputs/sdt_therapy/`
//! - SVG schematic for each top-5 design
//! - Head-to-head comparison: parametric rank-1 vs GA rank-1
//!
//! # Run with
//! ```bash
//! cargo run -p cfd-optim --example sdt_therapy
//! ```

use cfd_optim::{
    evo::GeneticOptimizer,
    save_comparison_svg, save_schematic_svg, save_top5_json,
    OptimMode, SdtOptimizer, SdtWeights,
};

fn main() {
    let out_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("outputs").join("sdt_therapy");
    std::fs::create_dir_all(&out_dir).expect("failed to create output directory");

    let weights = SdtWeights::default();
    let mode = OptimMode::SdtTherapy;

    print_header();

    // ── Part 1: Parametric sweep ──────────────────────────────────────────────
    println!("\n{}", "=".repeat(110));
    println!("  PART 1 — Parametric Sweep  |  Mode: {}", cfd_optim::score_description(mode));
    println!("{}", "=".repeat(110));
    println!("  Objectives: 35% 3-pop sep + 30% HI minimisation + 20% cavitation + 15% uniformity");
    println!("  Constraints: ΔP ≤ gauge pressure  |  main-channel shear ≤ 150 Pa");
    println!("{}", "-".repeat(110));

    let parametric_top5 = match SdtOptimizer::new(mode, weights).top_5() {
        Ok(d) => d,
        Err(e) => {
            eprintln!("  ERROR (parametric): {e}");
            return;
        }
    };

    print_table_header();
    for d in &parametric_top5 {
        print_row(d);
    }

    println!("\n  ── Detailed breakdown (parametric top-5) ──");
    for d in &parametric_top5 {
        println!();
        print_detailed(d);
    }

    // Save parametric results.
    let param_json = out_dir.join("parametric_top5.json");
    match save_top5_json(&parametric_top5, &param_json) {
        Ok(_)  => println!("\n  Saved: {}", param_json.display()),
        Err(e) => eprintln!("  WARN: could not save parametric JSON: {e}"),
    }

    // Save SVG schematics.
    for d in &parametric_top5 {
        let svg_path = out_dir.join(format!("param_{:02}_{}.svg", d.rank, safe_id(&d.candidate.id)));
        match save_schematic_svg(&d.candidate, &svg_path) {
            Ok(_)  => println!("  Saved: {}", svg_path.display()),
            Err(e) => eprintln!("  WARN: SVG failed for {}: {e}", d.candidate.id),
        }
    }

    // ── Part 2: Genetic algorithm ─────────────────────────────────────────────
    println!("\n{}", "=".repeat(110));
    println!("  PART 2 — Genetic Algorithm  |  Mode: {}", cfd_optim::score_description(mode));
    println!("{}", "=".repeat(110));
    println!("  Algorithm: real-coded GA  |  Population: 60  |  Generations: 100");
    println!("  Genome: 8 continuous genes (topology, Q, P_gauge, d_throat, w_ch, n_segs, L_seg, R_bend)");
    println!("{}", "-".repeat(110));

    let ga_optimizer = GeneticOptimizer::new(mode, weights)
        .with_population(60)
        .with_max_generations(100)
        .with_top_k(5);

    let ga_top5 = match ga_optimizer.run() {
        Ok(d) => d,
        Err(e) => {
            eprintln!("  ERROR (GA): {e}");
            return;
        }
    };

    print_table_header();
    for d in &ga_top5 {
        print_row(d);
    }

    println!("\n  ── Detailed breakdown (GA top-5) ──");
    for d in &ga_top5 {
        println!();
        print_detailed(d);
    }

    // Save GA results.
    let ga_json = out_dir.join("ga_top5.json");
    match save_top5_json(&ga_top5, &ga_json) {
        Ok(_)  => println!("\n  Saved: {}", ga_json.display()),
        Err(e) => eprintln!("  WARN: could not save GA JSON: {e}"),
    }

    // Save GA SVG schematics.
    for d in &ga_top5 {
        let svg_path = out_dir.join(format!("ga_{:02}_{}.svg", d.rank, safe_id(&d.candidate.id)));
        match save_schematic_svg(&d.candidate, &svg_path) {
            Ok(_)  => println!("  Saved: {}", svg_path.display()),
            Err(e) => eprintln!("  WARN: SVG failed for {}: {e}", d.candidate.id),
        }
    }

    // ── Part 3: Head-to-head comparison ──────────────────────────────────────
    println!("\n{}", "=".repeat(110));
    println!("  PART 3 — Head-to-Head: Parametric Rank #1 vs GA Rank #1");
    println!("{}", "=".repeat(110));

    if let (Some(p1), Some(g1)) = (parametric_top5.first(), ga_top5.first()) {
        let comparison = [p1.clone(), g1.clone()];
        let comp_svg = out_dir.join("head_to_head.svg");
        match save_comparison_svg(&comparison, &comp_svg, mode) {
            Ok(_)  => println!("  Comparison SVG saved: {}", comp_svg.display()),
            Err(e) => eprintln!("  WARN: comparison SVG failed: {e}"),
        }

        println!();
        println!("  {:.<50} {:^25} {:^25}", "", "Parametric #1", "GA #1");
        println!("  {:<50} {:>25} {:>25}", "Candidate ID",
            truncate(&p1.candidate.id, 24), truncate(&g1.candidate.id, 24));
        println!("  {:<50} {:>25} {:>25}", "Topology",
            p1.candidate.topology.name(), g1.candidate.topology.name());
        println!("  {:<50} {:>24.4} {:>24.4}", "SDT Therapy Score",
            p1.score, g1.score);

        let pm = &p1.metrics;
        let gm = &g1.metrics;

        println!("  {:<50} {:>24.4} {:>24.4}", "3-pop separation efficiency",
            pm.three_pop_sep_efficiency, gm.three_pop_sep_efficiency);
        println!("  {:<50} {:>24.2e} {:>24.2e}", "Haemolysis index / pass",
            pm.hemolysis_index_per_pass, gm.hemolysis_index_per_pass);
        let sigma_p = if pm.cavitation_number.is_finite() {
            format!("{:.4}", pm.cavitation_number)
        } else {
            "∞".to_string()
        };
        let sigma_g = if gm.cavitation_number.is_finite() {
            format!("{:.4}", gm.cavitation_number)
        } else {
            "∞".to_string()
        };
        println!("  {:<50} {:>25} {:>25}", "Cavitation number σ", sigma_p, sigma_g);
        println!("  {:<50} {:>24.4} {:>24.4}", "Flow uniformity",
            pm.flow_uniformity, gm.flow_uniformity);
        println!("  {:<50} {:>23.1} Pa {:>23.1} Pa", "Total pressure drop",
            pm.total_pressure_drop_pa, gm.total_pressure_drop_pa);
        println!("  {:<50} {:>24.0}% {:>24.0}%", "Well coverage",
            pm.well_coverage_fraction * 100.0, gm.well_coverage_fraction * 100.0);
        println!("  {:<50} {:>25} {:>25}",
            "FDA main-channel compliant",
            if pm.fda_main_compliant { "YES" } else { "NO" },
            if gm.fda_main_compliant { "YES" } else { "NO" });
        println!("  {:<50} {:>25} {:>25}",
            "Pressure feasible",
            if pm.pressure_feasible { "YES" } else { "NO" },
            if gm.pressure_feasible { "YES" } else { "NO" });
    }

    println!("\n{}", "=".repeat(110));
    println!("  Done.  Results in: {}", out_dir.display());
    println!("  Units: Pa, m³/s, m, s  |  Plate: ANSI/SLAS 1-2004 96-well (127.76 × 85.47 mm)");
    println!("{}", "=".repeat(110));
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn print_header() {
    println!("\n{}", "=".repeat(110));
    println!("  cfd-optim  |  SDT Therapy Pipeline  |  96-well plate, blood, FDA ≤ 150 Pa");
    println!("{}", "=".repeat(110));
    println!("  Combined therapy scoring: 35% 3-pop separation + 30% HI + 20% cavitation + 15% uniformity");
    println!("  Topology families: SingleVenturi … TrifurcationSerpentine (14 fixed + AdaptiveTree GA)");
    println!("{}", "=".repeat(110));
}

fn print_table_header() {
    println!(
        "  {:>4}  {:<38}  {:>7}  {:>8}  {:>8}  {:>7}  {:>7}  {:>7}  {:>7}",
        "#", "Candidate ID", "Score", "Sep3Eff", "HI/pass", "σ", "Unif", "Cov%", "ΔP kPa"
    );
    println!("  {}", "-".repeat(100));
}

fn print_row(d: &cfd_optim::RankedDesign) {
    let m = &d.metrics;
    let sigma = if m.cavitation_number.is_finite() {
        format!("{:>7.3}", m.cavitation_number)
    } else {
        format!("{:>7}", "∞")
    };
    println!(
        "  {:>4}  {:<38}  {:>7.4}  {:>8.4}  {:>8.2e}  {}  {:>7.4}  {:>6.0}%  {:>7.1}",
        d.rank,
        truncate(&d.candidate.id, 38),
        d.score,
        m.three_pop_sep_efficiency,
        m.hemolysis_index_per_pass,
        sigma,
        m.flow_uniformity,
        m.well_coverage_fraction * 100.0,
        m.total_pressure_drop_pa * 1e-3,
    );
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
            "    Throat     : width {:.0} μm  inlet {:.0} μm  L {:.0} μm",
            c.throat_diameter_m * 1e6,
            c.inlet_diameter_m * 1e6,
            c.throat_length_m * 1e6
        );
        let sigma_str = if m.cavitation_number.is_finite() {
            format!("{:.4}", m.cavitation_number)
        } else {
            "∞".to_string()
        };
        println!(
            "    Cavitation : σ = {}  |  potential = {:.4}",
            sigma_str, m.cavitation_potential
        );
    }
    println!(
        "    Main shear : {:.1} Pa  (FDA ≤ 150 Pa: {})",
        m.max_main_channel_shear_pa,
        if m.fda_main_compliant { "PASS" } else { "FAIL" }
    );
    println!(
        "    Haemolysis : HI/pass = {:.3e}  (limit {:.3e}  →  {})",
        m.hemolysis_index_per_pass,
        cfd_optim::constraints::HI_PASS_LIMIT,
        if m.hemolysis_index_per_pass <= cfd_optim::constraints::HI_PASS_LIMIT {
            "PASS"
        } else {
            "FAIL"
        }
    );
    println!(
        "    3-pop sep  : eff = {:.4}  WBC_center = {:.1}%  RBC_periph = {:.1}%",
        m.three_pop_sep_efficiency,
        m.wbc_center_fraction * 100.0,
        m.rbc_peripheral_fraction_three_pop * 100.0
    );
    println!(
        "    Cell eq.   : cancer x̃ = {:.3}  WBC x̃ = {:.3}  RBC x̃ = {:.3}",
        m.cancer_equilibrium_pos, m.wbc_equilibrium_pos, m.rbc_equilibrium_pos
    );
    println!(
        "    Coverage   : {:.0}%  |  Uniformity: {:.4}  |  Residence: {:.3} s",
        m.well_coverage_fraction * 100.0,
        m.flow_uniformity,
        m.mean_residence_time_s
    );
    println!(
        "    Pressure   : ΔP = {:.1} kPa  (gauge {:.1} kPa  →  {})",
        m.total_pressure_drop_pa * 1e-3,
        c.inlet_gauge_pa * 1e-3,
        if m.pressure_feasible { "FEASIBLE" } else { "OVER LIMIT" }
    );
    println!(
        "    Path length: {:.1} mm  |  SDT Therapy Score: {:.4}",
        m.total_path_length_mm, d.score
    );
}

fn truncate(s: &str, n: usize) -> &str {
    if s.len() <= n { s } else { &s[..n] }
}

fn safe_id(id: &str) -> String {
    id.chars()
        .map(|c| if c.is_alphanumeric() || c == '_' || c == '-' { c } else { '_' })
        .collect()
}
