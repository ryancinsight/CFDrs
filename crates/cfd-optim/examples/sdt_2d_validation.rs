//! 1D-vs-2D venturi throat validation for top SDT designs.
//!
//! Runs the parametric optimizer in [`OptimMode::SdtCavitation`] mode to obtain
//! the top 5 designs, then re-evaluates each design's venturi throat using the
//! 2D Navier-Stokes FVM solver ([`VenturiSolver2D`]).
//!
//! For each design that contains a venturi throat the example prints a
//! side-by-side comparison of:
//!
//! | Metric | 1D source | 2D source |
//! |--------|-----------|-----------|
//! | sigma (cavitation number) | `SdtMetrics::cavitation_number` | Bernoulli from 2D solution |
//! | dp_throat \[Pa\] | Bernoulli from candidate geometry | `VenturiFlowSolution::dp_throat` |
//!
//! # Run
//!
//! ```bash
//! cargo run -p cfd-optim --example sdt_2d_validation
//! ```

use cfd_2d::solvers::ns_fvm::BloodModel;
use cfd_2d::solvers::venturi_flow::{VenturiGeometry, VenturiSolver2D};
use cfd_optim::{OptimMode, SdtOptimizer, SdtWeights};

/// Atmospheric pressure [Pa].
const P_ATM_PA: f64 = 101_325.0;

/// Blood vapour pressure at 37 C [Pa].
const P_VAPOR_PA: f64 = 6_280.0;

/// Blood (water) density [kg/m^3].
const RHO: f64 = 1_060.0;

/// Water-equivalent dynamic viscosity for Newtonian 2D comparison [Pa s].
const MU_WATER: f64 = 8.9e-4;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // ── 1. Obtain top-5 SDT cavitation designs from the 1D optimizer ─────
    let optimizer = SdtOptimizer::new(OptimMode::SdtCavitation, SdtWeights::default());
    let top5 = optimizer.top_k(5)?;

    println!("=== 1D vs 2D Venturi Throat Comparison ===\n");
    println!(
        "{:>4}  {:<42} {:>10} {:>10} {:>10} {:>10}",
        "Rank", "ID", "sigma_1D", "sigma_2D", "dp_1D", "dp_2D"
    );
    println!("{}", "-".repeat(92));

    let mut validated = 0_usize;

    for design in &top5 {
        let c = &design.candidate;
        let m = &design.metrics;

        // Skip designs without a venturi throat.
        if !c.topology.has_venturi() {
            println!("{:>4}  {:<42}  (no venturi — skipped)", design.rank, c.id);
            continue;
        }

        // ── 2. Compute 1D Bernoulli dp_throat from candidate geometry ────
        //
        // The candidate uses rectangular venturi cross-sections (width × height).
        // The per-venturi flow rate accounts for selective-separation topologies
        // that split flow across multiple branches.
        let q_venturi = c.per_venturi_flow();
        let a_inlet = c.inlet_area_m2();
        let a_throat = c.throat_area_m2();
        let v_inlet_1d = q_venturi / a_inlet.max(1e-18);
        let v_throat_1d = q_venturi / a_throat.max(1e-18);
        let dp_throat_1d = 0.5 * RHO * (v_throat_1d * v_throat_1d - v_inlet_1d * v_inlet_1d);

        let sigma_1d = m.cavitation_number;

        // ── 3. Set up and solve the 2D venturi problem ───────────────────
        //
        // VenturiGeometry is a 2D (x, y) domain.  Channel widths map directly
        // from the candidate's rectangular cross-section dimensions.
        let geom = VenturiGeometry::<f64>::new(
            c.inlet_diameter_m,  // w_inlet
            c.throat_diameter_m, // w_throat
            0.003,               // l_inlet  (3 mm fixed approach)
            0.002,               // l_converge (2 mm linear taper)
            c.throat_length_m,   // l_throat
            0.004,               // l_diverge (4 mm recovery)
            c.channel_height_m,  // height
        );

        let blood = BloodModel::Newtonian(MU_WATER);
        // Adaptive grid: stretch y-cells toward centre to resolve throat.
        let cr = c.inlet_diameter_m / c.throat_diameter_m.max(1e-12);
        let ny_2d = (4.0 * cr).round().clamp(40.0, 200.0) as usize;
        let beta_2d =
            (1.0 - 4.0 * c.throat_diameter_m / c.inlet_diameter_m.max(1e-12)).clamp(0.0, 0.9);
        let mut solver = VenturiSolver2D::new_stretched(geom, blood, RHO, 60, ny_2d, beta_2d);

        // Inlet velocity for the 2D solver uses the rectangular cross-section.
        let area_2d = c.inlet_diameter_m * c.channel_height_m;
        let u_inlet_2d = q_venturi / area_2d.max(1e-18);

        let sol = match solver.solve(u_inlet_2d) {
            Ok(s) => s,
            Err(e) => {
                println!("{:>4}  {:<42}  2D solve failed: {}", design.rank, c.id, e);
                continue;
            }
        };

        // ── 4. Derive sigma from the 2D solution ────────────────────────
        //
        // sigma = (P_abs_throat - P_vapor) / (0.5 * rho * v_throat^2)
        //
        // The 2D solver reports gauge pressures relative to its own
        // reference.  We anchor the inlet to the candidate's absolute
        // pressure so that sigma is comparable to the 1D value.
        let p_abs_inlet = P_ATM_PA + c.inlet_gauge_pa;
        let p_abs_throat = p_abs_inlet + sol.dp_throat; // dp_throat is negative
        let dyn_p_2d = 0.5 * RHO * sol.u_throat * sol.u_throat;
        let sigma_2d = if dyn_p_2d > 1e-12 {
            (p_abs_throat - P_VAPOR_PA) / dyn_p_2d
        } else {
            f64::INFINITY
        };

        let dp_2d = sol.dp_throat.abs();

        println!(
            "{:>4}  {:<42} {:>10.4} {:>10.4} {:>10.1} {:>10.1}",
            design.rank, c.id, sigma_1d, sigma_2d, dp_throat_1d, dp_2d,
        );

        validated += 1;
    }

    println!(
        "\n{validated} / {} designs validated against 2D solver.",
        top5.len()
    );

    Ok(())
}
