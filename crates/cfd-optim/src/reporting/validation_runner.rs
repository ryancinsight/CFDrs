//! Shared Milestone 12 multi-fidelity venturi validation.

use std::path::Path;

use cfd_2d::solvers::ns_fvm::BloodModel;
use cfd_2d::solvers::venturi_flow::{VenturiGeometry, VenturiSolver2D};
use cfd_3d::venturi::{VenturiConfig3D, VenturiSolver3D};
use cfd_core::physics::fluid::blood::CarreauYasudaBlood;
use cfd_core::physics::fluid::CassonBlood;
use cfd_mesh::VenturiMeshBuilder;
use serde::Serialize;

use crate::constraints::{
    BLOOD_DENSITY_KG_M3 as RHO, BLOOD_VAPOR_PRESSURE_PA as P_VAPOR_PA, P_ATM_PA,
};
use crate::reporting::{pct_diff, Milestone12ReportDesign, ValidationRow};

#[derive(Debug, Serialize)]
struct Venturi2DResult {
    u_inlet_m_s: f64,
    u_throat_m_s: f64,
    dp_throat_pa: f64,
    dp_recovery_pa: f64,
    sigma_2d: f64,
    dp_bernoulli_1d_pa: f64,
}

#[derive(Debug, Serialize)]
struct Venturi3DResult {
    u_inlet_m_s: f64,
    u_throat_m_s: f64,
    dp_throat_pa: f64,
    dp_recovery_pa: f64,
    mass_error: f64,
    resolution: (usize, usize),
}

fn validate_venturi_candidate(
    track: &str,
    design: &Milestone12ReportDesign,
    out_dir: &Path,
) -> Result<ValidationRow, Box<dyn std::error::Error>> {
    let candidate = &design.candidate;
    let metrics = &design.metrics;
    let topology = candidate
        .topology_spec()
        .map_err(|error| error.to_string())?;
    let venturi = topology
        .venturi_placements
        .first()
        .ok_or_else(|| format!("validation requires venturi topology: {}", candidate.id))?;

    let parallel_paths = (metrics.active_venturi_throat_count
        / metrics.serial_venturi_stages_per_path.max(1))
    .max(1) as f64;
    let q = candidate.operating_point.flow_rate_m3_s * metrics.therapy_channel_fraction
        / parallel_paths;
    let inlet_width_m = venturi.throat_geometry.inlet_width_m;
    let throat_width_m = venturi.throat_geometry.throat_width_m;
    let channel_height_m = venturi.throat_geometry.throat_height_m;
    let throat_length_m = venturi.throat_geometry.throat_length_m;
    let a_in = inlet_width_m * channel_height_m;
    let a_th = throat_width_m * channel_height_m;
    let v_in_1d = q / a_in.max(1e-30);
    let v_th_1d = q / a_th.max(1e-30);
    let dp_bernoulli = 0.5 * RHO * (v_th_1d * v_th_1d - v_in_1d * v_in_1d);

    let geom2d = VenturiGeometry::<f64>::new(
        inlet_width_m,
        throat_width_m,
        3e-3,
        2e-3,
        throat_length_m,
        4e-3,
        channel_height_m,
    );
    let blood_2d = BloodModel::Casson(CassonBlood::<f64>::normal_blood());
    let cr = inlet_width_m / throat_width_m.max(1e-12);
    let ny_2d = (4.0 * cr).round().clamp(40.0, 200.0) as usize;
    let beta_2d = (1.0 - 4.0 * throat_width_m / inlet_width_m.max(1e-12)).clamp(0.0, 0.9);
    let mut solver2d = VenturiSolver2D::new_stretched(geom2d, blood_2d, RHO, 60, ny_2d, beta_2d);
    let u_inlet_2d = q / a_in.max(1e-30);
    let sol2d = solver2d
        .solve(u_inlet_2d)
        .map_err(|error| format!("2D FVM failed for {}: {error}", candidate.id))?;
    let p_abs_inlet = P_ATM_PA + candidate.operating_point.inlet_gauge_pa;
    let p_abs_throat = p_abs_inlet + sol2d.dp_throat;
    let dyn_p = 0.5 * RHO * sol2d.u_throat * sol2d.u_throat;
    let sigma_2d = if dyn_p > 1e-12 {
        (p_abs_throat - P_VAPOR_PA) / dyn_p
    } else {
        f64::INFINITY
    };
    let dp_2d = -sol2d.dp_throat;

    let result_2d = Venturi2DResult {
        u_inlet_m_s: u_inlet_2d,
        u_throat_m_s: sol2d.u_throat,
        dp_throat_pa: dp_2d,
        dp_recovery_pa: -sol2d.dp_recovery,
        sigma_2d,
        dp_bernoulli_1d_pa: dp_bernoulli,
    };
    std::fs::write(
        out_dir.join(format!("{}_2d_venturi.json", candidate.id)),
        serde_json::to_string_pretty(&result_2d)?,
    )?;

    let resolution = (60_usize, 10_usize);
    let builder3d = VenturiMeshBuilder::<f64>::new(
        inlet_width_m,
        throat_width_m,
        5.0 * inlet_width_m,
        3.0 * inlet_width_m,
        throat_length_m,
        7.0 * inlet_width_m,
        5.0 * inlet_width_m,
    )
    .with_resolution(resolution.0, resolution.1)
    .with_circular(false);

    let config3d = VenturiConfig3D::<f64> {
        inlet_flow_rate: q,
        resolution,
        circular: false,
        rect_height: Some(channel_height_m),
        ..Default::default()
    };
    let sol3d = VenturiSolver3D::new(builder3d, config3d)
        .solve(CarreauYasudaBlood::<f64>::normal_blood())
        .map_err(|error| format!("3D FEM failed for {}: {error}", candidate.id))?;
    let dp_3d = sol3d.dp_throat.abs();
    let mass_err_3d = sol3d.mass_error.abs();

    let result_3d = Venturi3DResult {
        u_inlet_m_s: sol3d.u_inlet,
        u_throat_m_s: sol3d.u_throat,
        dp_throat_pa: dp_3d,
        dp_recovery_pa: sol3d.dp_recovery.abs(),
        mass_error: mass_err_3d,
        resolution,
    };
    std::fs::write(
        out_dir.join(format!("{}_3d_result.json", candidate.id)),
        serde_json::to_string_pretty(&result_3d)?,
    )?;

    let row = ValidationRow {
        track: track.to_string(),
        id: candidate.id.clone(),
        topology: topology.short_code(),
        dp_1d_bernoulli_pa: dp_bernoulli,
        dp_2d_fvm_pa: dp_2d,
        dp_3d_fem_pa: dp_3d,
        agreement_1d_2d_pct: pct_diff(dp_bernoulli, dp_2d),
        agreement_2d_3d_pct: pct_diff(dp_2d, dp_3d),
        mass_error_3d_pct: mass_err_3d * 100.0,
        sigma_1d: metrics.cavitation_number,
        sigma_2d,
        score: design.score,
    };
    std::fs::write(
        out_dir.join(format!("{}_validation.json", candidate.id)),
        serde_json::to_string_pretty(&row)?,
    )?;
    Ok(row)
}

/// Run the canonical Milestone 12 venturi validation on the selected Option 2
/// and GA-ranked designs, writing per-design artifacts and the aggregate row set.
///
/// # Errors
/// Returns an error if the underlying 2D/3D solvers fail or any artifact cannot
/// be written.
pub fn run_milestone12_validation(
    out_dir: &Path,
    option2: &Milestone12ReportDesign,
    ga: &Milestone12ReportDesign,
) -> Result<Vec<ValidationRow>, Box<dyn std::error::Error>> {
    std::fs::create_dir_all(out_dir)?;
    let rows = vec![
        validate_venturi_candidate("Option 2 Combined Selective Venturi", option2, out_dir)?,
        validate_venturi_candidate("GA HydroSDT Venturi", ga, out_dir)?,
    ];
    std::fs::write(
        out_dir.join("milestone12_validation_rows.json"),
        serde_json::to_string_pretty(&rows)?,
    )?;
    Ok(rows)
}
