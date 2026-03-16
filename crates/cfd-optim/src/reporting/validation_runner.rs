//! Shared Milestone 12 multi-fidelity venturi validation.

use std::path::Path;

use cfd_2d::solvers::ns_fvm::{BloodModel, SIMPLEConfig};
use cfd_2d::solvers::venturi_flow::{VenturiGeometry, VenturiSolver2D};
use cfd_3d::venturi::{VenturiConfig3D, VenturiSolver3D};
use cfd_core::physics::fluid::blood::CarreauYasudaBlood;
use cfd_core::physics::fluid::CassonBlood;
use cfd_mesh::VenturiMeshBuilder;
use serde::Serialize;

use crate::constraints::{
    BLOOD_DENSITY_KG_M3 as RHO, BLOOD_VAPOR_PRESSURE_PA as P_VAPOR_PA,
    BLOOD_VISCOSITY_PA_S as MU, P_ATM_PA,
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
    converged: bool,
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
    // Hydraulic diameter for rectangular cross-section: Dh = 2wh/(w+h).
    // Used for both the Re-adaptive SIMPLE config and the Stokes-validity flag.
    let dh_throat_m =
        2.0 * throat_width_m * channel_height_m / (throat_width_m + channel_height_m).max(1e-30);

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
    let u_inlet_2d = q / a_in.max(1e-30);

    // Pre-compute Re_throat from 1D kinematics so we can pick SIMPLE settings
    // before constructing the solver.  At Re > 100 the default α_u=0.7/α_p=0.3
    // diverges; use conservative relaxation with more iterations instead.
    // Use hydraulic diameter Dh (not raw width) for rectangular channels.
    let re_throat_1d = RHO * v_th_1d * dh_throat_m / MU;
    let simple_config_2d = if re_throat_1d > 100.0 {
        SIMPLEConfig::new(5000, 1e-5_f64, 0.3_f64, 0.15_f64, 1.0_f64, 1)
    } else {
        SIMPLEConfig::default()
    };
    let mut solver2d = VenturiSolver2D::new_stretched_with_config(
        geom2d, blood_2d, RHO, 60, ny_2d, beta_2d, simple_config_2d,
    );
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

    // Velocity-continuity check: if less than 50% of the expected throat velocity
    // is captured, the SIMPLE solver did not converge properly (e.g. Re too high).
    let u_throat_expected = u_inlet_2d * inlet_width_m / throat_width_m.max(1e-12);
    let two_d_converged = sol2d.converged && sol2d.u_throat >= 0.5 * u_throat_expected;

    let result_2d = Venturi2DResult {
        u_inlet_m_s: u_inlet_2d,
        u_throat_m_s: sol2d.u_throat,
        dp_throat_pa: dp_2d,
        dp_recovery_pa: -sol2d.dp_recovery,
        sigma_2d,
        dp_bernoulli_1d_pa: dp_bernoulli,
        converged: two_d_converged,
    };
    std::fs::write(
        out_dir.join(format!("{}_2d_venturi.json", candidate.id)),
        serde_json::to_string_pretty(&result_2d)?,
    )?;

    // Resolution pyramid: try canonical (80,20) first, fall back to (40,10)
    // on solver failure.  Coarser grid has 8× fewer elements and DOF, so
    // the iterative solver converges more reliably on stiff geometries.
    let resolutions: [(usize, usize); 2] = [(80, 20), (40, 10)];
    let mut sol3d = None;
    let mut resolution = resolutions[0];
    for &res in &resolutions {
        resolution = res;
        let builder3d = VenturiMeshBuilder::<f64>::new(
            inlet_width_m,
            throat_width_m,
            5.0 * inlet_width_m,
            3.0 * inlet_width_m,
            throat_length_m,
            7.0 * inlet_width_m,
            5.0 * inlet_width_m,
        )
        .with_resolution(res.0, res.1)
        .with_circular(false);

        let config3d = VenturiConfig3D::<f64> {
            inlet_flow_rate: q,
            resolution: res,
            circular: false,
            rect_height: Some(channel_height_m),
            ..Default::default()
        };
        match VenturiSolver3D::new(builder3d, config3d)
            .solve(CarreauYasudaBlood::<f64>::normal_blood())
        {
            Ok(s) => {
                sol3d = Some(s);
                break;
            }
            Err(e) => {
                tracing::warn!(
                    resolution = ?res,
                    error = %e,
                    "3D FEM failed at resolution {:?}, trying coarser", res
                );
            }
        }
    }
    let sol3d = sol3d
        .ok_or_else(|| format!("3D FEM failed for {} at all resolutions", candidate.id))?;
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

    // Re_throat = ρ × v_throat × Dh / μ  (Dh = 2wh/(w+h) for rectangular section).
    // The 3D Stokes solver is only physically valid for Re << 1 (creeping flow).
    // Flag inertia-dominated cases so the report can explain the expected 1D/3D discrepancy.
    let re_throat = RHO * v_th_1d * dh_throat_m / MU;
    let high_re_stokes_mismatch = re_throat > 1.0;

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
        two_d_converged,
        high_re_stokes_mismatch,
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
