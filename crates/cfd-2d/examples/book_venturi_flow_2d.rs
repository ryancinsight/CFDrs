//! Book example: 2D Venturi flow with a compact solver + Bernoulli comparison.

use cfd_2d::solvers::ns_fvm::BloodModel;
use cfd_2d::solvers::venturi_flow::{BernoulliVenturi, VenturiGeometry, VenturiSolver2D};
use std::fs;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let output_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("book_venturi_flow_2d");
    fs::create_dir_all(&output_dir)?;

    let geometry = VenturiGeometry::<f64>::iso_5167_standard();
    let inlet_velocity = 0.1_f64;
    let density = 1060.0_f64;

    let mut solver = VenturiSolver2D::new(
        geometry.clone(),
        BloodModel::Newtonian(1.0e-3),
        density,
        60,
        30,
    );
    let solution = solver.solve(inlet_velocity)?;

    let bernoulli = BernoulliVenturi::new(geometry, inlet_velocity, solution.p_inlet, density);
    let velocity_rel_err =
        (solution.u_throat - bernoulli.velocity_throat()).abs() / bernoulli.velocity_throat();

    assert!(
        solution.u_throat.is_finite(),
        "throat velocity must be finite"
    );
    assert!(
        solution.p_throat.is_finite(),
        "throat pressure must be finite"
    );
    assert!(
        solution.dp_throat > 0.0,
        "Venturi throat pressure drop must be positive"
    );
    assert!(
        velocity_rel_err < 0.2,
        "expected throat velocity to stay within 20% of Bernoulli, got {velocity_rel_err}"
    );

    let summary_path = output_dir.join("summary.txt");
    let summary = format!(
        "u_inlet={:.6}\nu_throat={:.6}\ndp_throat={:.6}\nvelocity_rel_err={:.6}\n",
        solution.u_inlet, solution.u_throat, solution.dp_throat, velocity_rel_err
    );
    fs::write(summary_path, summary)?;

    Ok(())
}
