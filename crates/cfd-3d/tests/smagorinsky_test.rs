use approx::assert_relative_eq;
use cfd_3d::physics::turbulence::SmagorinskyModel;
use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;
use nalgebra::Vector3;

#[test]
fn test_smagorinsky_turbulent_viscosity() {
    let nx = 4;
    let ny = 4;
    let nz = 4;
    let mut flow_field = FlowField::<f64>::new(nx, ny, nz);
    let delta = 1.0 / nx as f64;

    // Set linear velocity field: u=x, v=y, w=z
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64 * delta;
                let y = j as f64 * delta;
                let z = k as f64 * delta;

                if let Some(vel) = flow_field.velocity.get_mut(i, j, k) {
                    vel.x = x;
                    vel.y = y;
                    vel.z = z;
                }
            }
        }
    }

    let cs = 0.1;
    let model = SmagorinskyModel::with_filter_width(cs, delta, delta, delta);
    let viscosity = model.turbulent_viscosity(&flow_field);

    assert_eq!(viscosity.len(), nx * ny * nz);

    // Check internal point (1, 1, 1)
    // S11 = du/dx = 1
    // S22 = dv/dy = 1
    // S33 = dw/dz = 1
    // Off-diagonals should be 0 because u only depends on x, v only on y, w only on z
    // |S| = sqrt(2 * (S11^2 + S22^2 + S33^2)) = sqrt(2 * 3) = sqrt(6)
    // nu_t = (Cs * Delta)^2 * |S|

    // Strain rate calculation in model uses central difference, so it should be exact for linear field.
    // du/dx = ( (i+1)d - (i-1)d ) / (2d) = 2d / 2d = 1.

    let expected_strain_mag = (6.0f64).sqrt();
    let expected_viscosity = (cs * delta).powi(2) * expected_strain_mag;

    // Index for (1, 1, 1)
    let idx = 1 * nx * ny + 1 * nx + 1;
    let computed = viscosity[idx];

    println!("Computed: {}, Expected: {}", computed, expected_viscosity);
    assert!((computed - expected_viscosity).abs() < 1e-10);
}

#[test]
fn test_smagorinsky_anisotropic_filter_width() {
    let nx = 3;
    let ny = 3;
    let nz = 3;
    let dx = 1.0;
    let dy = 2.0;
    let dz = 4.0;
    let mut flow_field = FlowField::<f64>::new(nx, ny, nz);

    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64 * dx;
                let y = j as f64 * dy;
                let z = k as f64 * dz;
                if let Some(vel) = flow_field.velocity.get_mut(i, j, k) {
                    *vel = Vector3::new(x, y, z);
                }
            }
        }
    }

    let cs = 0.1;
    let model = SmagorinskyModel::with_filter_width(cs, dx, dy, dz);
    let viscosity = model.turbulent_viscosity(&flow_field);
    let delta = (dx * dy * dz).cbrt();
    let expected_strain_mag = (6.0f64).sqrt();
    let expected_viscosity = (cs * delta).powi(2) * expected_strain_mag;
    let center = 1 * nx * ny + 1 * nx + 1;

    assert_relative_eq!(viscosity[center], expected_viscosity, epsilon = 1e-12);
}
