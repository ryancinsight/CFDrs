use cfd_3d::physics::turbulence::sigma::SigmaModel;
use cfd_3d::physics::turbulence::spalart_allmaras::SpalartAllmarasModel;
use cfd_3d::physics::turbulence::vreman::VremanModel;
use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;

#[test]
fn test_sigma_vanishing_solid_body_rotation() {
    // Sigma model must exactly vanish for solid body rotation
    let nx = 5;
    let ny = 5;
    let nz = 5;
    let mut flow_field = FlowField::<f64>::new(nx, ny, nz);
    let delta = 1.0;

    // Solid body rotation: u = -omega * y, v = omega * x, w = 0
    let omega = 2.0;
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64 * delta;
                let y = j as f64 * delta;

                if let Some(vel) = flow_field.velocity.get_mut(i, j, k) {
                    vel.x = -omega * y;
                    vel.y = omega * x;
                    vel.z = 0.0;
                }
            }
        }
    }

    let model = SigmaModel::with_filter_width(delta, delta, delta);
    let viscosity = model.turbulent_viscosity(&flow_field);

    // Check interior points
    for k in 1..nz - 1 {
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = k * nx * ny + j * nx + i;
                assert!(
                    viscosity[idx] < 1e-12,
                    "Sigma model MUST vanish for solid body rotation. Got {}",
                    viscosity[idx]
                );
            }
        }
    }
}

#[test]
fn test_vreman_vanishing_pure_shear() {
    // Vreman model must vanish for pure shear flow
    let nx = 5;
    let ny = 5;
    let nz = 5;
    let mut flow_field = FlowField::<f64>::new(nx, ny, nz);
    let delta = 1.0;

    // Pure shear: u = gamma * y, v = 0, w = 0
    let gamma = 3.0;
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let y = j as f64 * delta;

                if let Some(vel) = flow_field.velocity.get_mut(i, j, k) {
                    vel.x = gamma * y;
                    vel.y = 0.0;
                    vel.z = 0.0;
                }
            }
        }
    }

    let model = VremanModel::with_filter_width(delta, delta, delta);
    let viscosity = model.turbulent_viscosity(&flow_field);

    // Check interior points
    for k in 1..nz - 1 {
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = k * nx * ny + j * nx + i;
                assert!(
                    viscosity[idx] < 1e-12,
                    "Vreman model MUST vanish for pure shear flow. Got {:?}",
                    viscosity[idx]
                );
            }
        }
    }
}

#[test]
fn test_spalart_allmaras_math_bounds() {
    // SA model damping function testing
    // f_v1 = chi^3 / (chi^3 + c_v1^3)
    let nx = 2;
    let ny = 2;
    let nz = 2;
    let nu = 1e-6_f64;
    
    // Test wall damping limits implicitly
    let model = SpalartAllmarasModel::new(nx * ny * nz, nu, vec![1.0; nx * ny * nz]);
    let flow_field = FlowField::<f64>::new(nx, ny, nz);
    
    let viscosity = model.turbulent_viscosity(&flow_field);
    for v in viscosity {
        assert!(v >= 0.0, "Spalart Allmaras viscosity must be bounded positively");
    }
}
