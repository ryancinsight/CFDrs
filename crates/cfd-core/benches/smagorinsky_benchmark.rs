use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::turbulence::{SmagorinskyModel, TurbulenceModel};
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_smagorinsky(c: &mut Criterion) {
    let nx = 64;
    let ny = 64;
    let nz = 64;
    let mut flow_field = FlowField::<f64>::new(nx, ny, nz);
    let delta = 1.0 / nx as f64;

    // Initialize with some non-zero velocity to make calculations meaningful
    // Using simple linear gradients to ensure non-zero strain rate
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                if let Some(vel) = flow_field.velocity.get_mut(i, j, k) {
                    vel.x = (i as f64) * delta;
                    vel.y = (j as f64) * delta;
                    vel.z = (k as f64) * delta;
                }
            }
        }
    }

    let model = SmagorinskyModel::new(0.1);

    c.bench_function("smagorinsky_viscosity_64^3", |b| {
        b.iter(|| {
            model.turbulent_viscosity(&flow_field)
        })
    });
}

criterion_group!(benches, bench_smagorinsky);
criterion_main!(benches);
