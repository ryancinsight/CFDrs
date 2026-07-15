//! Benchmarks for Discontinuous Galerkin methods

use cfd_math::error::Result;
use cfd_math::high_order::dg::*;
use criterion::{criterion_group, criterion_main, Criterion};
use leto::{Array1, Array2};

fn dg_advection_benchmark(c: &mut Criterion) {
    let orders = [2, 4, 8];
    let num_elements_list = [16, 32, 64];

    for &order in &orders {
        for &num_elements in &num_elements_list {
            c.bench_function(
                &format!("DG Advection order={} elements={}", order, num_elements),
                |b| {
                    let params = DGOperatorParams::new()
                        .with_volume_flux(FluxType::Upwind)
                        .with_surface_flux(FluxType::Upwind);

                    let dg_op = DGOperator::new(order, 1, Some(params)).unwrap();
                    let u0 = |x: f64| Array1::from_shape_vec([1], vec![x]).unwrap();

                    let mut solver = DGSolver::new(
                        dg_op,
                        TimeIntegratorFactory::create(TimeIntegration::SSPRK3),
                        TimeIntegrationParams::new(TimeIntegration::SSPRK3)
                            .with_t_final(0.1)
                            .with_cfl(0.1),
                    );

                    solver.initialize(u0).unwrap();

                    b.iter(|| {
                        let f = |_: f64, u: &Array2<f64>| -> Result<Array2<f64>> {
                            Ok(Array2::from_shape_fn(u.shape(), |idx| -u[idx]))
                        };
                        solver
                            .step(&f, None::<&fn(f64, &Array2<f64>) -> Result<Array2<f64>>>)
                            .unwrap();
                    });
                },
            );
        }
    }
}

fn dg_burgers_benchmark(c: &mut Criterion) {
    let orders = [2, 3, 4];
    let num_elements_list = [16, 32, 64];

    for &order in &orders {
        for &num_elements in &num_elements_list {
            c.bench_function(
                &format!("DG Burgers order={} elements={}", order, num_elements),
                |b| {
                    let params = DGOperatorParams::new()
                        .with_volume_flux(FluxType::LaxFriedrichs)
                        .with_surface_flux(FluxType::LaxFriedrichs);

                    let dg_op = DGOperator::new(order, 1, Some(params)).unwrap();
                    let dg_op_clone = dg_op.clone();
                    let u0 = |x: f64| Array1::from_shape_vec([1], vec![0.5 + x]).unwrap();

                    let mut solver = DGSolver::new(
                        dg_op,
                        TimeIntegratorFactory::create(TimeIntegration::SSPRK3),
                        TimeIntegrationParams::new(TimeIntegration::SSPRK3)
                            .with_t_final(0.1)
                            .with_cfl(0.1),
                    );

                    solver.initialize(u0).unwrap();

                    b.iter(|| {
                        let f = |_: f64, u: &Array2<f64>| -> Result<Array2<f64>> {
                            let du_dx = dg_op_clone.compute_derivative(u)?;
                            Ok(Array2::from_shape_fn(u.shape(), |idx| -u[idx] * du_dx[idx]))
                        };
                        solver
                            .step(&f, None::<&fn(f64, &Array2<f64>) -> Result<Array2<f64>>>)
                            .unwrap();
                    });
                },
            );
        }
    }
}

criterion_group!(benches, dg_advection_benchmark, dg_burgers_benchmark);
criterion_main!(benches);
