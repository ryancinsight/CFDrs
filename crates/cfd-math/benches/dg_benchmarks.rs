//! Benchmarks for Discontinuous Galerkin methods

use cfd_math::high_order::dg::*;
use criterion::{criterion_group, criterion_main, Criterion};
use nalgebra::DVector;

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
                    let u0 = |x: f64| DVector::from_vec(vec![x]);
                    
                    let mut solver = DGSolver::new(
                        dg_op,
                        TimeIntegratorFactory::create(TimeIntegration::SSPRK3),
                        TimeIntegrationParams::new(TimeIntegration::SSPRK3)
                            .with_t_final(0.1)
                            .with_cfl(0.1),
                    ).unwrap();
                    
                    solver.initialize(u0).unwrap();
                    
                    b.iter(|| {
                        solver.step(|_, u| Ok(-u.clone())).unwrap();
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
                    let u0 = |x: f64| DVector::from_vec(vec![0.5 + x]);
                    
                    let mut solver = DGSolver::new(
                        dg_op,
                        TimeIntegratorFactory::create(TimeIntegration::SSPRK3),
                        TimeIntegrationParams::new(TimeIntegration::SSPRK3)
                            .with_t_final(0.1)
                            .with_cfl(0.1),
                    ).unwrap();
                    
                    solver.initialize(u0).unwrap();
                    
                    b.iter(|| {
                        solver.step(|_, u| {
                            let du_dx = solver.operator.compute_derivative(u)?;
                            Ok(-u.component_mul(&du_dx))
                        }).unwrap();
                    });
                },
            );
        }
    }
}

criterion_group!(
    benches,
    dg_advection_benchmark,
    dg_burgers_benchmark
);
criterion_main!(benches);
