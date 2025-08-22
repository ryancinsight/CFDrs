use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use cfd_core::domains::{
    fluid_dynamics::{FlowField, VelocityField, PressureField, FluidDynamicsService},
    numerical_methods::{DiscretizationScheme, TimeIntegrationScheme, finite_difference, time_integration},
};
use nalgebra::Vector3;
use std::collections::HashMap;

fn benchmark_flow_field_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("flow_field_operations");
    
    for size in [100, 1000, 10000].iter() {
        let flow_field = create_test_flow_field(*size);
        
        group.bench_with_input(
            BenchmarkId::new("divergence", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(flow_field.divergence())
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("vorticity", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(flow_field.vorticity())
                })
            },
        );
    }
    
    group.finish();
}

fn benchmark_numerical_schemes(c: &mut Criterion) {
    let mut group = c.benchmark_group("numerical_schemes");
    
    for size in [1000, 10000, 100000].iter() {
        let field: Vec<f64> = (0..*size).map(|i| (i as f64).sin()).collect();
        let dx = 0.01;
        
        let central_scheme = finite_difference::CentralDifference;
        let upwind_scheme = finite_difference::UpwindDifference;
        
        group.bench_with_input(
            BenchmarkId::new("central_difference", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(central_scheme.discretize(&field, dx))
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("upwind_difference", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(upwind_scheme.discretize(&field, dx))
                })
            },
        );
    }
    
    group.finish();
}

fn benchmark_time_integration(c: &mut Criterion) {
    let mut group = c.benchmark_group("time_integration");
    
    for size in [1000, 10000, 100000].iter() {
        let current: Vec<f64> = vec![1.0; *size];
        let derivative: Vec<f64> = vec![0.1; *size];
        let dt = 0.001;
        
        let forward_euler = time_integration::ForwardEuler;
        let rk4 = time_integration::RungeKutta4;
        
        group.bench_with_input(
            BenchmarkId::new("forward_euler", size),
            size,
            |b, _| {
                b.iter(|| {
                    let mut state = current.clone();
                    let _ = forward_euler.step(&mut state, 0.0, dt, |_, s| derivative(s));
                    black_box(state)
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("runge_kutta_4", size),
            size,
            |b, _| {
                b.iter(|| {
                    let mut state = current.clone();
                    let _ = rk4.step(&mut state, 0.0, dt, |_, s| derivative(s));
                    black_box(state)
                })
            },
        );
    }
    
    group.finish();
}

fn benchmark_mesh_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("mesh_operations");

    for size in [100, 500, 1000].iter() {
        // Create a simple test mesh for benchmarking
        let _vertices = create_test_vertices(*size);
        let _elements = create_test_elements(*size);

        group.bench_with_input(
            BenchmarkId::new("mesh_vertex_creation", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(create_test_vertices(*size))
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("mesh_element_creation", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(create_test_elements(*size))
                })
            },
        );
    }

    group.finish();
}

fn benchmark_reynolds_number_calculation(c: &mut Criterion) {
    let mut group = c.benchmark_group("reynolds_number");
    
    let service = FluidDynamicsService::<f64>::new();
    
    for count in [1000, 10000, 100000].iter() {
        group.bench_with_input(
            BenchmarkId::new("reynolds_calculation", count),
            count,
            |b, _| {
                b.iter(|| {
                    for i in 0..*count {
                        let velocity = 1.0 + (i as f64) * 0.001;
                        let length = 0.1;
                        let viscosity = 1e-6;
                        black_box(service.reynolds_number(velocity, length, viscosity));
                    }
                })
            },
        );
    }
    
    group.finish();
}

// Helper functions
fn create_test_flow_field(size: usize) -> FlowField<f64> {
    let velocity_components: Vec<Vector3<f64>> = (0..size)
        .map(|i| Vector3::new((i as f64).sin(), (i as f64).cos(), 0.0))
        .collect();
    
    let pressure_values: Vec<f64> = (0..size)
        .map(|i| 101325.0 + (i as f64) * 10.0)
        .collect();
    
    FlowField {
        velocity: VelocityField {
            components: velocity_components,
            dimensions: (size, 1, 1),
        },
        pressure: PressureField {
            values: pressure_values,
            dimensions: (size, 1, 1),
        },
        scalars: HashMap::new(),
    }
}

fn create_test_vertices(size: usize) -> Vec<nalgebra::Point3<f64>> {
    (0..size)
        .map(|i| {
            let x = (i as f64) / (size as f64);
            let y = ((i * 2) as f64) / (size as f64);
            let z = ((i * 3) as f64) / (size as f64);
            nalgebra::Point3::new(x, y, z)
        })
        .collect()
}

fn create_test_elements(size: usize) -> Vec<Vec<usize>> {
    (0..size/4)
        .map(|i| vec![i*4, i*4+1, i*4+2, i*4+3])
        .collect()
}

criterion_group!(
    benches,
    benchmark_flow_field_operations,
    benchmark_numerical_schemes,
    benchmark_time_integration,
    benchmark_mesh_operations,
    benchmark_reynolds_number_calculation
);
criterion_main!(benches);
