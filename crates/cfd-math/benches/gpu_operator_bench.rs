use criterion::{criterion_group, criterion_main, Criterion};

#[cfg(feature = "gpu")]
fn bench_gpu_operator(c: &mut Criterion) {
    use std::sync::Arc;
    use cfd_core::compute::gpu::GpuContext;
    use cfd_math::linear_solver::operators::gpu::{GpuLaplacianOperator2D, BoundaryType};
    use cfd_math::linear_solver::LinearOperator;
    use nalgebra::DVector;

    // Initialize GPU context once
    // Note: This might fail if run in parallel or repeatedly if context doesn't clean up well
    let context = match GpuContext::create() {
        Ok(ctx) => Arc::new(ctx),
        Err(e) => {
            eprintln!("Skipping benchmark: Failed to create GpuContext: {:?}", e);
            return;
        }
    };

    let mut group = c.benchmark_group("gpu_operator");

    let nx = 128;
    let ny = 128;
    let dx = 1.0;
    let dy = 1.0;
    let bc = BoundaryType::Dirichlet;

    let operator = GpuLaplacianOperator2D::new(
        context.clone(),
        nx,
        ny,
        dx,
        dy,
        bc
    );

    let size = nx * ny;
    let x = DVector::from_element(size, 1.0f32);
    let mut y = DVector::zeros(size);

    group.bench_function("apply_gpu_op", |b| {
        b.iter(|| {
            operator.apply(&x, &mut y).unwrap();
        });
    });

    group.finish();
}

#[cfg(not(feature = "gpu"))]
fn bench_gpu_operator(_c: &mut Criterion) {
    // No-op
}

criterion_group!(benches, bench_gpu_operator);
criterion_main!(benches);
