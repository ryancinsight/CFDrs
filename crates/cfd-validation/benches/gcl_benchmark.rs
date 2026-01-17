use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nalgebra::DMatrix;

fn benchmark_gcl_loop(c: &mut Criterion) {
    let mut group = c.benchmark_group("gcl_loop");
    // Use a large enough size to make allocation cost significant
    // 500x500 f64 = 250,000 * 8 bytes = 2 MB
    let nx = 500;
    let ny = 500;
    let steps = 10;
    let constant_value = 1.0;

    group.bench_function("baseline", |b| {
        b.iter(|| {
            // Setup similar to the test function
            let u = DMatrix::from_element(nx, ny, constant_value);
            let mut u_current = u.clone();

            for _ in 0..steps {
                // Baseline: Allocate new matrix every step
                let u_next = u_current.clone();

                // Simulate some access to prevent optimization
                black_box(&u_next[(0,0)]);

                u_current = u_next;
            }
        })
    });

    group.bench_function("optimized_copy_from", |b| {
        b.iter(|| {
            // Setup
            let u = DMatrix::from_element(nx, ny, constant_value);
            let mut u_current = u.clone();
            let mut u_next = u.clone(); // Pre-allocate

            for _ in 0..steps {
                // Optimized: Copy into existing buffer
                u_next.copy_from(&u_current);

                // Simulate some access
                black_box(&u_next[(0,0)]);

                // Swap buffers
                std::mem::swap(&mut u_current, &mut u_next);
            }
        })
    });

    group.bench_function("optimized_clone_from", |b| {
        b.iter(|| {
            // Setup
            let u = DMatrix::from_element(nx, ny, constant_value);
            let mut u_current = u.clone();
            let mut u_next = u.clone(); // Pre-allocate

            for _ in 0..steps {
                // Optimized: clone_from (should reuse capacity)
                u_next.clone_from(&u_current);

                // Simulate some access
                black_box(&u_next[(0,0)]);

                // Swap buffers
                std::mem::swap(&mut u_current, &mut u_next);
            }
        })
    });

    group.finish();
}

criterion_group!(benches, benchmark_gcl_loop);
criterion_main!(benches);
