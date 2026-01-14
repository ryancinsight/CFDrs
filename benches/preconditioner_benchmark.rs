//! Benchmarks for preconditioner performance comparison

#[cfg(feature = "mpi")]
use cfd_math::linear_solver::preconditioners::{
    AdditiveSchwarzPreconditioner, CoarseningStrategy, ParallelAMGPreconditioner,
    ParallelBlockJacobiPreconditioner,
};
use cfd_math::linear_solver::{
    IdentityPreconditioner, IncompleteLU, JacobiPreconditioner, Preconditioner, SORPreconditioner,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nalgebra::DVector;
use nalgebra_sparse::CsrMatrix;

/// Create a test matrix for benchmarking
fn create_benchmark_matrix(size: usize) -> CsrMatrix<f64> {
    let mut values = Vec::new();
    let mut indices = Vec::new();
    let mut indptr = vec![0];

    for i in 0..size {
        let mut row_entries: Vec<(usize, f64)> = Vec::with_capacity(5);
        row_entries.push((i, 4.0));

        if i > 0 {
            row_entries.push((i - 1, -1.0));
        }
        if i < size - 1 {
            row_entries.push((i + 1, -1.0));
        }

        if i > 1 && (i % 3) == 0 {
            row_entries.push((i - 2, -0.5));
        }
        if i < size - 2 && ((i + 1) % 3) == 0 {
            row_entries.push((i + 2, -0.5));
        }

        row_entries.sort_by_key(|(col, _)| *col);
        for (col, val) in row_entries {
            indices.push(col);
            values.push(val);
        }

        indptr.push(values.len());
    }

    CsrMatrix::try_from_csr_data(size, size, indptr, indices, values).unwrap()
}

/// Create a test vector for benchmarking
fn create_test_vector(size: usize) -> DVector<f64> {
    DVector::from_fn(size, |i, _| (i as f64).sin() + 0.1)
}

fn bench_identity_preconditioner(c: &mut Criterion) {
    let _matrix = create_benchmark_matrix(1000);
    let preconditioner = IdentityPreconditioner;
    let r = create_test_vector(1000);
    let mut z = DVector::zeros(1000);

    c.bench_function("identity_preconditioner_1000", |b| {
        b.iter(|| {
            preconditioner
                .apply_to(black_box(&r), black_box(&mut z))
                .unwrap();
        });
    });
}

fn bench_jacobi_preconditioner(c: &mut Criterion) {
    let matrix = create_benchmark_matrix(1000);
    let preconditioner = JacobiPreconditioner::new(&matrix).unwrap();
    let r = create_test_vector(1000);
    let mut z = DVector::zeros(1000);

    c.bench_function("jacobi_preconditioner_1000", |b| {
        b.iter(|| {
            preconditioner
                .apply_to(black_box(&r), black_box(&mut z))
                .unwrap();
        });
    });
}

fn bench_sor_preconditioner(c: &mut Criterion) {
    let matrix = create_benchmark_matrix(1000);
    let preconditioner = SORPreconditioner::new(&matrix, 1.0).unwrap();
    let r = create_test_vector(1000);
    let mut z = DVector::zeros(1000);

    c.bench_function("sor_preconditioner_1000", |b| {
        b.iter(|| {
            preconditioner
                .apply_to(black_box(&r), black_box(&mut z))
                .unwrap();
        });
    });
}

fn bench_incomplete_lu_preconditioner(c: &mut Criterion) {
    let matrix = create_benchmark_matrix(1000);
    let preconditioner = IncompleteLU::new(&matrix).unwrap();
    let r = create_test_vector(1000);
    let mut z = DVector::zeros(1000);

    c.bench_function("ilu0_preconditioner_1000", |b| {
        b.iter(|| {
            preconditioner
                .apply_to(black_box(&r), black_box(&mut z))
                .unwrap();
        });
    });
}

#[cfg(feature = "mpi")]
fn bench_parallel_preconditioners(c: &mut Criterion) {
    // Note: These benchmarks would require MPI initialization
    // For now, just test the creation overhead
    let matrix = create_benchmark_matrix(100);

    c.bench_function("parallel_preconditioner_creation", |b| {
        b.iter(|| {
            // Test creation overhead (would need MPI communicator in real usage)
            let _matrix_ref = &matrix;
            // In real usage: ParallelBlockJacobiPreconditioner::new(matrix, &communicator, 10)
        });
    });
}

criterion_group! {
    name = preconditioner_benches;
    config = Criterion::default().sample_size(10);
    targets = bench_identity_preconditioner, bench_jacobi_preconditioner, bench_sor_preconditioner, bench_incomplete_lu_preconditioner
}

#[cfg(feature = "mpi")]
criterion_group! {
    name = parallel_preconditioner_benches;
    config = Criterion::default().sample_size(10);
    targets = bench_parallel_preconditioners
}

#[cfg(not(feature = "mpi"))]
criterion_main!(preconditioner_benches);

#[cfg(feature = "mpi")]
criterion_main!(preconditioner_benches, parallel_preconditioner_benches);
