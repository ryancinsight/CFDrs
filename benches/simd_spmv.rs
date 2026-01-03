//! Sparse Matrix-Vector Multiply (SpMV) Benchmarks
//!
//! These benchmarks compare different SpMV implementations:
//! - **Scalar**: Baseline single-threaded implementation
//! - **SIMD** (Sprint 1.41.0): AVX2/SSE4.1 vectorization - **DEPRECATED** due to 27-32% regression
//! - **Parallel** (Sprint 1.67.0): Rayon-based parallelization - **RECOMMENDED**
//!
//! Performance findings (Sprint 1.55.0 validation):
//! - SIMD: 27-32% SLOWER than scalar (irregular memory access pattern)
//! - Parallel: 3-8x speedup on 4-8 cores (embarrassingly parallel rows)
//!
//! Test matrices:
//! - Tridiagonal (typical 1D diffusion, 3 nnz/row)
//! - Pentadiagonal (typical 2D Laplacian, 5 nnz/row)
//! - Large matrices (>1000 rows) for parallel benefit

use cfd_math::sparse::{spmv, spmv_parallel, SparseMatrixBuilder};
use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use nalgebra::DVector;
use nalgebra_sparse::CsrMatrix;

/// Create a tridiagonal CSR matrix (3 non-zeros per row)
/// Common in 1D heat/diffusion equations
fn create_tridiagonal_csr_f64(size: usize) -> CsrMatrix<f64> {
    let mut builder = SparseMatrixBuilder::new(size, size);

    for i in 0..size {
        builder.add_entry(i, i, 2.0).unwrap();
        if i > 0 {
            builder.add_entry(i, i - 1, -1.0).unwrap();
        }
        if i < size - 1 {
            builder.add_entry(i, i + 1, -1.0).unwrap();
        }
    }

    builder.build().unwrap()
}

/// Create a tridiagonal CSR matrix (f32 for SIMD)
#[allow(dead_code)]
fn create_tridiagonal_csr_f32(size: usize) -> CsrMatrix<f32> {
    let mut builder = SparseMatrixBuilder::new(size, size);

    for i in 0..size {
        builder.add_entry(i, i, 2.0f32).unwrap();
        if i > 0 {
            builder.add_entry(i, i - 1, -1.0f32).unwrap();
        }
        if i < size - 1 {
            builder.add_entry(i, i + 1, -1.0f32).unwrap();
        }
    }

    builder.build().unwrap()
}

/// Create a pentadiagonal sparse matrix (5 non-zeros per row)
/// Common in 2D Laplacian operators
fn create_pentadiagonal_csr_f64(n: usize) -> CsrMatrix<f64> {
    let size = n * n;
    let mut builder = SparseMatrixBuilder::new(size, size);

    for i in 0..size {
        let row = i / n;
        let col = i % n;

        if col > 0 {
            builder.add_entry(i, i - 1, -1.0).unwrap();
        }
        if row > 0 {
            builder.add_entry(i, i - n, -1.0).unwrap();
        }
        builder.add_entry(i, i, 4.0).unwrap();
        if row < n - 1 {
            builder.add_entry(i, i + n, -1.0).unwrap();
        }
        if col < n - 1 {
            builder.add_entry(i, i + 1, -1.0).unwrap();
        }
    }

    builder.build().unwrap()
}

/// Create a pentadiagonal sparse matrix (f32 for SIMD)
fn create_pentadiagonal_csr_f32(n: usize) -> CsrMatrix<f32> {
    let size = n * n;
    let mut builder = SparseMatrixBuilder::new(size, size);

    for i in 0..size {
        let row = i / n;
        let col = i % n;

        if col > 0 {
            builder.add_entry(i, i - 1, -1.0f32).unwrap();
        }
        if row > 0 {
            builder.add_entry(i, i - n, -1.0f32).unwrap();
        }
        builder.add_entry(i, i, 4.0f32).unwrap();
        if row < n - 1 {
            builder.add_entry(i, i + n, -1.0f32).unwrap();
        }
        if col < n - 1 {
            builder.add_entry(i, i + 1, -1.0f32).unwrap();
        }
    }

    builder.build().unwrap()
}

/// Benchmark scalar SpMV baseline (f64)
fn bench_scalar_spmv(c: &mut Criterion) {
    let mut group = c.benchmark_group("spmv_scalar_f64");

    // Small: 100x100 (300 non-zeros)
    let matrix = create_tridiagonal_csr_f64(100);
    let x = DVector::from_element(100, 1.0);
    let mut y = DVector::zeros(100);
    group.throughput(Throughput::Elements(matrix.nnz() as u64));
    group.bench_function("tridiagonal_100", |b| {
        b.iter(|| {
            spmv(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });

    // Medium: 500x500 (1500 non-zeros)
    let matrix = create_tridiagonal_csr_f64(500);
    let x = DVector::from_element(500, 1.0);
    let mut y = DVector::zeros(500);
    group.throughput(Throughput::Elements(matrix.nnz() as u64));
    group.bench_function("tridiagonal_500", |b| {
        b.iter(|| {
            spmv(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });

    // Large: 2000x2000 (6000 non-zeros)
    let matrix = create_tridiagonal_csr_f64(2000);
    let x = DVector::from_element(2000, 1.0);
    let mut y = DVector::zeros(2000);
    group.throughput(Throughput::Elements(matrix.nnz() as u64));
    group.bench_function("tridiagonal_2000", |b| {
        b.iter(|| {
            spmv(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });

    group.finish();
}

/// Benchmark pentadiagonal matrices (denser, more realistic CFD)
fn bench_pentadiagonal(c: &mut Criterion) {
    let mut group = c.benchmark_group("spmv_pentadiagonal");

    // 32x32 grid = 1024 unknowns
    let matrix = create_pentadiagonal_csr_f64(32);
    let x = DVector::from_element(1024, 1.0);
    let mut y = DVector::zeros(1024);
    group.throughput(Throughput::Elements(matrix.nnz() as u64));
    group.bench_function("scalar_32x32", |b| {
        b.iter(|| {
            spmv(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        let matrix_f32 = create_pentadiagonal_csr_f32(32);
        let x_f32 = DVector::from_element(1024, 1.0f32);
        let mut y_f32 = DVector::zeros(1024);
        group.bench_function("simd_32x32_deprecated", |b| {
            b.iter(|| {
                // SIMD implementation was removed due to 27-32% regression
                // Replaced with scalar f32 for comparison
                spmv(black_box(&matrix_f32), black_box(&x_f32), black_box(&mut y_f32));
            });
        });
    }

    // 64x64 grid = 4096 unknowns
    let matrix = create_pentadiagonal_csr_f64(64);
    let x = DVector::from_element(4096, 1.0);
    let mut y = DVector::zeros(4096);
    group.throughput(Throughput::Elements(matrix.nnz() as u64));
    group.bench_function("scalar_64x64", |b| {
        b.iter(|| {
            spmv(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        let matrix_f32 = create_pentadiagonal_csr_f32(64);
        let x_f32 = DVector::from_element(4096, 1.0f32);
        let mut y_f32 = DVector::zeros(4096);
        group.bench_function("simd_64x64_deprecated", |b| {
            b.iter(|| {
                // SIMD implementation was removed due to 27-32% regression
                // Replaced with scalar f32 for comparison
                spmv(black_box(&matrix_f32), black_box(&x_f32), black_box(&mut y_f32));
            });
        });
    }

    group.finish();
}

/// Benchmark parallel SpMV (rayon-based, Sprint 1.67.0)
fn bench_parallel_spmv(c: &mut Criterion) {
    let mut group = c.benchmark_group("spmv_parallel");

    // Medium: 1000x1000 tridiagonal (3000 non-zeros)
    // Parallel benefit starts to show
    let matrix = create_tridiagonal_csr_f64(1000);
    let x = DVector::from_element(1000, 1.0);
    let mut y = DVector::zeros(1000);
    group.throughput(Throughput::Elements(matrix.nnz() as u64));
    group.bench_function("tridiagonal_1000_scalar", |b| {
        b.iter(|| {
            spmv(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });
    group.bench_function("tridiagonal_1000_parallel", |b| {
        b.iter(|| {
            spmv_parallel(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });

    // Large: 5000x5000 tridiagonal (15000 non-zeros)
    // Expected 3-8x speedup on 4-8 cores
    let matrix = create_tridiagonal_csr_f64(5000);
    let x = DVector::from_element(5000, 1.0);
    let mut y = DVector::zeros(5000);
    group.throughput(Throughput::Elements(matrix.nnz() as u64));
    group.bench_function("tridiagonal_5000_scalar", |b| {
        b.iter(|| {
            spmv(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });
    group.bench_function("tridiagonal_5000_parallel", |b| {
        b.iter(|| {
            spmv_parallel(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });

    // XL: 10000x10000 tridiagonal (30000 non-zeros)
    // Target: demonstrate scaling benefit
    let matrix = create_tridiagonal_csr_f64(10000);
    let x = DVector::from_element(10000, 1.0);
    let mut y = DVector::zeros(10000);
    group.throughput(Throughput::Elements(matrix.nnz() as u64));
    group.bench_function("tridiagonal_10000_scalar", |b| {
        b.iter(|| {
            spmv(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });
    group.bench_function("tridiagonal_10000_parallel", |b| {
        b.iter(|| {
            spmv_parallel(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });

    group.finish();
}

/// Benchmark parallel pentadiagonal (realistic CFD matrices)
fn bench_parallel_pentadiagonal(c: &mut Criterion) {
    let mut group = c.benchmark_group("spmv_parallel_pentadiagonal");

    // 50x50 grid = 2500 unknowns (12500 non-zeros)
    let matrix = create_pentadiagonal_csr_f64(50);
    let x = DVector::from_element(2500, 1.0);
    let mut y = DVector::zeros(2500);
    group.throughput(Throughput::Elements(matrix.nnz() as u64));
    group.bench_function("50x50_scalar", |b| {
        b.iter(|| {
            spmv(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });
    group.bench_function("50x50_parallel", |b| {
        b.iter(|| {
            spmv_parallel(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });

    // 100x100 grid = 10000 unknowns (50000 non-zeros)
    // Target workload for CFD simulations
    let matrix = create_pentadiagonal_csr_f64(100);
    let x = DVector::from_element(10000, 1.0);
    let mut y = DVector::zeros(10000);
    group.throughput(Throughput::Elements(matrix.nnz() as u64));
    group.bench_function("100x100_scalar", |b| {
        b.iter(|| {
            spmv(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });
    group.bench_function("100x100_parallel", |b| {
        b.iter(|| {
            spmv_parallel(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });

    // 200x200 grid = 40000 unknowns (200000 non-zeros)
    // Large-scale CFD simulation
    let matrix = create_pentadiagonal_csr_f64(200);
    let x = DVector::from_element(40000, 1.0);
    let mut y = DVector::zeros(40000);
    group.throughput(Throughput::Elements(matrix.nnz() as u64));
    group.bench_function("200x200_scalar", |b| {
        b.iter(|| {
            spmv(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });
    group.bench_function("200x200_parallel", |b| {
        b.iter(|| {
            spmv_parallel(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });

    group.finish();
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
criterion_group!(
    benches,
    bench_scalar_spmv,
    bench_pentadiagonal,
    bench_parallel_spmv,
    bench_parallel_pentadiagonal
);

#[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
criterion_group!(
    benches,
    bench_scalar_spmv,
    bench_pentadiagonal,
    bench_parallel_spmv,
    bench_parallel_pentadiagonal
);

criterion_main!(benches);
