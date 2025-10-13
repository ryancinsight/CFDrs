//! SIMD Sparse Matrix-Vector Multiply (SpMV) Benchmarks
//!
//! These benchmarks validate the SIMD optimization implemented in Sprint 1.41.0.
//! Expected performance gains:
//! - AVX2 (256-bit): 2-4x speedup vs scalar
//! - SSE4.1 (128-bit): 1.5-2x speedup vs scalar
//!
//! Test matrices:
//! - Tridiagonal (typical 1D diffusion)
//! - Pentadiagonal (typical 2D diffusion)

use cfd_math::sparse::{spmv, spmv_f32_simd, SparseMatrixBuilder};
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

/// Benchmark SIMD SpMV (AVX2/SSE4.1 with runtime dispatch)
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
fn bench_simd_spmv(c: &mut Criterion) {
    let mut group = c.benchmark_group("spmv_simd_f32");
    
    // Small: 100x100 tridiagonal
    let matrix = create_tridiagonal_csr_f32(100);
    let x = DVector::from_element(100, 1.0f32);
    let mut y = DVector::zeros(100);
    group.throughput(Throughput::Elements(matrix.nnz() as u64));
    group.bench_function("tridiagonal_100", |b| {
        b.iter(|| {
            spmv_f32_simd(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });
    
    // Medium: 500x500 tridiagonal
    let matrix = create_tridiagonal_csr_f32(500);
    let x = DVector::from_element(500, 1.0f32);
    let mut y = DVector::zeros(500);
    group.throughput(Throughput::Elements(matrix.nnz() as u64));
    group.bench_function("tridiagonal_500", |b| {
        b.iter(|| {
            spmv_f32_simd(black_box(&matrix), black_box(&x), black_box(&mut y));
        });
    });
    
    // Large: 2000x2000 tridiagonal
    let matrix = create_tridiagonal_csr_f32(2000);
    let x = DVector::from_element(2000, 1.0f32);
    let mut y = DVector::zeros(2000);
    group.throughput(Throughput::Elements(matrix.nnz() as u64));
    group.bench_function("tridiagonal_2000", |b| {
        b.iter(|| {
            spmv_f32_simd(black_box(&matrix), black_box(&x), black_box(&mut y));
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
        group.bench_function("simd_32x32", |b| {
            b.iter(|| {
                spmv_f32_simd(black_box(&matrix_f32), black_box(&x_f32), black_box(&mut y_f32));
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
        group.bench_function("simd_64x64", |b| {
            b.iter(|| {
                spmv_f32_simd(black_box(&matrix_f32), black_box(&x_f32), black_box(&mut y_f32));
            });
        });
    }
    
    group.finish();
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
criterion_group!(
    benches,
    bench_scalar_spmv,
    bench_simd_spmv,
    bench_pentadiagonal
);

#[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
criterion_group!(
    benches,
    bench_scalar_spmv,
    bench_pentadiagonal
);

criterion_main!(benches);
