#![allow(missing_docs, clippy::doc_markdown, clippy::ignored_unit_patterns)]
//! Criterion: AMG smoother dispatch on the preconditioner hot path.
//!
//! `AlgebraicMultigrid` holds one smoother per level and calls it on every
//! pre- and post-smoothing sweep of every V-cycle, so the smoother's
//! dispatch shape is on the hottest path a Krylov solve touches. This bench
//! measures the two ends of that path for each configured smoother type:
//! hierarchy construction (`AlgebraicMultigrid::new`, which builds the
//! per-level smoothers) and one preconditioner application (`apply_to`).
//!
//! Workload: the 2-D five-point Poisson stencil on a 64x64 grid — the
//! canonical SPD multigrid problem, sized to keep the bench inside the
//! committed suite budget while exercising a real multi-level hierarchy.

use cfd_math::linear_solver::preconditioners::multigrid::SmootherType;
use cfd_math::linear_solver::{AMGConfig, AlgebraicMultigrid};
use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use leto::Array1;
use leto_ops::CsrMatrix;

/// Five-point 2-D Poisson matrix on an `n` by `n` interior grid.
fn poisson_matrix(n: usize) -> CsrMatrix<f64> {
    let size = n * n;
    let mut values = Vec::with_capacity(size * 5);
    let mut col_indices = Vec::with_capacity(size * 5);
    let mut row_offsets = Vec::with_capacity(size + 1);
    row_offsets.push(0usize);
    let mut entries = 0usize;
    for row in 0..size {
        let (i, j) = (row / n, row % n);
        // Entries in strictly increasing column order, as the CSR builder
        // requires: south (row-n) < west (row-1) < diagonal < east (row+1)
        // < north (row+n).
        if i > 0 {
            values.push(-1.0);
            col_indices.push(row - n);
            entries += 1;
        }
        if j > 0 {
            values.push(-1.0);
            col_indices.push(row - 1);
            entries += 1;
        }
        values.push(4.0);
        col_indices.push(row);
        entries += 1;
        if j + 1 < n {
            values.push(-1.0);
            col_indices.push(row + 1);
            entries += 1;
        }
        if i + 1 < n {
            values.push(-1.0);
            col_indices.push(row + n);
            entries += 1;
        }
        row_offsets.push(entries);
    }
    CsrMatrix::from_parts(values, col_indices, row_offsets, size, size).expect("expected value")
}

fn config_for(smoother: SmootherType) -> AMGConfig {
    AMGConfig {
        smoother_type: smoother,
        ..AMGConfig::default()
    }
}

fn smoother_cases() -> [(String, SmootherType); 5] {
    [
        ("gauss_seidel".to_owned(), SmootherType::GaussSeidel),
        (
            "symmetric_gauss_seidel".to_owned(),
            SmootherType::SymmetricGaussSeidel,
        ),
        ("jacobi".to_owned(), SmootherType::Jacobi),
        ("sor".to_owned(), SmootherType::SOR),
        ("chebyshev".to_owned(), SmootherType::Chebyshev),
    ]
}

fn bench_amg_smoother(criterion: &mut Criterion) {
    let matrix = poisson_matrix(64);
    let residual = Array1::<f64>::from_shape_fn([matrix.nrows()], |[i]| {
        // Smooth right-hand side: a low-frequency mode the smoother damps.
        (std::f64::consts::TAU * (i as f64) / matrix.nrows() as f64).cos()
    });

    let mut construction = criterion.benchmark_group("amg_hierarchy_construction");
    for (name, smoother) in smoother_cases() {
        construction.bench_with_input(
            BenchmarkId::new("poisson_64", name),
            &smoother,
            |bench, smoother| {
                bench.iter(|| {
                    let amg =
                        AlgebraicMultigrid::<f64>::new(black_box(&matrix), config_for(*smoother))
                            .expect("expected value");
                    black_box(amg);
                });
            },
        );
    }
    construction.finish();

    let mut application = criterion.benchmark_group("amg_v_cycle_apply");
    for (name, smoother) in smoother_cases() {
        let amg =
            AlgebraicMultigrid::<f64>::new(&matrix, config_for(smoother)).expect("expected value");
        application.bench_with_input(BenchmarkId::new("poisson_64", name), &amg, |bench, amg| {
            let mut correction = Array1::<f64>::zeros([matrix.nrows()]);
            bench.iter(|| {
                amg.apply_to(black_box(&residual), black_box(&mut correction))
                    .expect("expected value");
            });
        });
    }
    application.finish();
}

criterion_group!(benches, bench_amg_smoother);
criterion_main!(benches);
