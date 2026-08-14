//! Criterion comparison: CSR-shaped SIMPLE coupling store vs the historical
//! jagged `Vec<Vec<(usize, f64)>>` layout (ATLAS-ARCH-008, CFDrs slice).
//!
//! `SimplePreconditioner` traverses the divergence and gradient coupling
//! blocks once per application — every Krylov iteration. The flat store
//! removes `n_pressure` per-row allocations at construction and replaces
//! pointer chasing with one contiguous pass. The historical jagged layout is
//! kept inline as the baseline, and a parity gate pins both layouts to the
//! same outputs before any measurement is trusted.

#![expect(
    missing_docs,
    reason = "ratchet CFDRS-DOCS-1: criterion macro-generated items"
)]
// The jagged baseline below intentionally mirrors the historical
// implementation loop-for-loop so the comparison measures layout, not
// algorithm; the indexed loops and slice copies are part of that faithful
// transcription.
#![expect(
    clippy::needless_range_loop,
    reason = "baseline mirrors the historical jagged implementation verbatim"
)]
#![expect(
    clippy::manual_memcpy,
    reason = "baseline mirrors the historical jagged implementation verbatim"
)]
#![expect(
    clippy::unwrap_used,
    reason = "ratchet CFDRS-UNWRAP-1: pre-existing debt"
)]

use cfd_math::linear_solver::SimplePreconditioner;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use leto::Array1;
use leto_ops::CsrMatrix;

/// Historical jagged SIMPLE layout, kept inline as the comparison baseline.
struct JaggedSimplePreconditioner {
    momentum_inv: Vec<f64>,
    schur_diag_inv: Vec<f64>,
    divergence_rows: Vec<Vec<(usize, f64)>>,
    gradient_columns: Vec<Vec<(usize, f64)>>,
    n_velocity: usize,
    n_pressure: usize,
}

fn diagonal(matrix: &CsrMatrix<f64>, row: usize) -> f64 {
    let row_range = matrix.row(row);
    for (column, &value) in row_range.col_indices().iter().zip(row_range.values()) {
        if *column == row {
            return value;
        }
    }
    0.0
}

impl JaggedSimplePreconditioner {
    fn new(matrix: &CsrMatrix<f64>, n_velocity: usize, n_pressure: usize) -> Self {
        let eps = 1e-14;
        let mut momentum_inv = vec![0.0; n_velocity];
        for index in 0..n_velocity {
            let d = diagonal(matrix, index);
            momentum_inv[index] = if d.abs() > eps { 1.0 / d } else { 1.0 };
        }

        let mut divergence_rows = Vec::with_capacity(n_pressure);
        for i in 0..n_pressure {
            let global_row = n_velocity + i;
            let row = matrix.row(global_row);
            let entries: Vec<(usize, f64)> = row
                .col_indices()
                .iter()
                .zip(row.values())
                .filter(|(&c, _)| c < n_velocity)
                .map(|(&c, &v)| (c, v))
                .collect();
            divergence_rows.push(entries);
        }
        let mut gradient_columns = vec![Vec::new(); n_pressure];
        for velocity_row in 0..n_velocity {
            let row = matrix.row(velocity_row);
            for (&column, &value) in row.col_indices().iter().zip(row.values()) {
                if let Some(pressure_column) = column
                    .checked_sub(n_velocity)
                    .filter(|&index| index < n_pressure)
                {
                    gradient_columns[pressure_column].push((velocity_row, value));
                }
            }
        }

        let mut schur_diag_inv = vec![0.0; n_pressure];
        for i in 0..n_pressure {
            let mut s_ii = diagonal(matrix, n_velocity + i);
            for &(velocity_index, divergence_value) in &divergence_rows[i] {
                if let Some(&(_, gradient_value)) = gradient_columns[i]
                    .iter()
                    .find(|&&(index, _)| index == velocity_index)
                {
                    s_ii -= divergence_value * momentum_inv[velocity_index] * gradient_value;
                }
            }
            schur_diag_inv[i] = if s_ii.abs() > eps { 1.0 / s_ii } else { 1.0 };
        }

        Self {
            momentum_inv,
            schur_diag_inv,
            divergence_rows,
            gradient_columns,
            n_velocity,
            n_pressure,
        }
    }

    fn apply(&self, b: &Array1<f64>) -> Array1<f64> {
        // Faithful transcription of the pre-change library apply: the same
        // jagged coupling store and the same `leto::Array1` working vectors,
        // so the comparison isolates the coupling-store layout.
        let mut u_star = Array1::zeros([self.n_velocity]);
        for index in 0..self.n_velocity {
            u_star[index] = b[index] * self.momentum_inv[index];
        }

        let mut rhs_p = Array1::zeros([self.n_pressure]);
        for index in 0..self.n_pressure {
            rhs_p[index] = b[self.n_velocity + index];
        }
        for i in 0..self.n_pressure {
            let mut b_u = 0.0;
            for &(velocity_index, divergence_value) in &self.divergence_rows[i] {
                b_u += divergence_value * u_star[velocity_index];
            }
            rhs_p[i] -= b_u;
        }
        let mut p = Array1::zeros([self.n_pressure]);
        for index in 0..self.n_pressure {
            p[index] = rhs_p[index] * self.schur_diag_inv[index];
        }

        let mut u_corrected = u_star;
        for i in 0..self.n_pressure {
            for &(velocity_index, gradient_value) in &self.gradient_columns[i] {
                u_corrected[velocity_index] -=
                    self.momentum_inv[velocity_index] * gradient_value * p[i];
            }
        }

        let mut x = Array1::zeros([self.n_velocity + self.n_pressure]);
        for index in 0..self.n_velocity {
            x[index] = u_corrected[index];
        }
        for index in 0..self.n_pressure {
            x[self.n_velocity + index] = p[index];
        }
        x
    }
}

/// CSR stores transcribed from the matrix with plain slices, mirroring the
/// library's two-pass build.
struct FlatStores {
    divergence_offsets: Vec<usize>,
    divergence_entries: Vec<(usize, f64)>,
    gradient_offsets: Vec<usize>,
    gradient_entries: Vec<(usize, f64)>,
    momentum_inv: Vec<f64>,
    schur_diag_inv: Vec<f64>,
}

fn build_flat_stores(matrix: &CsrMatrix<f64>, n_velocity: usize, n_pressure: usize) -> FlatStores {
    let eps = 1e-14;
    let mut momentum_inv = vec![0.0; n_velocity];
    for index in 0..n_velocity {
        let d = diagonal(matrix, index);
        momentum_inv[index] = if d.abs() > eps { 1.0 / d } else { 1.0 };
    }

    // Divergence block: rows of the lower block filtered to velocity columns.
    let mut lengths = Vec::with_capacity(n_pressure);
    for i in 0..n_pressure {
        let count = matrix
            .row(n_velocity + i)
            .col_indices()
            .iter()
            .filter(|&&column| column < n_velocity)
            .count();
        lengths.push(count);
    }
    let mut divergence_offsets = Vec::with_capacity(n_pressure + 1);
    let mut running = 0usize;
    divergence_offsets.push(0);
    for &length in &lengths {
        running += length;
        divergence_offsets.push(running);
    }
    let mut divergence_entries = Vec::with_capacity(running);
    for i in 0..n_pressure {
        let row = matrix.row(n_velocity + i);
        for (&column, &value) in row.col_indices().iter().zip(row.values()) {
            if column < n_velocity {
                divergence_entries.push((column, value));
            }
        }
    }

    // Gradient block: per-pressure-column transposed scan.
    let mut column_lengths = vec![0usize; n_pressure];
    for velocity_row in 0..n_velocity {
        let row = matrix.row(velocity_row);
        for &column in row.col_indices() {
            if let Some(pressure_column) = column
                .checked_sub(n_velocity)
                .filter(|&index| index < n_pressure)
            {
                column_lengths[pressure_column] += 1;
            }
        }
    }
    let mut gradient_offsets = Vec::with_capacity(n_pressure + 1);
    let mut running = 0usize;
    gradient_offsets.push(0);
    for &length in &column_lengths {
        running += length;
        gradient_offsets.push(running);
    }
    let mut gradient_entries = vec![(0usize, 0.0); running];
    let mut cursors = gradient_offsets.clone();
    for velocity_row in 0..n_velocity {
        let row = matrix.row(velocity_row);
        for (&column, &value) in row.col_indices().iter().zip(row.values()) {
            if let Some(pressure_column) = column
                .checked_sub(n_velocity)
                .filter(|&index| index < n_pressure)
            {
                let position = cursors[pressure_column];
                gradient_entries[position] = (velocity_row, value);
                cursors[pressure_column] += 1;
            }
        }
    }

    // Diagonal Schur complement, transcribed from the library.
    let mut schur_diag_inv = vec![0.0; n_pressure];
    for i in 0..n_pressure {
        let mut s_ii = diagonal(matrix, n_velocity + i);
        for (velocity_index, divergence_value) in
            &divergence_entries[divergence_offsets[i]..divergence_offsets[i + 1]]
        {
            if let Some((_, gradient_value)) = gradient_entries
                [gradient_offsets[i]..gradient_offsets[i + 1]]
                .iter()
                .find(|&&(index, _)| index == *velocity_index)
            {
                s_ii -= divergence_value * momentum_inv[*velocity_index] * gradient_value;
            }
        }
        schur_diag_inv[i] = if s_ii.abs() > eps { 1.0 / s_ii } else { 1.0 };
    }

    FlatStores {
        divergence_offsets,
        divergence_entries,
        gradient_offsets,
        gradient_entries,
        momentum_inv,
        schur_diag_inv,
    }
}

/// CSR apply transcribed with plain slices and explicit offset-range loops.
///
/// This isolates the layout from the `row()` iterator helper and the
/// `Array1`/`DiagonalPreconditioner` wrappers so the comparison attributes
/// any difference to the store shape itself. Entries stay AoS-packed
/// `(index, value)` pairs: the split SoA form loses ~14-35% on the
/// sequential store read, while the flat AoS form holds parity with the
/// jagged baseline.
fn csr_flat_apply(
    divergence_offsets: &[usize],
    divergence_entries: &[(usize, f64)],
    gradient_offsets: &[usize],
    gradient_entries: &[(usize, f64)],
    momentum_inv: &[f64],
    schur_diag_inv: &[f64],
    b: &[f64],
) -> Vec<f64> {
    let n_velocity = momentum_inv.len();
    let n_pressure = schur_diag_inv.len();
    let mut u_star = vec![0.0; n_velocity];
    for index in 0..n_velocity {
        u_star[index] = b[index] * momentum_inv[index];
    }

    let mut rhs_p = vec![0.0; n_pressure];
    for index in 0..n_pressure {
        rhs_p[index] = b[n_velocity + index];
    }
    for i in 0..n_pressure {
        let mut b_u = 0.0;
        for (velocity_index, divergence_value) in
            &divergence_entries[divergence_offsets[i]..divergence_offsets[i + 1]]
        {
            b_u += divergence_value * u_star[*velocity_index];
        }
        rhs_p[i] -= b_u;
    }
    let mut p = vec![0.0; n_pressure];
    for index in 0..n_pressure {
        p[index] = rhs_p[index] * schur_diag_inv[index];
    }

    let mut u_corrected = u_star;
    for i in 0..n_pressure {
        for (velocity_index, gradient_value) in
            &gradient_entries[gradient_offsets[i]..gradient_offsets[i + 1]]
        {
            u_corrected[*velocity_index] -= momentum_inv[*velocity_index] * gradient_value * p[i];
        }
    }

    let mut x = vec![0.0; n_velocity + n_pressure];
    x[..n_velocity].copy_from_slice(&u_corrected);
    x[n_velocity..].copy_from_slice(&p);
    x
}

/// 2-D Taylor–Hood-like saddle system: two velocity DOFs per cell plus one
/// pressure DOF, nearest-neighbour momentum coupling and cell-wise
/// divergence/gradient entries.
fn taylor_hood_saddle(nx: usize, ny: usize) -> (CsrMatrix<f64>, usize, usize) {
    let n_cells = nx * ny;
    let n_velocity = 2 * n_cells;
    let n_pressure = n_cells;
    let n = n_velocity + n_pressure;

    let mut row_offsets = Vec::with_capacity(n + 1);
    let mut col_indices = Vec::new();
    let mut values = Vec::new();
    row_offsets.push(0);

    for row in 0..n {
        let mut entries: Vec<(usize, f64)> = Vec::with_capacity(8);
        if row < n_velocity {
            let cell = row / 2;
            let (cx, cy) = (cell % nx, cell / nx);
            entries.push((row, 4.0));
            if cx > 0 {
                entries.push((2 * (cell - 1), -1.0));
            }
            if cx + 1 < nx {
                entries.push((2 * (cell + 1), -1.0));
            }
            if cy > 0 {
                entries.push((2 * (cell - nx), -1.0));
            }
            if cy + 1 < ny {
                entries.push((2 * (cell + nx), -1.0));
            }
            entries.push((n_velocity + cell, 1.0));
        } else {
            let cell = row - n_velocity;
            entries.push((2 * cell, 1.0));
            entries.push((2 * cell + 1, 1.0));
            entries.push((row, 1.0));
        }
        // leto's CSR contract requires strictly increasing column indices per
        // row; the momentum rows above are built in neighbour-visit order.
        entries.sort_unstable_by_key(|&(column, _)| column);
        for (column, value) in entries {
            col_indices.push(column);
            values.push(value);
        }
        row_offsets.push(col_indices.len());
    }

    let matrix = CsrMatrix::from_parts(values, col_indices, row_offsets, n, n).unwrap();
    (matrix, n_velocity, n_pressure)
}

fn bench_simple_preconditioner(c: &mut Criterion) {
    let fixtures = [16usize, 32, 64]
        .iter()
        .map(|&nx| {
            let (matrix, n_velocity, n_pressure) = taylor_hood_saddle(nx, nx);
            let n = n_velocity + n_pressure;
            let b = vec![1.0f64; n];
            let b_array = Array1::from_shape_vec([n], b.clone()).unwrap();
            let csr = SimplePreconditioner::new(&matrix, n_velocity, n_pressure).unwrap();
            let jagged = JaggedSimplePreconditioner::new(&matrix, n_velocity, n_pressure);

            let flat = build_flat_stores(&matrix, n_velocity, n_pressure);

            // Parity gate: the CSR store and the raw flat transcription must
            // reproduce the jagged outputs before any measurement is trusted.
            let csr_output = csr.apply(&b_array).unwrap();
            let jagged_output = jagged.apply(&b_array);
            let flat_output = csr_flat_apply(
                &flat.divergence_offsets,
                &flat.divergence_entries,
                &flat.gradient_offsets,
                &flat.gradient_entries,
                &flat.momentum_inv,
                &flat.schur_diag_inv,
                &b,
            );
            for index in 0..n {
                let lhs = csr_output[index];
                let mid = jagged_output[index];
                let rhs = flat_output[index];
                assert!(
                    (lhs - rhs).abs() < 1e-12 && (mid - rhs).abs() < 1e-12,
                    "SIMPLE parity drift between CSR, flat, and jagged layouts"
                );
            }

            (
                nx, matrix, n_velocity, n_pressure, b, b_array, csr, jagged, flat,
            )
        })
        .collect::<Vec<_>>();

    let mut build = c.benchmark_group("simple_coupling_build");
    for (nx, matrix, n_velocity, n_pressure, ..) in &fixtures {
        let label = format!("{nx}x{nx}");
        build.bench_with_input(BenchmarkId::new("csr", &label), &(), |bencher, _| {
            bencher.iter(|| {
                black_box(
                    SimplePreconditioner::new(black_box(matrix), *n_velocity, *n_pressure).unwrap(),
                );
            });
        });
        build.bench_with_input(BenchmarkId::new("jagged", &label), &(), |bencher, _| {
            bencher.iter(|| {
                black_box(JaggedSimplePreconditioner::new(
                    black_box(matrix),
                    *n_velocity,
                    *n_pressure,
                ));
            });
        });
    }
    build.finish();

    let mut apply = c.benchmark_group("simple_coupling_apply");
    for (nx, _, _, _, b, b_array, csr, jagged, flat) in &fixtures {
        let label = format!("{nx}x{nx}");
        apply.bench_with_input(BenchmarkId::new("csr", &label), &(), |bencher, _| {
            bencher.iter(|| {
                black_box(csr.apply(black_box(b_array)).unwrap());
            });
        });
        apply.bench_with_input(BenchmarkId::new("csr_flat", &label), &(), |bencher, _| {
            bencher.iter(|| {
                black_box(csr_flat_apply(
                    black_box(&flat.divergence_offsets),
                    black_box(&flat.divergence_entries),
                    black_box(&flat.gradient_offsets),
                    black_box(&flat.gradient_entries),
                    black_box(&flat.momentum_inv),
                    black_box(&flat.schur_diag_inv),
                    black_box(b),
                ));
            });
        });
        apply.bench_with_input(BenchmarkId::new("jagged", &label), &(), |bencher, _| {
            bencher.iter(|| {
                black_box(jagged.apply(black_box(b_array)));
            });
        });
    }
    apply.finish();
}

criterion_group!(benches, bench_simple_preconditioner);
criterion_main!(benches);
