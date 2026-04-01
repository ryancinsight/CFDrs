//! Criterion benchmarks for 3D Delaunay tetrahedralization.
//!
//! Measures the Bowyer-Watson insertion throughput at increasing point counts
//! to validate the expected O(n log n) scaling of the `cfd-mesh` implementation.
//!
//! # Analytical Reference
//!
//! The Bowyer-Watson algorithm has expected O(n log n) time complexity for
//! uniformly distributed point sets in ℝ³ (de Berg et al., *Computational
//! Geometry*, 3rd ed., Springer, 2008, Chapter 9).

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use nalgebra::Point3;

/// Generate a deterministic set of BCC lattice points in the unit cube.
///
/// Uses a body-centred cubic arrangement for uniform spatial distribution,
/// matching the lattice generator used in `cfd-mesh::application::delaunay::dim3::lattice`.
fn generate_bcc_points(n_side: usize) -> Vec<Point3<f64>> {
    let h = 1.0 / (n_side as f64);
    let mut points = Vec::with_capacity(n_side * n_side * n_side * 2);

    for iz in 0..n_side {
        for iy in 0..n_side {
            for ix in 0..n_side {
                let x = (ix as f64 + 0.5) * h;
                let y = (iy as f64 + 0.5) * h;
                let z = (iz as f64 + 0.5) * h;
                points.push(Point3::new(x, y, z));

                // BCC offset point
                let xo = x + 0.5 * h;
                let yo = y + 0.5 * h;
                let zo = z + 0.5 * h;
                if xo < 1.0 && yo < 1.0 && zo < 1.0 {
                    points.push(Point3::new(xo, yo, zo));
                }
            }
        }
    }
    points
}

fn bench_bowyer_watson_3d(c: &mut Criterion) {
    let mut group = c.benchmark_group("delaunay_3d_bowyer_watson");

    // Increasing point counts via BCC lattice side length
    for &n_side in &[4, 6, 8, 12] {
        let points = generate_bcc_points(n_side);
        let n_points = points.len();

        group.bench_with_input(
            BenchmarkId::new("insert", n_points),
            &points,
            |b, pts| {
                b.iter(|| {
                    let min = Point3::new(-0.5, -0.5, -0.5);
                    let max = Point3::new(1.5, 1.5, 1.5);
                    let mut bw =
                        cfd_mesh::application::delaunay::dim3::tetrahedralize::BowyerWatson3D::<
                            f64,
                        >::new(min, max);

                    for p in pts {
                        bw.insert_point(*p);
                    }

                    let (vertices, tets) = bw.finalize();
                    black_box((vertices.len(), tets.len()));
                });
            },
        );
    }

    group.finish();
}

fn bench_bcc_lattice_generation(c: &mut Criterion) {
    let mut group = c.benchmark_group("delaunay_3d_lattice_gen");

    for &n_side in &[4, 8, 16, 32] {
        let n_points = n_side * n_side * n_side * 2;
        group.bench_with_input(
            BenchmarkId::new("bcc_generate", n_points),
            &n_side,
            |b, &n| {
                b.iter(|| {
                    let pts = generate_bcc_points(n);
                    black_box(pts.len());
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_bowyer_watson_3d, bench_bcc_lattice_generation);
criterion_main!(benches);
