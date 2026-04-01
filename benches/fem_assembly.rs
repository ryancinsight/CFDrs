//! Criterion benchmarks for FEM stiffness matrix assembly.
//!
//! Exercises real `cfd-3d::fem` code paths: mesh construction via
//! `StructuredGridBuilder`, vertex-index extraction, and full Stokes
//! FEM solve via `FemSolver::solve`.
//!
//! # Analytical Reference
//!
//! For a structured hexahedral mesh of N³ cells decomposed into 6N³
//! tetrahedra, P1 FEM assembly scales as O(N³) in element count.
//! The FEM solve (GMRES) scales as O(N³ · k) where k is the Krylov
//! iteration count (typically O(N) for ill-conditioned saddle-point
//! systems without AMG preconditioning).

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use std::collections::HashMap;

/// Build a unit-cube structured hex mesh of given resolution and return
/// the vertex count for P1 elements.
fn build_structured_cube(
    n: usize,
) -> (cfd_mesh::IndexedMesh<f64>, usize) {
    let mesh = cfd_mesh::domain::grid::StructuredGridBuilder::new(n, n, n)
        .build()
        .expect("StructuredGridBuilder::build must succeed for n >= 2");
    let n_corner = mesh.vertex_count();
    (mesh, n_corner)
}

fn bench_extract_vertex_indices(c: &mut Criterion) {
    let mut group = c.benchmark_group("fem_extract_vertex_indices");

    for &n in &[3, 5, 8] {
        let (mesh, n_corner) = build_structured_cube(n);

        group.bench_with_input(
            BenchmarkId::new("extract", mesh.cell_count()),
            &(mesh, n_corner),
            |b, (mesh, n_corner)| {
                b.iter(|| {
                    let mut total_indices = 0usize;
                    for cell in &mesh.cells {
                        let indices =
                            cfd_3d::fem::mesh_utils::extract_vertex_indices(cell, mesh, *n_corner)
                                .expect("extract_vertex_indices must succeed");
                        total_indices += indices.len();
                    }
                    black_box(total_indices);
                });
            },
        );
    }

    group.finish();
}

fn bench_fem_solve_structured(c: &mut Criterion) {
    let mut group = c.benchmark_group("fem_stokes_solve");
    group.sample_size(10); // FEM solves are expensive

    for &n in &[3, 4] {
        let (mesh, n_corner) = build_structured_cube(n);
        let n_cells = mesh.cell_count();

        // Set up a trivial Stokes problem: uniform inlet velocity, zero-pressure outlet,
        // no-slip walls on all other boundary faces.
        let mut boundary_conditions = HashMap::new();
        let face_sets = cfd_3d::fem::AxialBoundaryClassifier::new(&mesh, n).classify();

        for &v in &face_sets.inlet_nodes {
            boundary_conditions.insert(
                v,
                cfd_core::physics::boundary::BoundaryCondition::VelocityInlet {
                    velocity: nalgebra::Vector3::new(0.0, 0.0, 0.01),
                },
            );
        }
        for &v in &face_sets.outlet_nodes {
            boundary_conditions.insert(
                v,
                cfd_core::physics::boundary::BoundaryCondition::PressureOutlet { pressure: 0.0 },
            );
        }
        for &v in &face_sets.wall_nodes {
            boundary_conditions.entry(v).or_insert(
                cfd_core::physics::boundary::BoundaryCondition::Dirichlet {
                    value: 0.0,
                    component_values: Some(vec![Some(0.0), Some(0.0), Some(0.0), None]),
                },
            );
        }

        let fluid = cfd_core::physics::fluid::ConstantPropertyFluid {
            name: "benchmark_water".to_string(),
            density: 998.0,
            viscosity: 1.002e-3,
            specific_heat: 4182.0,
            thermal_conductivity: 0.598,
            speed_of_sound: 1481.0,
        };

        let problem = cfd_3d::fem::StokesFlowProblem::new(
            mesh,
            fluid,
            boundary_conditions,
            n_corner,
        );

        group.bench_with_input(
            BenchmarkId::new("solve", n_cells),
            &problem,
            |b, problem| {
                b.iter(|| {
                    let config = cfd_3d::fem::FemConfig::<f64>::default();
                    let mut solver = cfd_3d::fem::FemSolver::new(config);
                    let result = solver.solve(problem, None);
                    black_box(result.is_ok());
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_extract_vertex_indices, bench_fem_solve_structured);
criterion_main!(benches);
