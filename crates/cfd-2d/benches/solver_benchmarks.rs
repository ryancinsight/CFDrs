use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use cfd_2d::{
    grid::{StructuredGrid, BoundaryCondition},
    fdm::{PoissonSolver, AdvectionDiffusionSolver},
    fvm::FvmSolver,
    lbm::{LbmSolver, D2Q9},
    pressure_velocity_coupling::{PressureVelocityCouplerSolver, PressureVelocityCouplingConfig},
};
use nalgebra::Vector2;

fn benchmark_fdm_solvers(c: &mut Criterion) {
    let mut group = c.benchmark_group("fdm_solvers");
    
    for size in [50, 100, 200].iter() {
        let grid = create_test_grid(*size, *size);
        let mut poisson_solver = PoissonSolver::new(grid.clone());
        let mut advection_diffusion_solver = AdvectionDiffusionSolver::new(grid.clone());
        
        // Setup initial conditions
        let initial_field: Vec<Vec<f64>> = (0..*size)
            .map(|i| (0..*size).map(|j| ((i + j) as f64).sin()).collect())
            .collect();
        
        let velocity_field: Vec<Vec<Vector2<f64>>> = (0..*size)
            .map(|i| (0..*size).map(|j| Vector2::new(1.0, 0.5)).collect())
            .collect();
        
        group.bench_with_input(
            BenchmarkId::new("poisson_solver", size),
            size,
            |b, _| {
                b.iter(|| {
                    let mut field = initial_field.clone();
                    black_box(poisson_solver.solve(&mut field, 0.01).unwrap())
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("advection_diffusion_solver", size),
            size,
            |b, _| {
                b.iter(|| {
                    let mut field = initial_field.clone();
                    black_box(advection_diffusion_solver.solve(&mut field, &velocity_field, 0.01, 0.001).unwrap())
                })
            },
        );
    }
    
    group.finish();
}

fn benchmark_fvm_solver(c: &mut Criterion) {
    let mut group = c.benchmark_group("fvm_solver");
    
    for size in [50, 100, 200].iter() {
        let grid = create_test_grid(*size, *size);
        let mut fvm_solver = FvmSolver::new(grid);
        
        let initial_field: Vec<Vec<f64>> = (0..*size)
            .map(|i| (0..*size).map(|j| ((i + j) as f64).sin()).collect())
            .collect();
        
        group.bench_with_input(
            BenchmarkId::new("fvm_diffusion_solve", size),
            size,
            |b, _| {
                b.iter(|| {
                    let mut field = initial_field.clone();
                    black_box(fvm_solver.solve_diffusion(&mut field, 0.01, 0.001).unwrap())
                })
            },
        );
    }
    
    group.finish();
}

fn benchmark_lbm_solver(c: &mut Criterion) {
    let mut group = c.benchmark_group("lbm_solver");
    
    for size in [50, 100, 200].iter() {
        let mut lbm_solver = LbmSolver::<f64, D2Q9>::new(*size, *size, 0.1);
        
        group.bench_with_input(
            BenchmarkId::new("lbm_single_step", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(lbm_solver.step().unwrap())
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("lbm_equilibrium_calculation", size),
            size,
            |b, _| {
                b.iter(|| {
                    for i in 0..*size {
                        for j in 0..*size {
                            let density = 1.0;
                            let velocity = Vector2::new(0.1, 0.05);
                            for q in 0..9 {
                                black_box(lbm_solver.equilibrium_distribution(q, &density, &velocity));
                            }
                        }
                    }
                })
            },
        );
    }
    
    group.finish();
}

fn benchmark_simple_solver(c: &mut Criterion) {
    let mut group = c.benchmark_group("simple_solver");
    
    for size in [25, 50, 100].iter() {
        let config = SimpleConfig::default();
        let mut simple_solver = SimpleSolver::new(config, *size, *size);
        
        group.bench_with_input(
            BenchmarkId::new("simple_single_iteration", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(simple_solver.solve_step().unwrap())
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("simple_pressure_correction", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(simple_solver.solve_pressure_correction().unwrap())
                })
            },
        );
    }
    
    group.finish();
}

fn benchmark_grid_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("grid_operations");
    
    for size in [100, 500, 1000].iter() {
        let grid = create_test_grid(*size, *size);
        
        group.bench_with_input(
            BenchmarkId::new("grid_creation", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(create_test_grid(*size, *size))
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("neighbor_iteration", size),
            size,
            |b, _| {
                b.iter(|| {
                    let mut count = 0;
                    for i in 0..*size {
                        for j in 0..*size {
                            for neighbor in grid.neighbors(i, j) {
                                count += neighbor.0 + neighbor.1;
                            }
                        }
                    }
                    black_box(count)
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("boundary_detection", size),
            size,
            |b, _| {
                b.iter(|| {
                    let mut boundary_count = 0;
                    for i in 0..*size {
                        for j in 0..*size {
                            if grid.is_boundary(i, j) {
                                boundary_count += 1;
                            }
                        }
                    }
                    black_box(boundary_count)
                })
            },
        );
    }
    
    group.finish();
}

fn benchmark_memory_access_patterns(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_access");
    
    for size in [100, 500, 1000].iter() {
        let field: Vec<Vec<f64>> = (0..*size)
            .map(|i| (0..*size).map(|j| (i + j) as f64).collect())
            .collect();
        
        group.bench_with_input(
            BenchmarkId::new("row_major_access", size),
            size,
            |b, _| {
                b.iter(|| {
                    let mut sum = 0.0;
                    for i in 0..*size {
                        for j in 0..*size {
                            sum += field[i][j];
                        }
                    }
                    black_box(sum)
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("column_major_access", size),
            size,
            |b, _| {
                b.iter(|| {
                    let mut sum = 0.0;
                    for j in 0..*size {
                        for i in 0..*size {
                            sum += field[i][j];
                        }
                    }
                    black_box(sum)
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("stencil_operations", size),
            size,
            |b, _| {
                b.iter(|| {
                    let mut result = vec![vec![0.0; *size]; *size];
                    for i in 1..(*size-1) {
                        for j in 1..(*size-1) {
                            result[i][j] = 0.25 * (
                                field[i-1][j] + field[i+1][j] + 
                                field[i][j-1] + field[i][j+1]
                            );
                        }
                    }
                    black_box(result)
                })
            },
        );
    }
    
    group.finish();
}

// Helper functions
fn create_test_grid(nx: usize, ny: usize) -> StructuredGrid<f64> {
    let mut grid = StructuredGrid::new(nx, ny, 1.0 / nx as f64, 1.0 / ny as f64);
    
    // Set boundary conditions
    for i in 0..nx {
        grid.set_boundary_condition(i, 0, BoundaryCondition::Dirichlet(0.0));
        grid.set_boundary_condition(i, ny-1, BoundaryCondition::Dirichlet(0.0));
    }
    for j in 0..ny {
        grid.set_boundary_condition(0, j, BoundaryCondition::Dirichlet(0.0));
        grid.set_boundary_condition(nx-1, j, BoundaryCondition::Dirichlet(1.0));
    }
    
    grid
}

criterion_group!(
    benches,
    benchmark_fdm_solvers,
    benchmark_fvm_solver,
    benchmark_lbm_solver,
    benchmark_simple_solver,
    benchmark_grid_operations,
    benchmark_memory_access_patterns
);
criterion_main!(benches);
