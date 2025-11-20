use super::advection_diffusion::AdvectionDiffusionSolver;
use super::config::FdmConfig;
use crate::grid::structured::StructuredGrid2D;
use std::collections::HashMap;

#[test]
fn test_advection_diffusion_mms_dirichlet_linear_solution() {
    let nx = 32usize;
    let ny = 32usize;
    let grid = StructuredGrid2D::unit_square(nx, ny).expect("grid");
    let dx = 1.0f64 / (nx as f64 - 1.0);
    let dy = 1.0f64 / (ny as f64 - 1.0);

    let mut vel_x: HashMap<(usize, usize), f64> = HashMap::new();
    let mut vel_y: HashMap<(usize, usize), f64> = HashMap::new();
    for (i, j) in grid.iter() {
        vel_x.insert((i, j), 0.7);
        vel_y.insert((i, j), 0.3);
    }

    let diffusivity = 0.1f64;
    let mut source: HashMap<(usize, usize), f64> = HashMap::new();
    let mut boundary: HashMap<(usize, usize), f64> = HashMap::new();

    for (i, j) in grid.iter() {
        let x = i as f64 * dx;
        let y = j as f64 * dy;
        let phi = x + y;
        if grid.is_boundary(i, j) {
            boundary.insert((i, j), phi);
        } else {
            source.insert((i, j), 0.7 + 0.3);
        }
    }

    let cfg = FdmConfig::default();
    let solver = AdvectionDiffusionSolver::new(cfg);
    let result = solver
        .solve_steady(&grid, &vel_x, &vel_y, diffusivity, &source, &boundary)
        .expect("solve");

    let mut max_err = 0.0f64;
    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            let x = i as f64 * dx;
            let y = j as f64 * dy;
            let phi_exact = x + y;
            let phi_num = *result.get(&(i, j)).unwrap();
            let err = (phi_num - phi_exact).abs();
            if err > max_err { max_err = err; }
        }
    }
    assert!(max_err < 5.0e-3, "MMS error too large: {}", max_err);
}