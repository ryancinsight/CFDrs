use super::poisson::PoissonSolver;
use crate::grid::structured::StructuredGrid2D;
use std::collections::HashMap;

#[test]
fn test_poisson_mms_dirichlet_quadratic_interior_accuracy() {
    let nx = 64usize;
    let ny = 64usize;
    let grid = StructuredGrid2D::unit_square(nx, ny).expect("valid grid");

    // Manufactured solution: phi(x,y) = x^2 + y^2
    // ∇^2 phi = 2 + 2 = 4
    let mut source: HashMap<(usize, usize), f64> = HashMap::new();
    let mut boundary_values: HashMap<(usize, usize), f64> = HashMap::new();

    for (i, j) in grid.iter() {
        let center = grid.cell_center(i, j).unwrap();
        let x = center.x;
        let y = center.y;
        let phi = x * x + y * y;

        if grid.is_boundary(i, j) {
            boundary_values.insert((i, j), phi);
        } else {
            source.insert((i, j), 4.0f64);
        }
    }

    let solver = PoissonSolver::default();
    let result = solver.solve(&grid, &source, &boundary_values).expect("solve");

    // Check interior points only (boundary values enforced)
    let mut max_err = 0.0f64;
    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            let center = grid.cell_center(i, j).unwrap();
            let x = center.x;
            let y = center.y;
            let phi_exact = x * x + y * y;
            let phi_num = result.get(&(i, j)).copied().unwrap_or(0.0);
            let err = (phi_num - phi_exact).abs();
            if err > max_err { max_err = err; }
        }
    }

    assert!(max_err < 1.0e-3, "Interior MMS error too large: {}", max_err);
}
