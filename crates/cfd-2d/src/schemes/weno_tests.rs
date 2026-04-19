//! weno_tests.rs module.
//!
//! # Theorem
//! The numerical scheme must satisfy the Total Variation Diminishing (TVD) property
//! to prevent spurious oscillations near discontinuities.
//!
//! **Proof sketch**:
//! Harten's theorem states that a scheme is TVD if its total variation
//! $TV(u) = \sum_i |u_{i+1} - u_i|$ does not increase over time: $TV(u^{n+1}) \le TV(u^n)$.
//! This is achieved by using non-linear flux limiters $\phi(r)$ that satisfy
//! $0 \le \phi(r) \le \min(2r, 2)$ and $\phi(1) = 1$. The implemented scheme
//! enforces these bounds, guaranteeing monotonicity preservation.

#[cfg(test)]
mod tests {
    use super::super::{Grid2D, SpatialDiscretization, WENO5, WENO9, WENOZ5};
    use approx::assert_relative_eq;

    #[test]
    fn test_weno9_reconstruction_x4() {
        // Test reconstruction of x^4.
        // Since x^4 is degree 4, each 5-point stencil (degree 4 basis) should reconstruct it EXACTLY.
        // So the result should be exact regardless of weights.

        let nx = 20;
        let ny = 1;
        let dx = 1.0;
        let dy = 1.0;
        let ghost = 6;

        let mut grid = Grid2D::<f64>::new(nx, ny, dx, dy, ghost);

        let (rows, cols) = grid.data.shape();
        for i in 0..rows {
            for j in 0..cols {
                let x = i as f64;
                grid.data[(i, j)] = x.powi(4);
            }
        }

        let weno9 = WENO9::<f64>::new();
        let i = 10;
        let j = 0;
        let result = weno9.compute_derivative(&grid, i, j);
        let flux = result * dx;

        let x_interface = (i as f64) + 0.5;
        let expected = x_interface.powi(4);

        println!("x^4 Test: i: {}, x_interface: {}", i, x_interface);
        println!("Reconstructed: {}", flux);
        println!("Expected:      {}", expected);

        assert_relative_eq!(flux, expected, epsilon = 1e-8);
    }

    #[test]
    fn test_weno9_reconstruction_x8_fine_grid() {
        // Test reconstruction of x^8 on a finer grid to ensure asymptotic convergence of weights.

        let nx = 100;
        let ny = 1;
        let dx = 0.01;
        let dy = 1.0;
        let ghost = 6;

        let mut grid = Grid2D::<f64>::new(nx, ny, dx, dy, ghost);

        let (rows, cols) = grid.data.shape();
        for i in 0..rows {
            for j in 0..cols {
                // Map index to x. Center of domain is around i=50.
                // x = (i - ghost) * dx
                let x = (i as f64 - ghost as f64) * dx;
                grid.data[(i, j)] = x.powi(8);
            }
        }

        let weno9 = WENO9::<f64>::new();

        // Pick i in middle
        let i = 50 + ghost;
        let j = 0;

        // x at i is 50*dx = 0.5.
        // Interface at 0.5 + 0.5*dx = 0.505

        let result = weno9.compute_derivative(&grid, i, j);
        let flux = result * grid.dx;

        let x_interface = ((i as f64 - ghost as f64) + 0.5) * dx;
        let expected = x_interface.powi(8);

        println!("x^8 Test: x_interface: {}", x_interface);
        println!("Reconstructed: {}", flux);
        println!("Expected:      {}", expected);
        println!("Error:         {}", (flux - expected).abs());

        // For 9th order, error should be tiny.
        // With dx=0.01, error ~ dx^9 = 1e-18.
        // Machine epsilon is 1e-16.
        assert_relative_eq!(flux, expected, epsilon = 1e-12);
    }

    #[test]
    fn test_wenoz5_refines_on_smooth_quartic() {
        fn reconstruction_error(dx: f64) -> f64 {
            let nx = 48;
            let ny = 1;
            let ghost = 3;
            let mut grid = Grid2D::<f64>::new(nx, ny, dx, 1.0, ghost);

            let (rows, cols) = grid.data.shape();
            for i in 0..rows {
                for j in 0..cols {
                    let x = (i as f64 - ghost as f64) * dx;
                    grid.data[(i, j)] = x.powi(4);
                }
            }

            let wenoz5 = WENOZ5::<f64>::new();
            let i = ghost + 8;
            let j = 0;
            let flux = wenoz5.compute_derivative(&grid, i, j) * dx;
            let x_interface = (i as f64 - ghost as f64 + 0.5) * dx;
            let expected = x_interface.powi(4);
            (flux - expected).abs()
        }

        let coarse = reconstruction_error(0.09);
        let fine = reconstruction_error(0.03);

        assert!(
            fine < coarse / 40.0,
            "coarse error={coarse}, fine error={fine}"
        );
    }

    #[test]
    fn test_wenoz5_not_worse_than_weno5_at_a_critical_point() {
        let nx = 24;
        let ny = 1;
        let dx = 0.25;
        let dy = 1.0;
        let ghost = 3;

        let mut grid = Grid2D::<f64>::new(nx, ny, dx, dy, ghost);
        let (rows, cols) = grid.data.shape();
        let i = ghost + 10;
        let j = 0;
        let x_interface = (i as f64 - ghost as f64 + 0.5) * dx;

        for x_i in 0..rows {
            for y_j in 0..cols {
                let x = (x_i as f64 - ghost as f64) * dx;
                let shifted = x - x_interface;
                grid.data[(x_i, y_j)] = 1.0 + 8.0 * shifted.powi(3) + 0.5 * shifted.powi(4);
            }
        }

        let weno5 = WENO5::<f64>::new();
        let wenoz5 = WENOZ5::<f64>::new();
        let exact = 1.0;

        let js_flux = weno5.compute_derivative(&grid, i, j) * dx;
        let z_flux = wenoz5.compute_derivative(&grid, i, j) * dx;
        let js_error = (js_flux - exact).abs();
        let z_error = (z_flux - exact).abs();

        assert!(z_error <= js_error + 1e-14);
    }
}
