#[cfg(test)]
mod tests {
    use cfd_math::linear_solver::preconditioners::multigrid::ruge_stueben_coarsening;
    use nalgebra::DMatrix;
    use nalgebra_sparse::CsrMatrix;

    fn create_disconnected_matrix() -> CsrMatrix<f64> {
        // Create a matrix with some disconnected components or weak connections
        // that might leave points undecided if the first pass isn't perfect.
        let n = 6;
        let mut matrix = DMatrix::zeros(n, n);

        // Block 1: 0-1-2 strongly connected
        matrix[(0, 0)] = 2.0;
        matrix[(0, 1)] = -1.0;
        matrix[(1, 0)] = -1.0;
        matrix[(1, 1)] = 2.0;
        matrix[(1, 2)] = -1.0;
        matrix[(2, 1)] = -1.0;
        matrix[(2, 2)] = 2.0;

        // Block 2: 3-4-5 strongly connected
        matrix[(3, 3)] = 2.0;
        matrix[(3, 4)] = -1.0;
        matrix[(4, 3)] = -1.0;
        matrix[(4, 4)] = 2.0;
        matrix[(4, 5)] = -1.0;
        matrix[(5, 4)] = -1.0;
        matrix[(5, 5)] = 2.0;

        // Weak connection between blocks
        matrix[(2, 3)] = -0.01;
        matrix[(3, 2)] = -0.01;

        CsrMatrix::from(&matrix)
    }

    fn create_matrix_with_island() -> CsrMatrix<f64> {
        // Create a matrix with an isolated point
        let n = 3;
        let mut matrix = DMatrix::zeros(n, n);

        matrix[(0, 0)] = 2.0;
        matrix[(0, 1)] = -1.0;
        matrix[(1, 0)] = -1.0;
        matrix[(1, 1)] = 2.0;

        // Point 2 is isolated (diagonal only)
        matrix[(2, 2)] = 1.0;

        CsrMatrix::from(&matrix)
    }

    #[test]
    fn test_undecided_points() {
        let matrix = create_disconnected_matrix();
        // Use a threshold that ignores the weak connection
        let result = ruge_stueben_coarsening(&matrix, 0.25).unwrap();

        // Check if all points are decided (either in coarse_points or mapped in fine_to_coarse_map)
        // Note: fine_to_coarse_map maps F points to C points.
        // C points map to themselves (or rather, their index in coarse_points).
        // Let's verify status coverage.

        let n = matrix.nrows();
        // In the current implementation, C-points are in `coarse_points`.
        // F-points are those in `fine_to_coarse_map` that are NOT C-points but have a mapping?
        // Actually, `fine_to_coarse_map` is populated for C-points too:
        // `fine_to_coarse_map[i] = Some(coarse_points.len() - 1);` when i is made C-point.
        // And for F-points: `fine_to_coarse_map[i] = Some(idx);` in Step 4.

        // So checking if `fine_to_coarse_map[i]` is Some should be enough to say it was handled...
        // UNLESS Step 4 fails to find a C-point for an F-point.

        // But the issue is "undecided points". If they are undecided after Step 3, they enter Step 4 with `fine_to_coarse_map[i] == None`.
        // Step 4 tries to map them to a neighbor C-point.
        // If an undecided point becomes F-point implicitly but has no strong C-neighbor, it might remain None?
        // Or does "Undecided" mean `status == 0`?

        // In the code:
        // Step 4 iterates 0..n. `if fine_to_coarse_map[i].is_none()`.
        // It tries to find a neighbor C-point.
        // If it finds one, it maps.
        // If it doesn't find one, it remains None.

        // So we should check if any point has `fine_to_coarse_map` as None.

        let mut unmapped_count = 0;
        for i in 0..n {
            if result.fine_to_coarse_map[i].is_none() {
                unmapped_count += 1;
                println!("Point {} is unmapped", i);
            }
        }

        assert_eq!(
            unmapped_count, 0,
            "Found {} unmapped points. Second pass missing?",
            unmapped_count
        );
    }

    #[test]
    fn test_island_points() {
        let matrix = create_matrix_with_island();
        let result = ruge_stueben_coarsening(&matrix, 0.25).unwrap();

        let n = matrix.nrows();
        let mut unmapped_count = 0;
        for i in 0..n {
            if result.fine_to_coarse_map[i].is_none() {
                unmapped_count += 1;
                println!("Island point {} is unmapped", i);
            }
        }

        assert_eq!(
            unmapped_count, 0,
            "Found {} unmapped island points. Second pass missing?",
            unmapped_count
        );
    }
}
