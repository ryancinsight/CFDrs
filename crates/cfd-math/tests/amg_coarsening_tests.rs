#[cfg(test)]
mod tests {
    use cfd_math::linear_solver::preconditioners::multigrid::ruge_stueben_coarsening;
    use leto_ops::CsrMatrix;

    fn create_disconnected_matrix() -> CsrMatrix<f64> {
        CsrMatrix::from_parts(
            vec![
                2.0, -1.0, -1.0, 2.0, -1.0, -1.0, 2.0, -0.01, -0.01, 2.0, -1.0, -1.0, 2.0, -1.0,
                -1.0, 2.0,
            ],
            vec![0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5],
            vec![0, 2, 5, 8, 11, 14, 16],
            6,
            6,
        )
        .unwrap()
    }

    fn create_matrix_with_island() -> CsrMatrix<f64> {
        CsrMatrix::from_parts(
            vec![2.0, -1.0, -1.0, 2.0, 1.0],
            vec![0, 1, 0, 1, 2],
            vec![0, 2, 4, 5],
            3,
            3,
        )
        .unwrap()
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
            }
        }

        assert_eq!(
            unmapped_count, 0,
            "Found {} unmapped island points. Second pass missing?",
            unmapped_count
        );
    }
}
