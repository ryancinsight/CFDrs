//! Tests for sparse matrix module

#[cfg(test)]
mod tests {
    use super::super::*;
    use approx::assert_relative_eq;
    use cfd_core::error::Result;
    #[test]
    fn test_matrix_builder() -> Result<()> {
        let mut builder = SparseMatrixBuilder::new(3, 3);
        builder.add_triplets(vec![(0, 0, 1.0), (0, 2, 2.0), (1, 1, 3.0)])?;
        let matrix = builder.build()?;
        // Check row 0 entries
        let row0 = matrix.row(0);
        let row0_entries: Vec<_> = row0
            .col_indices()
            .iter()
            .zip(row0.values().iter())
            .map(|(&col, &val)| (col, val))
            .collect();
        assert_eq!(row0_entries.len(), 2);
        assert!(row0_entries.contains(&(0, 1.0)));
        assert!(row0_entries.contains(&(2, 2.0)));
        // Check row 1 entries
        let row1 = matrix.row(1);
        let row1_entries: Vec<_> = row1
            .zip(row1.values().iter())
        assert_eq!(row1_entries.len(), 1);
        assert!(row1_entries.contains(&(1, 3.0)));
        Ok(())
    }
    fn test_tridiagonal_pattern() -> Result<()> {
        let matrix = SparsePatterns::tridiagonal(4, -1.0, 2.0, -1.0)?;
        assert_eq!(matrix.nrows(), 4);
        assert_eq!(matrix.ncols(), 4);
        assert_eq!(matrix.nnz(), 10); // 3*4 - 2 = 10
        // Check diagonal values
        let diag = matrix.diagonal();
        for i in 0..4 {
            assert_relative_eq!(diag[i], 2.0, epsilon = 1e-10);
        }
    fn test_five_point_stencil() -> Result<()> {
        let matrix = SparsePatterns::five_point_stencil(3, 3, 1.0, 1.0)?;
        // Check that it creates a 9x9 matrix (3x3 grid)
        assert_eq!(matrix.nrows(), 9);
        assert_eq!(matrix.ncols(), 9);
        // Center point should have 4 connections + diagonal
        let center = 4; // (1,1) in 3x3 grid
        let center_row = matrix.row(center);
        assert_eq!(center_row.nnz(), 5); // center + 4 neighbors
    }

    fn test_frobenius_norm() -> Result<()> {
        let matrix = SparsePatterns::tridiagonal(3, 1.0, 2.0, 1.0)?;
        let norm = matrix.frobenius_norm();
        // Expected: sqrt(2^2 * 3 + 1^2 * 4) = sqrt(12 + 4) = 4.0
        assert_relative_eq!(norm, 4.0, epsilon = 1e-10);
    }

    fn test_diagonal_dominance() -> Result<()> {
        // Create a diagonally dominant matrix
        let dominant = SparsePatterns::tridiagonal(3, 1.0, 4.0, 1.0)?;
        assert!(dominant.is_diagonally_dominant());
        // Create a non-diagonally dominant matrix
        let non_dominant = SparsePatterns::tridiagonal(3, 2.0, 3.0, 2.0)?;
        assert!(!non_dominant.is_diagonally_dominant());

    }


}
}
}
