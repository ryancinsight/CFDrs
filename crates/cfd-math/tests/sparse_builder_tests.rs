use cfd_math::sparse::SparseMatrixBuilder;

#[test]
fn test_sparse_builder_summation() {
    let rows = 2;
    let cols = 2;

    // Case 1: allow_duplicates = false (Current default, uses HashMap)
    let mut builder1 = SparseMatrixBuilder::<f64>::new(rows, cols).allow_duplicates(false);
    builder1.add_entry(0, 0, 1.0).unwrap();
    builder1.add_entry(0, 0, 2.0).unwrap(); // Duplicate
    builder1.add_entry(1, 1, 3.0).unwrap();
    let csr1 = builder1.build().unwrap();

    // Case 2: allow_duplicates = true (Optimization candidate, skips HashMap)
    let mut builder2 = SparseMatrixBuilder::<f64>::new(rows, cols).allow_duplicates(true);
    builder2.add_entry(0, 0, 1.0).unwrap();
    builder2.add_entry(0, 0, 2.0).unwrap(); // Duplicate
    builder2.add_entry(1, 1, 3.0).unwrap();
    let csr2 = builder2.build().unwrap();

    // Check if results are identical
    assert_eq!(csr1.nnz(), csr2.nnz(), "Number of non-zeros should match");
    assert_eq!(csr1.row_offsets(), csr2.row_offsets(), "Row offsets should match");
    assert_eq!(csr1.col_indices(), csr2.col_indices(), "Column indices should match");
    assert_eq!(csr1.values(), csr2.values(), "Values should match");

    // Specifically check the summed value
    assert_eq!(csr1.values()[0], 3.0);
    assert_eq!(csr2.values()[0], 3.0);
}
